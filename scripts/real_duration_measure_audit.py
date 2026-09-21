"""Replay paired transit samples under two explicit density measures.

This sensitivity audit leaves canonical posteriors untouched. The transformed
measure assumes a log-uniform duration proposal and a density prior independent
of eccentricity. It is not asserted to be Sagear's implementation.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.stats import qmc

from duration_measure_control import density_and_derivative
from extract_eccentricity_posteriors_direct import (
    _sample_frame, density_log_likelihood, match_planets_by_period,
    nested_sample_weights, paired_period_samples, read_alderaan_planets,
    weighted_quantile,
)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(samples, period, index, rho, hi, lo, seed, draws, e_grid=None, proposal_mode='random',
            e_max=.95, density_error_mode='symmetric-average'):
    if not 0 < e_max < 1:
        raise ValueError('Invalid eccentricity proposal support')
    rng = np.random.default_rng(seed)
    weights = nested_sample_weights(samples['LN_WT'].to_numpy(float), 'dynesty')
    if proposal_mode == 'random':
        selected = rng.choice(len(samples), size=draws, p=weights)
        e = rng.uniform(0, e_max, draws)
        w = rng.uniform(-np.pi/2, 3*np.pi/2, draws)
    elif proposal_mode == 'sobol':
        if draws < 2 or draws & (draws-1):
            raise ValueError('Sobol draw count must be a power of two')
        unit = qmc.Sobol(d=3, scramble=True, seed=seed).random_base2(int(np.log2(draws)))
        cumulative = np.cumsum(weights)
        cumulative[-1] = 1.
        selected = np.searchsorted(cumulative, unit[:, 0], side='right')
        e = e_max*unit[:, 1]
        w = -np.pi/2+2*np.pi*unit[:, 2]
    else:
        raise ValueError('Unknown proposal mode')
    t = samples[f'DUR14_{index}'].to_numpy(float)[selected]
    r = samples[f'ROR_{index}'].to_numpy(float)[selected]
    b = samples[f'IMPACT_{index}'].to_numpy(float)[selected]
    p = period[selected]
    angle = np.pi*t/p*(1+e*np.sin(w))/np.sqrt(1-e*e)
    with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
        implied, derivative = density_and_derivative(t, e, w, p, r, b)
        loglike = density_log_likelihood(implied, rho, hi, lo, density_error_mode)
    valid = ((angle > 0) & (angle < np.pi/2) & (b >= 0) & (b < 1+r)
             & np.isfinite(loglike) & np.isfinite(derivative) & (derivative > 0)
             & (t > 0) & (implied > 0))
    if not valid.any():
        raise ValueError('No physical duration proposals')
    rows = []
    modes = ['density_likelihood_only', 'log_duration_measure', 'uniform_focal_radius_only']
    if e_grid is not None:
        modes.append('uniform_all_radii')
    for mode in modes:
        logweight = loglike[valid].copy()
        if mode == 'log_duration_measure':
            logweight += np.log(t[valid]) + np.log(derivative[valid])
        elif mode == 'uniform_focal_radius_only':
            # Equal support gives p_uniform(r) / p_loguniform(r) proportional to r.
            logweight += np.log(r[valid])
        elif mode == 'uniform_all_radii':
            for column in samples:
                if column.startswith('ROR_'):
                    radius = samples[column].to_numpy(float)[selected][valid]
                    if not np.all(np.isfinite(radius) & (radius >= 1e-5) & (radius <= .99)):
                        raise ValueError('Invalid joint radius-prior support')
                    logweight += np.log(radius)
        iw = np.exp(logweight-logweight.max())
        iw /= iw.sum()
        q = weighted_quantile(e[valid], iw, [.16, .5, .84])
        rows.append(dict(mode=mode, seed=seed, e16=q[0], e50=q[1], e84=q[2],
                         mean_e=float(iw @ e[valid]), ess=float(1/(iw @ iw)),
                         valid_fraction=float(valid.mean())))
        if e_grid is not None:
            edges = np.r_[-np.inf, (e_grid[:-1]+e_grid[1:])/2, np.inf]
            transit = (1+e[valid]*np.sin(w[valid]))/(1-e[valid]**2)
            rows[-1]['_masses'] = {
                selection: np.histogram(e[valid], bins=edges, weights=iw*factor)[0]
                for selection, factor in [('none', np.ones(len(iw))),
                    ('manuscript_reciprocal', 1/transit), ('legacy_forward_norm', transit)]}
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--fits-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--draws', type=int, default=250000)
    parser.add_argument('--seeds', type=int, default=3)
    args = parser.parse_args()
    catalog = pd.read_csv(args.summary)
    if catalog['kepoi_name'].duplicated().any():
        raise ValueError('Duplicate planet identifiers')
    rows, manifest = [], []
    for path in sorted(args.fits_dir.glob('*-results.fits')):
        target = path.name.removesuffix('-results.fits')
        group = catalog.loc[catalog.koi_target == target].reset_index(drop=True)
        if group.empty:
            manifest.append(dict(target=target, planet='', status='absent_from_density_summary',
                                 fits_sha256=sha256(path)))
            continue
        with fits.open(path) as hdul:
            samples = _sample_frame(hdul['SAMPLES'].data)
            matches = match_planets_by_period(group, read_alderaan_planets(hdul), .01)
            matched = {i for i, *_ in matches}
            for i, item in group.iterrows():
                if i not in matched:
                    manifest.append(dict(target=target, planet=item.kepoi_name, status='unmatched_period'))
            for i, index, _, mismatch in matches:
                item = group.iloc[i]
                rho, hi, lo = (float(item[k]) for k in ['rho_true_solar', 'rho_err_hi_solar', 'rho_err_lo_solar'])
                if not np.all(np.isfinite([rho, hi, lo])) or min(rho, hi, lo) <= 0:
                    manifest.append(dict(target=target, planet=item.kepoi_name, status='invalid_density'))
                    continue
                periods = paired_period_samples(samples, hdul, index)
                for seed in range(args.seeds):
                    for result in compare(samples, periods, index, rho, hi, lo, seed, args.draws):
                        rows.append(dict(target=target, planet=item.kepoi_name,
                                         disk=item.disk, system=item.system,
                                         density_fractional_error=(hi+lo)/(2*rho), **result))
                manifest.append(dict(target=target, planet=item.kepoi_name, status='compared',
                                     fits_sha256=sha256(path), period_relative_error=mismatch))
        print(f'{target}: {len(matches)} matched planets', flush=True)
    args.output.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(args.output/'comparison.csv', index=False)
    averages = frame.groupby(['planet', 'mode']).agg(
        e50_mean=('e50', 'mean'), e50_seed_std=('e50', 'std'),
        minimum_ess=('ess', 'min'), density_fractional_error=('density_fractional_error', 'first'))
    averages.to_csv(args.output/'planet_summary.csv')
    pd.DataFrame(manifest).to_csv(args.output/'manifest.csv', index=False)
    (args.output/'scope.json').write_text(json.dumps(dict(
        summary_sha256=sha256(args.summary), draws=args.draws, seeds=args.seeds,
        e_max=.95, density_error_mode='symmetric-average',
        transit_selection='not applied; per-planet sensitivity only',
        limitation='Selected validation systems; not a representative population fit. All modes enforce physical branch. Prior bounds remain those of the input fits. Uniform focal-radius reweighting requires adequate original posterior support and leaves companion-radius and timing priors unchanged.'
    ), indent=2)+'\n')


if __name__ == '__main__':
    main()
