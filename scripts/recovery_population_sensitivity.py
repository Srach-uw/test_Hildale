"""Compare inference choices on fixed local recovery membership, not a release fit."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from extract_eccentricity_posteriors_direct import (
    _sample_frame, match_planets_by_period, paired_period_samples, read_alderaan_planets,
)
from hierarchical_rayleigh import fit_from_mass_matrix
from real_duration_measure_audit import compare, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--fits-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--draws', type=int, default=50000)
    parser.add_argument('--seeds', type=int, default=2)
    parser.add_argument('--disk', choices=['thin', 'thick'])
    parser.add_argument('--system', choices=['single', 'multi'])
    parser.add_argument('--planet', help='Restrict extraction after full-system period matching')
    parser.add_argument('--proposal-mode', choices=['random', 'sobol'], default='random')
    parser.add_argument('--source-settings', action='store_true',
                        help='Require and replay stored proposal support and density-error mode')
    args = parser.parse_args()
    catalog = pd.read_csv(args.summary)
    if catalog.kepoi_name.duplicated().any():
        raise ValueError('Duplicate catalog planets')
    flags = catalog.qc_primary_exclude.astype(str).str.lower()
    if not flags.isin(['true', 'false']).all():
        raise ValueError('Missing or ambiguous QC flags')
    if args.disk:
        catalog = catalog.loc[catalog.disk == args.disk]
    if args.system:
        catalog = catalog.loc[catalog.system == args.system]
    if args.planet:
        target_rows = catalog.loc[catalog.kepoi_name == args.planet, 'koi_target']
        if len(target_rows) != 1:
            raise ValueError('Requested planet must occur exactly once')
        catalog = catalog.loc[catalog.koi_target == target_rows.iloc[0]]
    paths = {}
    for path in args.fits_root.rglob('*-results.fits'):
        target = path.name.removesuffix('-results.fits')
        if target in paths:
            raise ValueError(f'Duplicate FITS target: {target}')
        paths[target] = path
    args.output.mkdir(parents=True, exist_ok=True)
    egrid = np.linspace(0, .95, 300)
    sigmas = np.linspace(.002, 1, 600)
    records, ledger, matrices = [], [], {}
    for number, (target, group) in enumerate(catalog.groupby('koi_target')):
        group = group.reset_index(drop=True)
        if target not in paths:
            for item in group.itertuples():
                ledger.append(dict(planet=item.kepoi_name, target=target, status='no_recovery_FITS'))
            continue
        path = paths[target]
        digest = sha256(path)
        with fits.open(path) as hdul:
            samples = _sample_frame(hdul['SAMPLES'].data)
            matches = match_planets_by_period(group, read_alderaan_planets(hdul), .01)
            mapping = {i: index for i, index, _, _ in matches}
            for i, item in group.iterrows():
                if args.planet and item.kepoi_name != args.planet:
                    continue
                entry = dict(planet=item.kepoi_name, target=target, fits_sha256=digest)
                if str(item.qc_primary_exclude).lower() == 'true':
                    ledger.append(dict(**entry, status='recorded_QC_exclusion'))
                    continue
                if i not in mapping:
                    ledger.append(dict(**entry, status='period_unmatched'))
                    continue
                rho, hi, lo = (float(item[k]) for k in ['rho_true_solar', 'rho_err_hi_solar', 'rho_err_lo_solar'])
                if not np.all(np.isfinite([rho, hi, lo])) or min(rho, hi, lo) <= 0:
                    ledger.append(dict(**entry, status='invalid_density'))
                    continue
                index = mapping[i]
                periods = paired_period_samples(samples, hdul, index)
                emax = float(item.proposal_e_max) if args.source_settings else .95
                density_mode = str(item.density_error_mode) if args.source_settings else 'symmetric-average'
                if args.source_settings and (item.nested_weight_mode != 'dynesty'
                        or item.period_sampling_mode != 'paired_alderaan'
                        or item.density_sampling_mode != 'fixed_central'
                        or item.impact_mode != 'alderaan' or float(item.e_max) != .95):
                    raise ValueError('Unsupported source settings; refusing implicit conversion')
                for seed in range(args.seeds):
                    results = compare(samples, periods, index, rho, hi, lo, seed, args.draws, egrid,
                                      args.proposal_mode, emax, density_mode)
                    for result in results:
                        masses = result.pop('_masses')
                        records.append(dict(planet=item.kepoi_name, disk=item.disk, system=item.system, **result))
                        for selection, mass in masses.items():
                            key = (seed, item.disk, item.system, result['mode'], selection)
                            matrices.setdefault(key, []).append((item.kepoi_name, mass))
                ledger.append(dict(**entry, status='included'))
        if len(ledger) % 20 < 4:
            print(f'Processed {target}; {len(records)} posterior summaries', flush=True)
    pd.DataFrame(records).to_csv(args.output/'planet_diagnostics.csv', index=False)
    pd.DataFrame(ledger).to_csv(args.output/'membership.csv', index=False)
    results = []
    for key, members in matrices.items():
        seed, disk, system, mode, selection = key
        matrix = np.array([mass for _, mass in members])
        mass_dir = args.output/'mass_matrices'
        mass_dir.mkdir(exist_ok=True)
        np.savez_compressed(mass_dir/f'{seed}_{disk}_{system}_{mode}_{selection}.npz',
                            mass_matrix=matrix, planets=np.array([p for p, _ in members]),
                            e_grid=egrid, sigma_grid=sigmas)
        fit = fit_from_mass_matrix(matrix, egrid, sigmas,
                                   selection != 'none', selection)
        results.append(dict(seed=seed, disk=disk, system=system, mode=mode,
                            n=len(members), membership='|'.join(p for p, _ in members), **fit))
    frame = pd.DataFrame(results)
    for _, group in frame.groupby(['seed', 'disk', 'system']):
        if group.membership.nunique() != 1:
            raise ValueError('Sensitivity membership changed')
    frame.to_csv(args.output/'population_diagnostics.csv', index=False)
    (args.output/'scope.json').write_text(json.dumps(dict(
        summary_sha256=sha256(args.summary), draws=args.draws, seeds=args.seeds,
        disk_filter=args.disk, system_filter=args.system, planet_filter=args.planet,
        proposal_mode=args.proposal_mode,
        source_settings=args.source_settings,
        e_grid_nodes=len(egrid), e_max=.95, sigma_min=.002, sigma_max=1.,
        membership='Recovery FITS overlap, summary QC flag false, period match, positive density errors',
        limits='Diagnostic subset, not canonical release QC or the paper sample. ESS and seed sensitivity must be inspected. No detection-efficiency model. Treat boundary fits as limits.'
    ), indent=2)+'\n')
    print(f'Finished {len(results)} fixed-membership population diagnostics', flush=True)


if __name__ == '__main__':
    main()
