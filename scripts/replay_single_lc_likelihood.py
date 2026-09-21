"""Replay the pinned ALDERAAN single-planet LC likelihood from saved inputs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from astropy.io import fits


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['directory', 'ephemeris-source', 'dependency-path', 'output']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--target', default='K02712')
    p.add_argument('--profile', action='store_true')
    p.add_argument('--leave-out', action='store_true', help='Refit low/high impact after omitting quarters or seasons')
    p.add_argument('--residual-report', action='store_true', help='Report best stored-sample photometric residual diagnostics after replay passes')
    a = p.parse_args()
    sys.path.insert(0, str(a.dependency_path))
    from batman import _quadratic_ld
    spec = importlib.util.spec_from_file_location('replay_ephemeris', a.ephemeris_source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    base = a.directory / a.target
    paths = [Path(str(base)+s) for s in ['_transit_parameters.csv', '_lc_filtered.fits', '_00_quick.ttvs', '-results.fits']]
    cat = pd.read_csv(paths[0])
    if len(cat) != 1 or int(cat.iloc[0].npl) != 1:
        raise ValueError('Only single-planet LC replay is supported')
    if Path(str(base)+'_sc_filtered.fits').exists():
        raise ValueError('SC input requires a separate replay implementation')
    pars = cat.iloc[0]
    with fits.open(paths[1]) as h:
        time, flux, error, quarter = [np.asarray(h[k].data).ravel() for k in ['TIME', 'FLUX', 'ERROR', 'QUARTER']]
    if not np.isfinite(np.column_stack([time, flux, error])).all() or (error <= 0).any():
        raise ValueError('Invalid photometry; refusing an undocumented mask')
    tt = np.loadtxt(paths[2], ndmin=2)
    ephem = module.Ephemeris(tt[:, 0].astype(int), tt[:, 2])
    exposure = 29.4243885/60/24
    duration = pars.duration/24
    depth = pars.depth*1e-6
    tau = 13*(pars.period/365.25)**(1/3)*np.sqrt(depth)/24
    sigma = np.mean(error/flux)*.04
    n = int(np.ceil(np.sqrt((depth/tau)*(exposure/8/sigma))))
    oversample = max(n+(n % 2+1), 7)
    mask = np.zeros(len(time), dtype=bool)
    for center in tt[:, 2]:
        if time.min() <= center <= time.max():
            mask |= abs(time-center) < max(1/24, 1.5*duration)
    groups = []
    for q in np.unique(quarter):
        use = mask & (quarter == q)
        if use.sum() <= max(1, int(np.floor(duration/exposure))):
            continue
        warped, index = ephem._warp_times(time[use], return_inds=True)
        legx = (index-int(tt[-1, 0])//2)/(tt[-1, 0]/2)
        groups.append((int(q), time[use], warped, legx, flux[use], error[use]))
    with fits.open(paths[3]) as h:
        samples = h['SAMPLES'].data.copy()
    chosen = np.unique(np.concatenate([np.linspace(0, len(samples)-1, 12, dtype=int),
                                       np.argsort(samples['LN_WT'])[-12:]]))
    offsets = np.linspace(-exposure/2, exposure/2, oversample)
    def residual(s, excluded=(), return_by_quarter=False):
        q1, q2 = s['LD_Q1'], s['LD_Q2']
        u1, u2 = 2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)
        pieces = [np.array([(u1-pars.limbdark_1)/.1, (u2-pars.limbdark_2)/.1])]
        r, b, dur = s['ROR_0'], s['IMPACT_0'], s['DUR14_0']
        quarter_pieces = []
        for q, original_time, warped, legx, f, err in groups:
            if q in excluded:
                continue
            t = (warped+s['C0_0']+s['C1_0']*legx)[:, None]+offsets
            phase = np.fmod(t.ravel(), ephem.period)
            phase[phase > ephem.period/2] -= ephem.period
            z = np.sqrt(b*b+4*((1+r)**2-b*b)*(phase/dur)**2)
            z[abs(phase) > dur] = 100.
            model = _quadratic_ld._quadratic_ld(z, abs(r), u1, u2, 1).reshape(-1, oversample).mean(axis=1)
            photometric = (model-f)/err
            pieces.append(photometric)
            quarter_pieces.append((q, original_time, photometric))
        if return_by_quarter:
            return quarter_pieces
        return np.concatenate(pieces)
    rows = []
    for idx in chosen:
        s = samples[idx]
        res = residual(s)
        value = -.5*np.dot(res, res)
        rows.append(dict(sample_index=int(idx), stored_lnlike=float(s['LN_LIKE']),
                         replay_lnlike=float(value), difference=float(value-s['LN_LIKE'])))
    a.output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(a.output/'replay.csv', index=False)
    passed = all(np.isclose(r['stored_lnlike'], r['replay_lnlike'], rtol=1e-12, atol=1e-7) for r in rows)
    summary = dict(input_hashes={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths+[a.ephemeris_source]},
                   oversample=oversample, selected_points=sum(len(g[4]) for g in groups), quarters=[g[0] for g in groups],
                   max_absolute_difference=max(abs(r['difference']) for r in rows),
                   replay_gate_passed=passed, absolute_tolerance=1e-7, relative_tolerance=1e-12,
                   scope='Single-planet LC white-noise replay. Current fork chord formula and saved Ephemeris source; historical runtime not attested. No posterior replacement.')
    (a.output/'provenance.json').write_text(json.dumps(summary, indent=2))
    print(pd.DataFrame(rows).to_string(index=False))
    print(json.dumps(summary, indent=2))
    if not passed:
        raise ValueError('Saved likelihood replay failed; do not proceed to inference')
    if a.residual_report:
        from native_residual_metrics import write_residual_report
        best_index = int(np.argmax(samples['LN_LIKE']))
        best = samples[best_index]
        quarter_path, global_path, report = write_residual_report(
            residual(best, return_by_quarter=True), a.output, exposure=exposure,
            sample_index=best_index, stored_lnlike=best['LN_LIKE'])
        summary['residual_report'] = dict(
            quarter_metrics=str(quarter_path), global_metrics=str(global_path),
            sample_index=best_index, scope=report['scope'])
        (a.output/'provenance.json').write_text(json.dumps(summary, indent=2))
        print(f'wrote residual diagnostics to {quarter_path} and {global_path}')
    if a.profile or a.leave_out:
        from scipy.optimize import least_squares
        from scipy.special import logsumexp
        names = ['C0_0', 'C1_0', 'ROR_0', 'DUR14_0', 'LD_Q1', 'LD_Q2']
        best = samples[np.argmax(samples['LN_LIKE'])]
        weights = np.exp(samples['LN_WT']-logsumexp(samples['LN_WT']))
        if a.leave_out:
            cases = [('all', [])]+[(f'without_q{q}', [q]) for q in [3, 6, 10]]
            cases += [(f'without_season{s}', [g[0] for g in groups if g[0] % 4 == s]) for s in range(4)]
            results = []
            for label, excluded in cases:
                for impact in [0., .96]:
                    def objective(x):
                        s = dict(zip(names, x)); s['IMPACT_0'] = impact
                        return residual(s, excluded)
                    candidates = []
                    for hours in [4.5, 5.5, 6.5]:
                        initial = np.array([best[k] for k in names]); initial[3] = hours/24
                        fit = least_squares(objective, initial,
                            bounds=([-1., -1., 1e-5, 58.848777/86400, 1e-8, 1e-8],
                                    [1., 1., .99, 3*duration, 1-1e-8, 1-1e-8]),
                            diff_step=1e-4, x_scale='jac', max_nfev=300)
                        candidates.append(fit)
                    fit = min(candidates, key=lambda f:f.cost)
                    row = dict(case=label, excluded_quarters=','.join(map(str, excluded)), impact=impact,
                               points=sum(len(g[4]) for g in groups if g[0] not in excluded),
                               lnlike=float(-fit.cost), success=bool(fit.success),
                               converged_starts=sum(bool(f.success) for f in candidates))
                    row.update(dict(zip(names, map(float, fit.x))))
                    results.append(row)
                    print(row, flush=True)
                    pd.DataFrame(results).to_csv(a.output/'leave_out.csv', index=False)
            if not a.profile:
                return
        profiles = []
        for impact in [0., .5, .8, .96]:
            def fixed(x):
                s = dict(zip(names, x)); s['IMPACT_0'] = impact
                return residual(s)
            for start_hours in [4.5, 5.5, 6.5]:
                initial = np.array([best[k] for k in names])
                initial[3] = start_hours/24
                fit = least_squares(fixed, initial,
                    bounds=([-1., -1., 1e-5, 58.848777/86400, 1e-8, 1e-8],
                            [1., 1., .99, 3*duration, 1-1e-8, 1-1e-8]),
                    diff_step=1e-4, x_scale='jac', max_nfev=300)
                row = dict(impact=impact, start_hours=start_hours, success=bool(fit.success),
                           evaluations=int(fit.nfev), lnlike=float(-fit.cost),
                           delta_from_stored_best=float(-fit.cost-best['LN_LIKE']),
                           archived_mass_within_005=float(weights[abs(samples['IMPACT_0']-impact)<.05].sum()))
                row.update(dict(zip(names, map(float, fit.x))))
                profiles.append(row)
                print(row, flush=True)
                pd.DataFrame(profiles).to_csv(a.output/'profile.csv', index=False)
        high = max((r for r in profiles if r['impact'] == .96), key=lambda r:r['lnlike'])
        low = max((r for r in profiles if r['impact'] == 0.), key=lambda r:r['lnlike'])
        residuals = []
        for solution in [low, high]:
            s = {k:solution[k] for k in names}; s['IMPACT_0'] = solution['impact']
            residuals.append(residual(s))
        contributions = []
        start = 0
        for label, length in [('ld_penalty', 2)]+[(str(g[0]), len(g[4])) for g in groups]:
            stop = start+length
            ll_low, ll_high = [-.5*np.dot(r[start:stop], r[start:stop]) for r in residuals]
            contributions.append(dict(component=label, points=length, low_lnlike=float(ll_low),
                                      high_lnlike=float(ll_high), high_minus_low=float(ll_high-ll_low)))
            start = stop
        pd.DataFrame(contributions).to_csv(a.output/'quarter_contrast.csv', index=False)
        if not np.isclose(sum(r['high_minus_low'] for r in contributions), high['lnlike']-low['lnlike'], atol=1e-7):
            raise ValueError('Quarter contrast does not reconstruct total likelihood difference')


if __name__ == '__main__':
    main()
