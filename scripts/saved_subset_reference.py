"""Fit saved source posteriors on a recovery overlap and expose likelihood tensions."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_population_seed_leverage import contrast
from hierarchical_rayleigh import load_population_masses, fit_from_mass_matrix
from real_duration_measure_audit import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--fits-root', type=Path, required=True)
    parser.add_argument('--disk', required=True, choices=['thin', 'thick'])
    parser.add_argument('--system', required=True, choices=['single', 'multi'])
    parser.add_argument('--low-mean', type=float, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    summary = pd.read_csv(args.summary)
    flags = summary.qc_primary_exclude.astype(str).str.lower()
    if not flags.isin(['true', 'false']).all() or summary.kepoi_name.duplicated().any():
        raise ValueError('Invalid membership ledger')
    targets = {p.name.removesuffix('-results.fits') for p in args.fits_root.rglob('*-results.fits')}
    selected = summary.loc[(summary.disk == args.disk) & (summary.system == args.system)
                           & summary.koi_target.isin(targets) & (flags == 'false')].copy()
    if selected.empty:
        raise ValueError('Empty overlap')
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for mode in ['manuscript_reciprocal', 'legacy_forward_norm', 'none']:
        matrix, grid = load_population_masses(selected, mode != 'none', mode)
        fit = fit_from_mass_matrix(matrix, grid, np.linspace(.002, 1, 600), mode != 'none', mode)
        rows.append(dict(n=len(selected), **fit))
        scores = contrast(matrix, grid, args.low_mean/np.sqrt(np.pi/2), fit['sigma_rayleigh'], mode)
        columns = [c for c in ['kepoi_name', 'koi_period', 'rho_true_solar', 'rho_err_hi_solar',
                               'rho_err_lo_solar', 'rho_circular_to_catalog_ratio',
                               'rho_circular_median_solar', 'e16', 'e50', 'e84'] if c in selected]
        table = selected[columns].copy()
        table['log_likelihood_high_over_low'] = scores
        table['reference_low_mean'] = args.low_mean
        table['reference_high_mean'] = fit['expected_e']
        table.sort_values('log_likelihood_high_over_low', ascending=False).to_csv(
            args.output/f'planet_contrast_{mode}.csv', index=False)
    pd.DataFrame(rows).to_csv(args.output/'population.csv', index=False)
    selected[['kepoi_name', 'koi_target', 'disk', 'system']].to_csv(args.output/'membership.csv', index=False)
    (args.output/'scope.json').write_text(json.dumps(dict(summary_sha256=sha256(args.summary),
        n=len(selected), disk=args.disk, system=args.system,
        interpretation='Recovery overlap, saved source posteriors; likelihood contrasts are not posterior odds or grounds for exclusion.'), indent=2)+'\n')


if __name__ == '__main__':
    main()
