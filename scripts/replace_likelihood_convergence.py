"""Test one planet's integration convergence without changing cohort membership."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from audit_population_seed_leverage import contrast
from hierarchical_rayleigh import fit_from_mass_matrix


def replace_row(planets, matrix, planet, replacement):
    indices = np.flatnonzero(planets == planet)
    if len(indices) != 1 or replacement.shape != matrix[indices[0]].shape:
        raise ValueError('Replacement must identify exactly one existing planet and match its grid')
    result = matrix.copy()
    result[indices[0]] = replacement
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-dir', type=Path, required=True)
    parser.add_argument('--replacement-dir', type=Path, required=True)
    parser.add_argument('--planet', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for base_path in sorted(args.base_dir.glob('*_density_likelihood_only_*.npz')):
        suffix = base_path.name.split('_', 1)[1]
        for replacement_path in sorted(args.replacement_dir.glob('*_'+suffix)):
            with np.load(base_path) as base, np.load(replacement_path) as replacement:
                np.testing.assert_array_equal(base['e_grid'], replacement['e_grid'])
                if replacement['planets'].tolist() != [args.planet]:
                    raise ValueError('Replacement archive does not contain exactly the requested planet')
                mode = suffix.split('density_likelihood_only_', 1)[1].removesuffix('.npz')
                matrix = replace_row(base['planets'], base['mass_matrix'], args.planet,
                                     replacement['mass_matrix'][0])
                fit = fit_from_mass_matrix(matrix, base['e_grid'], base['sigma_grid'], mode != 'none', mode)
                delta = contrast(replacement['mass_matrix'], base['e_grid'], .06, .11, mode)[0]
                rows.append(dict(base_seed=base_path.name.split('_')[0],
                                 replacement_seed=replacement_path.name.split('_')[0],
                                 n=len(base['planets']), replaced_planet=args.planet,
                                 planet_log_contrast=delta, **fit))
    if not rows:
        raise ValueError('No corresponding saved mass matrices')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(pd.DataFrame(rows).groupby('selection_mode').expected_e.agg(['min', 'median', 'max']).to_string())


if __name__ == '__main__':
    main()
