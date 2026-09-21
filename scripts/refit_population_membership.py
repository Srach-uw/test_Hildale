"""Apply an existing membership ledger to saved diagnostic likelihoods."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from hierarchical_rayleigh import fit_from_mass_matrix
from real_duration_measure_audit import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mass-dir', type=Path, required=True)
    parser.add_argument('--membership', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    allowed = set(pd.read_csv(args.membership).kepoi_name)
    rows, ledger = [], []
    for path in sorted(args.mass_dir.glob('*.npz')):
        label = path.stem
        mode = next(m for m in ['manuscript_reciprocal', 'legacy_forward_norm', 'none'] if label.endswith('_'+m))
        with np.load(path) as data:
            keep = np.array([p in allowed for p in data['planets']])
            if not keep.any():
                raise ValueError('No planets remain')
            fit = fit_from_mass_matrix(data['mass_matrix'][keep], data['e_grid'], data['sigma_grid'], mode != 'none', mode)
            rows.append(dict(source=label, n=int(keep.sum()), membership_sha256=sha256(args.membership), **fit))
            for planet, retained in zip(data['planets'], keep):
                ledger.append(dict(source=label, planet=planet, retained=bool(retained)))
    args.output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output/'population.csv', index=False)
    pd.DataFrame(ledger).to_csv(args.output/'membership.csv', index=False)
    print(pd.DataFrame(rows)[['source', 'n', 'expected_e']].to_string(index=False))


if __name__ == '__main__':
    main()
