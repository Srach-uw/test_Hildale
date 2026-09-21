"""Locate Monte Carlo changes in a fixed-cohort population likelihood."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp

def contrast(matrix, egrid, low, high, mode):
    if mode not in {'none', 'legacy_forward_norm', 'manuscript_reciprocal'}:
        raise ValueError('Unsupported selection mode')
    if min(low, high) <= 0 or np.any(matrix < 0) or not np.isfinite(matrix).all():
        raise ValueError('Invalid likelihood inputs')
    widths = np.gradient(egrid)
    widths[[0, -1]] *= .5
    positive = egrid > 0
    with np.errstate(divide='ignore'):
        logmass = np.log(matrix)
    logterms = []
    for sigma in [low, high]:
        logray = np.full(len(egrid), -np.inf)
        logray[positive] = np.log(egrid[positive])-2*np.log(sigma)-.5*(egrid[positive]/sigma)**2
        logray -= logsumexp(logray+np.log(widths))
        term = logsumexp(logmass+logray[None, :], axis=1)
        if mode == 'legacy_forward_norm':
            term -= logsumexp(logray+np.log(widths)-np.log1p(-egrid**2))
        logterms.append(term)
    if not np.isfinite(logterms).all():
        raise ValueError('Zero or invalid likelihood support; contrast is undefined')
    return logterms[1]-logterms[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=Path, required=True)
    parser.add_argument('--second', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sigma-low', type=float, default=.06)
    parser.add_argument('--sigma-high', type=float, default=.11)
    parser.add_argument('--selection-mode', default='manuscript_reciprocal')
    args = parser.parse_args()
    with np.load(args.first) as first, np.load(args.second) as second:
        if not np.array_equal(first['planets'], second['planets']):
            raise ValueError('Planet membership or ordering differs')
        np.testing.assert_array_equal(first['e_grid'], second['e_grid'])
        a = contrast(first['mass_matrix'], first['e_grid'], args.sigma_low, args.sigma_high, args.selection_mode)
        b = contrast(second['mass_matrix'], second['e_grid'], args.sigma_low, args.sigma_high, args.selection_mode)
        rows = pd.DataFrame(dict(planet=first['planets'], first_log_high_over_low=a,
                                  second_log_high_over_low=b, contrast_change=b-a,
                                  absolute_change=abs(b-a)))
    rows = rows.sort_values('absolute_change', ascending=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(args.output, index=False)
    print(rows.head(10).to_string(index=False))
    print(f'Total log-contrast: first={a.sum():.6f}, second={b.sum():.6f}')


if __name__ == '__main__':
    main()
