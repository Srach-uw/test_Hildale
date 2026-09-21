"""Compare clipped and log-domain likelihoods without changing saved inputs."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from hierarchical_rayleigh import rayleigh_grid, posterior_weights_from_ll, weighted_quantile


def stable_terms(matrix, grid, sigmas, mode):
    if mode not in {'none', 'manuscript_reciprocal', 'legacy_forward_norm'}:
        raise ValueError('Unsupported mode')
    if (matrix.ndim != 2 or matrix.shape[1] != len(grid)
            or len(grid) < 2 or not np.isfinite(matrix).all()
            or np.any(matrix < 0) or not np.isfinite(grid).all()
            or np.any(np.diff(grid) <= 0) or grid[0] < 0 or grid[-1] >= 1
            or not np.isfinite(sigmas).all() or np.any(sigmas <= 0)):
        raise ValueError('Invalid inputs')
    widths = np.gradient(grid)
    widths[[0, -1]] *= .5
    with np.errstate(divide='ignore'):
        logmass = np.log(matrix)
        loge = np.log(grid)
    terms = []
    for sigma in sigmas:
        logr = loge - 2*np.log(sigma) - .5*(grid/sigma)**2
        logr -= logsumexp(logr + np.log(widths))
        term = logsumexp(logmass + logr, axis=1)
        if mode == 'legacy_forward_norm':
            term -= logsumexp(logr + np.log(widths) - np.log1p(-grid**2))
        if not np.isfinite(term).all():
            raise ValueError('Zero likelihood support')
        terms.append(term)
    return np.asarray(terms).T


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    mode = 'manuscript_reciprocal'
    with np.load(args.matrix) as data:
        matrix, grid = data['mass_matrix'], data['e_grid']
    sigmas = np.linspace(.002, 1, 600)
    stable = stable_terms(matrix, grid, sigmas, mode)
    rays, _ = rayleigh_grid(grid, sigmas, True, mode)
    legacy = np.log(np.clip(matrix @ rays, 1e-300, None))
    result = dict(input_sha256=hashlib.sha256(args.matrix.read_bytes()).hexdigest(),
                  input_file=str(args.matrix), n=len(matrix), selection_mode=mode)
    for name, terms in [('legacy', legacy), ('stable', stable)]:
        weights = posterior_weights_from_ll(sigmas, terms.sum(axis=0))
        result[name+'_mean_quantiles'] = list(weighted_quantile(
            sigmas*np.sqrt(np.pi/2), weights, [.16, .5, .84]))
    cold = np.array([.022/np.sqrt(np.pi/2)])
    exact = stable_terms(matrix, grid, cold, mode)[:, 0]
    rays, _ = rayleigh_grid(grid, cold, True, mode)
    clipped = np.log(np.clip((matrix @ rays)[:, 0], 1e-300, None))
    result['cold_terms_below_floor'] = int((exact < np.log(1e-300)).sum())
    result['cold_loglikelihood_overstatement'] = float((clipped-exact).sum())
    result['maximum_loglikelihood_overstatement_on_grid'] = float((legacy-stable).sum(axis=0).max())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
