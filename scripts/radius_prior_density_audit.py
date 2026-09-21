"""Compare transit-density posteriors under explicit radius-prior reweighting.

No stellar catalog or eccentricity extraction is used. Existing nested samples
cannot reveal modes they missed; effective sample size is only a local check.
"""

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from duration_measure_control import density_and_derivative
from extract_eccentricity_posteriors_direct import (
    _sample_frame, nested_sample_weights, paired_period_samples, weighted_quantile,
)


def prior_weights(weights, radii):
    """Apply the joint uniform/log-uniform ratio on identical bounds."""
    radii = np.asarray(radii, float)
    if radii.ndim != 2 or radii.shape[0] != len(weights):
        raise ValueError('Radius samples must be row-paired')
    if not np.all(np.isfinite(radii) & (radii >= 1e-5) & (radii <= .99)):
        raise ValueError('Radius samples outside the declared prior support')
    logratio = np.log(radii).sum(axis=1)
    result = weights * np.exp(logratio-logratio.max())
    return result/result.sum()


def audit(path):
    rows = []
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with fits.open(path) as hdul:
        samples = _sample_frame(hdul['SAMPLES'].data)
        indices = sorted(int(c[4:]) for c in samples if c.startswith('ROR_'))
        base = nested_sample_weights(samples.LN_WT.to_numpy(float), 'dynesty')
        all_radii = samples[[f'ROR_{i}' for i in indices]].to_numpy(float)
        joint = prior_weights(base, all_radii)
        for index in indices:
            r = samples[f'ROR_{index}'].to_numpy(float)
            b = samples[f'IMPACT_{index}'].to_numpy(float)
            t = samples[f'DUR14_{index}'].to_numpy(float)
            period = paired_period_samples(samples, hdul, index)
            with np.errstate(invalid='ignore', divide='ignore'):
                density, _ = density_and_derivative(t, 0., 0., period, r, b)
            valid = ((t > 0) & (t < period/2) & (b >= 0) & (b < 1+r)
                     & np.isfinite(density) & (density > 0))
            for mode, weights in [('loguniform', base),
                                  ('uniform_focal', prior_weights(base, r[:, None])),
                                  ('uniform_all_radii', joint)]:
                retained = weights[valid].sum()
                if retained <= 0:
                    raise ValueError(f'No valid density mass: {path.name}, planet {index}')
                w = weights[valid]/retained
                quantiles = weighted_quantile(np.log10(density[valid]), w, [.025, .16, .5, .84, .975])
                rows.append(dict(target=path.name.removesuffix('-results.fits'),
                                 planet_index=index, period_days=float(np.median(period)),
                                 mode=mode, log_rho_p025=quantiles[0], log_rho_p16=quantiles[1],
                                 log_rho_p50=quantiles[2], log_rho_p84=quantiles[3],
                                 log_rho_p975=quantiles[4], halfwidth68=(quantiles[3]-quantiles[1])/2,
                                 width95=quantiles[4]-quantiles[0], ess=float(1/(w@w)),
                                 retained_weight=float(retained), fits_sha256=digest))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fits-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for path in sorted(args.fits_dir.glob('*-results.fits')):
        rows.extend(audit(path))
    if not rows:
        raise ValueError('No FITS samples found')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(f'Wrote {len(rows)} rows for {len(rows)//3} fitted planets')


if __name__ == '__main__':
    main()
