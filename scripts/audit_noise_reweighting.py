"""Conditional common-error-scale sensitivity, preserving limb-darkening penalty."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.special import logsumexp

from audit_influential_transit_shapes import quantiles
from duration_measure_control import density_and_derivative
from extract_eccentricity_posteriors_direct import _sample_frame, paired_period_samples


def reweight(logwt, loglike, ld_penalty, scale):
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('Invalid error scale')
    arrays = [np.asarray(x, float) for x in [logwt, loglike, ld_penalty]]
    if any(a.shape != arrays[0].shape for a in arrays) or not all(np.isfinite(a).all() for a in arrays):
        raise ValueError('Invalid sample arrays')
    phot_like = arrays[1]-arrays[2]
    updated = arrays[0]+(scale**-2-1)*phot_like
    return np.exp(updated-logsumexp(updated))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    path = args.directory/'K02712-results.fits'
    catalog = args.directory/'K02712_transit_parameters.csv'
    pars = pd.read_csv(catalog).iloc[0]
    with fits.open(path) as hdul:
        samples = _sample_frame(hdul['SAMPLES'].data)
        period = paired_period_samples(samples, hdul, 0)
    q1, q2 = samples.LD_Q1.to_numpy(float), samples.LD_Q2.to_numpy(float)
    u1, u2 = 2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)
    ld = -.5*((u1-pars.limbdark_1)**2+(u2-pars.limbdark_2)**2)/.01
    t, r, b = [samples[key].to_numpy(float) for key in ['DUR14_0','ROR_0','IMPACT_0']]
    rho, _ = density_and_derivative(t, 0, 0, period, r, b)
    if not np.isfinite(rho).all():
        raise ValueError('Invalid geometry in saved samples')
    rows = []
    for scale in [1., 1.2, 1.405731, 2., 3.]:
        w = reweight(samples.LN_WT, samples.LN_LIKE, ld, scale)
        row = dict(error_scale=scale, ess=float(1/(w@w)), max_weight=float(w.max()),
                   probability_b_below_0p8=float(w[b < .8].sum()))
        for name, values in [('impact',b), ('duration_hours',24*t), ('rho_circ_solar',rho)]:
            for label, value in zip(['p16','p50','p84'],quantiles(values,w)):
                row[name+'_'+label] = float(value)
        rows.append(row)
    args.output.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output/'sensitivity.csv',index=False)
    provenance = {'inputs': {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [path,catalog]},
      'assumptions': 'LN_LIKE is white photometric likelihood plus sigma=0.1 limb-darkening penalty, as in pinned helper. Common error scale, not GP or quarter-dependent noise.',
      'limitations': 'Reuses nested quadrature. ESS is diagnostic, not proof of coverage or convergence. No evidence comparison across scales; Gaussian normalization constants omitted because they cancel within each fixed scale.'}
    (args.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()
