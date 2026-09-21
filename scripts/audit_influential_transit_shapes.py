"""Compare weighted paired transit geometry with frozen KOI catalog fits."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from duration_measure_control import density_and_derivative, duration_from_density
from extract_eccentricity_posteriors_direct import (
    _sample_frame, nested_sample_weights, paired_period_samples,
)


def quantiles(values, weights):
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cdf = (np.cumsum(weights)-weights/2)/weights.sum()
    return np.interp([.16, .5, .84], cdf, values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['summary', 'crosscheck', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--comparison-root', type=Path)
    args = parser.parse_args()
    summary = pd.read_csv(args.summary)
    cross = pd.read_csv(args.crosscheck)
    chosen = cross[['kepoi_name', 'koi_duration', 'koi_impact']].merge(
        summary, on='kepoi_name', validate='one_to_one')
    if len(chosen) != len(cross):
        raise ValueError('Missing summary match')
    rows, hashes = [], {}
    for _, row in chosen.iterrows():
        path = Path(row.raw_result_file)
        hashes[row.kepoi_name] = hashlib.sha256(path.read_bytes()).hexdigest()
        index = int(row.alderaan_planet_index)
        with fits.open(path, memmap=False) as hdul:
            samples = _sample_frame(hdul['SAMPLES'].data)
            period = paired_period_samples(samples, hdul, index)
        weight = nested_sample_weights(samples.LN_WT.to_numpy(float), 'dynesty')
        t, r, b = [samples[f'{key}_{index}'].to_numpy(float) for key in ['DUR14', 'ROR', 'IMPACT']]
        rho, _ = density_and_derivative(t, 0, 0, period, r, b)
        circular_t = duration_from_density(row.rho_true_solar, 0, 0, period, r, b)
        valid = np.isfinite(rho) & np.isfinite(circular_t) & (t > 0) & (b >= 0) & (b < 1+r)
        out = dict(kepoi_name=row.kepoi_name, catalog_duration_hours=row.koi_duration,
                   catalog_impact=row.koi_impact, valid_weight=float(weight[valid].sum()),
                   nested_ess=float(1/np.sum(weight**2)),
                   grazing_weight=float(weight[b > 1-r].sum()))
        for name, values in [('duration_hours', t*24), ('impact', b), ('ror', r),
                             ('circular_density_solar', rho),
                             ('duration_over_circular', t/circular_t)]:
            for label, value in zip(['p16', 'p50', 'p84'], quantiles(values[valid], weight[valid])):
                out[name+'_'+label] = float(value)
        out['duration_ratio_to_catalog'] = out['duration_hours_p50']/row.koi_duration
        rows.append(out)
    args.output.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(args.output/'transit_shapes.csv', index=False)
    provenance = {'fits_sha256': hashes, 'inputs': {
        name: hashlib.sha256(getattr(args, name).read_bytes()).hexdigest()
        for name in ['summary', 'crosscheck']},
        'scope': 'Five influential objects. Frozen catalog comparison, not fit validation or exclusion.'}
    (args.output/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    if args.comparison_root:
        comparisons = []
        for path in sorted(args.comparison_root.rglob('K02712-results.fits')):
            with fits.open(path, memmap=False) as hdul:
                sample = _sample_frame(hdul['SAMPLES'].data)
            weights = nested_sample_weights(sample.LN_WT.to_numpy(float), 'dynesty')
            entry = {'file': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
            for key, scale in [('DUR14_0', 24), ('IMPACT_0', 1)]:
                for label, value in zip(['p16', 'p50', 'p84'], quantiles(sample[key].to_numpy(float)*scale, weights)):
                    entry[key+'_'+label] = float(value)
            comparisons.append(entry)
        pd.DataFrame(comparisons).to_csv(args.output/'K02712_arm_geometry.csv', index=False)
        print(pd.DataFrame(comparisons).drop(columns='sha256').to_string(index=False))
    print(result[['kepoi_name', 'catalog_duration_hours', 'duration_hours_p50', 'catalog_impact',
                  'impact_p50', 'grazing_weight', 'duration_over_circular_p50']].to_string(index=False))


if __name__ == '__main__':
    main()
