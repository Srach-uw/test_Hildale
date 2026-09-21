"""Cross-check influential planets against frozen catalog inputs, without cuts."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['contrast', 'catalog', 'inventory', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    scores = pd.read_csv(args.contrast).sort_values('log_likelihood_high_over_low', ascending=False).head(5)
    catalog = pd.read_csv(args.catalog, comment='#')
    inventory = pd.read_csv(args.inventory)
    catcols = ['kepoi_name', 'kepid', 'kepler_name', 'koi_disposition', 'koi_pdisposition',
               'koi_prad', 'koi_srad', 'koi_slogg', 'koi_duration', 'koi_impact', 'koi_model_snr']
    invcols = ['kepoi_name', 'berger_logg', 'berger_rad', 'berger2018_evol', 'berger2018_bin']
    rows = scores.merge(catalog[catcols], on='kepoi_name', validate='one_to_one', how='left')
    rows = rows.merge(inventory[invcols], on='kepoi_name', validate='one_to_one', how='left')
    if rows.kepid.isna().any() or rows.berger_rad.isna().any():
        raise ValueError('Missing catalog matches')
    rows['berger_to_koi_stellar_radius'] = rows.berger_rad/rows.koi_srad
    # Estimating a diagnostic density from catalog gravity and radius.
    rows['koi_gravity_radius_density_solar'] = 10**(rows.koi_slogg-4.438)/rows.koi_srad
    rows['circular_to_koi_gravity_radius_density'] = rows.rho_circular_median_solar/rows.koi_gravity_radius_density_solar
    rows['adopted_to_koi_density'] = rows.rho_true_solar/rows.koi_gravity_radius_density_solar
    if not np.isfinite(rows.adopted_to_koi_density).all():
        raise ValueError('Invalid density comparison')
    args.output.mkdir(parents=True, exist_ok=True)
    rows.to_csv(args.output/'catalog_crosscheck.csv', index=False)
    provenance = {key: {'path': str(getattr(args, key)),
                       'sha256': hashlib.sha256(getattr(args, key).read_bytes()).hexdigest()}
                  for key in ['contrast', 'catalog', 'inventory']}
    provenance['interpretation'] = ('Top five by saved likelihood contrast; not a rejection rule. '
        'KOI gravity-radius density is a diagnostic without propagated errors, not a replacement prior. '
        'Frozen 2024 catalog dispositions are not a current disposition audit.')
    (args.output/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(rows[['kepoi_name', 'berger_to_koi_stellar_radius',
                'rho_circular_to_catalog_ratio', 'circular_to_koi_gravity_radius_density']].to_string(index=False))


if __name__ == '__main__':
    main()
