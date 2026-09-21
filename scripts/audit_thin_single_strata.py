"""Describe thin-single tension without changing selection or refitting means."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def audit(summary, inventory, leverage):
    for frame in (summary, inventory, leverage):
        if frame.kepoi_name.isna().any() or frame.kepoi_name.duplicated().any():
            raise ValueError('Missing or duplicate planet identifiers')
    columns = ['kepoi_name', 'berger_logg', 'berger_rad', 'berger2018_evol',
               'koi_model_snr', 'koi_prad', 'transit_fit_source']
    data = leverage[['kepoi_name', 'log_contrast_fit_over_paper']].merge(
        summary, on='kepoi_name', validate='one_to_one', how='left', indicator=True)
    if not data['_merge'].eq('both').all():
        raise ValueError('Leverage rows absent from posterior summary')
    data = data.drop(columns='_merge').merge(
        inventory[columns], on='kepoi_name', validate='one_to_one', how='left', indicator=True)
    if not data['_merge'].eq('both').all():
        raise ValueError('Posterior rows absent from inventory')
    if not (data.disk.eq('thin') & data.system.eq('single')).all():
        raise ValueError('Expected retained thin singles only')
    data = data.drop(columns='_merge')
    ratio = pd.to_numeric(data.rho_circular_to_catalog_ratio, errors='coerce')
    if not (np.isfinite(ratio) & ratio.gt(0)).all():
        raise ValueError('Invalid density ratio')
    data['delta_dex'] = np.log10(ratio)
    data['radius_proxy_earth'] = data.ror50 * data.berger_rad * 109.1
    data['logg_bin'] = pd.cut(data.berger_logg, [-np.inf, 4, 4.4, np.inf],
                              labels=['<=4', '4-4.4', '>4.4']).astype(str)
    data['radius_bin'] = pd.cut(data.radius_proxy_earth, [-np.inf, 3.5, 8, np.inf],
                                labels=['<=3.5', '3.5-8', '>8']).astype(str)
    data['snr_bin'] = pd.cut(data.koi_model_snr, [-np.inf, 20, 50, 100, np.inf],
                             labels=['<=20', '20-50', '50-100', '>100']).astype(str)
    rows = []
    for field in ['logg_bin', 'radius_bin', 'snr_bin', 'berger2018_evol', 'transit_fit_source']:
        for label, group in data.groupby(field, dropna=False):
            rows.append(dict(stratum=field, label=str(label), n=len(group),
                delta_median=group.delta_dex.median(),
                abs_delta_median=group.delta_dex.abs().median(),
                fraction_ratio_above_one=group.delta_dex.gt(0).mean(),
                individual_e50_median=group.e50.median(),
                net_log_contrast=group.log_contrast_fit_over_paper.sum(),
                positive_contrast_fraction=group.log_contrast_fit_over_paper.gt(0).mean()))
    return data, pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['summary', 'inventory', 'leverage', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    inputs = {name: getattr(args, name) for name in ['summary', 'inventory', 'leverage']}
    data, groups = audit(*(pd.read_csv(path) for path in inputs.values()))
    args.output.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.output / 'planet_diagnostics.csv', index=False)
    groups.to_csv(args.output / 'strata.csv', index=False)
    scope = dict(n=len(data), inputs={name: dict(path=str(path),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for name, path in inputs.items()},
        caveats=['Fit-source labels are historical inventory labels, not verified raw FITS provenance.',
                 'Radius proxy uses posterior median radius ratio and Berger radius; no uncertainty propagation.',
                 'Individual eccentricity medians are not hierarchical population means.',
                 'Strata are diagnostic only; no new selection cuts or canonical writes.'])
    (args.output / 'scope.json').write_text(json.dumps(scope, indent=2), encoding='utf-8')
    print(groups.to_string(index=False))


if __name__ == '__main__':
    main()
