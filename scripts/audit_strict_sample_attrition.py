"""Reconcile strict-summary omissions against exact-ID quality records."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd


def boolean(series):
    values = series.astype(str).str.lower()
    if not values.isin(['true', 'false']).all():
        raise ValueError('Missing or invalid boolean flags')
    return values.eq('true')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['inventory', 'summary', 'qc', 'output']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    inv, summary, qc = [pd.read_csv(path) for path in [a.inventory,a.summary,a.qc]]
    for frame in [inv, summary, qc]:
        if frame.kepoi_name.isna().any() or frame.kepoi_name.duplicated().any():
            raise ValueError('Invalid identifiers')
    if set(inv.kepoi_name) != set(qc.kepoi_name):
        raise ValueError('Inventory and QC identities differ')
    if not set(summary.kepoi_name) <= set(inv.kepoi_name):
        raise ValueError('Summary outside inventory')
    ledger = inv[['kepoi_name','disk','system']].merge(qc, on='kepoi_name', validate='one_to_one')
    grazing = boolean(ledger.gilbert_grazing_exclude)
    radius = boolean(ledger.gilbert_radius_precision_exclude)
    available = boolean(ledger.gilbert_qc_available)
    present = ledger.kepoi_name.isin(summary.kepoi_name)
    if not available.all() or not (present == ~(grazing | radius)).all():
        raise ValueError('Summary omissions not fully explained by QC')
    excluded_ids = set(summary.loc[boolean(summary.qc_primary_exclude), 'kepoi_name'])
    ledger['stage'] = 'retained'
    ledger.loc[grazing & ~radius, 'stage'] = 'grazing_only'
    ledger.loc[~grazing & radius, 'stage'] = 'radius_precision_only'
    ledger.loc[grazing & radius, 'stage'] = 'grazing_and_radius_precision'
    ledger.loc[ledger.kepoi_name.isin(excluded_ids), 'stage'] = 'summary_qc_excluded'
    counts = pd.crosstab([ledger.disk,ledger.system], ledger.stage)
    counts['inventory'] = counts.sum(axis=1)
    a.output.mkdir(parents=True, exist_ok=True)
    ledger.to_csv(a.output/'planet_ledger.csv', index=False)
    counts.to_csv(a.output/'counts.csv')
    (a.output/'scope.json').write_text(json.dumps(dict(
        inputs={str(path):hashlib.sha256(path.read_bytes()).hexdigest()
                for path in [a.inventory,a.summary,a.qc]},
        warning='Available QC metadata does not independently verify raw FITS integrity or justify cuts as Sagear-equivalent.'),indent=2))
    print(counts.to_string())


if __name__ == '__main__':
    main()
