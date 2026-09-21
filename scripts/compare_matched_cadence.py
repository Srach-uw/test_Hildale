"""Compare conditional geometry grids without treating them as posteriors."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    tables, hashes = [], {}
    for label, folder in [
        ('raw_lc', 'K02712_matched_lc_20260919'),
        ('raw_sc', 'K02712_matched_sc_20260919'),
        ('detrended_lc', 'K02712_matched_detrended_lc_20260919'),
    ]:
        path = args.root / folder / 'fits.csv'
        table = pd.read_csv(path)
        if not table.success.all() or table.event_indices.nunique() != 1:
            raise ValueError(f'Invalid fit or event coverage: {label}')
        table['source'] = label
        table['delta_chi_square'] = table.chi_square - table.groupby('timing').chi_square.transform('min')
        tables.append(table)
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    combined = pd.concat(tables, ignore_index=True)
    if combined.event_indices.nunique() != 1:
        raise ValueError('Event sets do not agree across data sources')
    combined.to_csv(args.root / 'K02712_same_event_comparison_20260919.csv', index=False)
    caveat = ('Same events, not identical timestamps across sources. Detrended LC has '
              '549 points versus raw LC 531. Fixed geometry grids, not posterior '
              'probabilities; do not compare absolute chi-square between sources.')
    (args.root / 'K02712_same_event_comparison_20260919.json').write_text(
        json.dumps(dict(hashes=hashes, caveat=caveat), indent=2))
    print(combined[['source', 'timing', 'fixed_impact', 'delta_chi_square', 'rho_circ_solar', 'points']].to_string(index=False))


if __name__ == '__main__':
    main()
