"""Verify identical posterior inputs before interpreting selection sensitivity."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['before', 'after', 'output']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    before, after = pd.read_csv(a.before), pd.read_csv(a.after)
    if any(d.kepoi_name.duplicated().any() or d.kepoi_name.isna().any() for d in [before,after]):
        raise ValueError('Invalid IDs')
    matched = after.merge(before,on='kepoi_name',suffixes=('_after','_before'),validate='one_to_one')
    if len(matched) != len(after):
        raise ValueError('After is not an exact subset')
    checked = []
    for c in before.columns.intersection(after.columns):
        if c == 'kepoi_name':
            continue
        x,y = matched[c+'_before'], matched[c+'_after']
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            equal = np.allclose(x,y,rtol=1e-12,atol=1e-15,equal_nan=True)
        else:
            equal = x.fillna('').astype(str).equals(y.fillna('').astype(str))
        if not equal:
            raise ValueError('Input difference: '+c)
        checked.append(c)
    paths = sorted(set(matched.posterior_file_before))
    hashes = {}
    for path in paths:
        hashes[path] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    result = dict(before_rows=len(before), after_rows=len(after), matched_rows=len(matched),
        checked_columns=checked,shared_posterior_hashes=hashes,
        input_hashes={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in [a.before,a.after]})
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2))
    print(f'PASS: {len(matched)} matched rows; {len(checked)} fields; {len(hashes)} posterior hashes')


if __name__ == '__main__':
    main()
