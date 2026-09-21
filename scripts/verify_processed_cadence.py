"""Validate processed cadence evidence before accepting a cadence experiment."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits


def inspect(path, cadence):
    with fits.open(path) as hdul:
        arrays = {name: np.asarray(hdul[name].data, dtype=float).ravel()
                  for name in ['TIME','FLUX','ERROR']}
    lengths = {len(v) for v in arrays.values()}
    if len(lengths) != 1 or min(lengths) < 3:
        raise ValueError('Mismatched or insufficient photometry')
    valid = np.logical_and.reduce([np.isfinite(v) for v in arrays.values()])
    valid &= arrays['ERROR'] > 0
    time = np.sort(arrays['TIME'][valid])
    if len(time) < 3:
        raise ValueError('Insufficient finite positive-error photometry')
    delta = np.diff(time)
    if np.any(delta <= 0):
        raise ValueError('Duplicate timestamps')
    median = float(np.median(delta)*86400)
    lower, upper = (30,120) if cadence == 'short' else (1200,2400)
    if not lower <= median <= upper:
        raise ValueError(f'{cadence} cadence mismatch: median spacing {median:.1f}s')
    return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                rows=len(valid),usable_rows=int(valid.sum()),median_spacing_seconds=median)


def verify_products(result_dir, target, mode, stages):
    records={}
    for stage in stages:
        records[stage]={}
        for short,cadence in [('lc','long'),('sc','short')]:
            path=result_dir/f'{target}_{short}_{stage}.fits'
            if path.exists():
                records[stage][cadence]=inspect(path,cadence)
        if mode=='both' and 'short' not in records[stage]:
            raise ValueError(f'SC treatment lacks verified {stage} short-cadence data')
        if mode=='long' and ('long' not in records[stage] or 'short' in records[stage]):
            raise ValueError('LC control must contain LC and no SC product; use a fresh directory')
    return records


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--result-dir',type=Path,required=True)
    p.add_argument('--target',required=True)
    p.add_argument('--mode',choices=['long','both'],required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--stage',choices=['detrended','filtered','both'],default='both')
    a=p.parse_args()
    stages=['detrended','filtered'] if a.stage=='both' else [a.stage]
    records=verify_products(a.result_dir,a.target,a.mode,stages)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(dict(target=a.target,mode=a.mode,products=records,
        limitation='Confirms processed data cadence, not actual sampler ingestion, transit coverage, or fit convergence.'),indent=2))
    print(json.dumps(records,indent=2))


if __name__=='__main__':
    main()
