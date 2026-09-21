"""Measure raw SC coverage around saved transit times without fitting flux."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-dir',type=Path,required=True)
    p.add_argument('--timings',type=Path,required=True)
    p.add_argument('--duration-hours',type=float,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.duration_hours<=0:
        raise ValueError('Positive duration required')
    timings=np.loadtxt(a.timings,ndmin=2)
    times=[]; clean=[]; hashes={}; headers=[]
    for path in sorted(a.raw_dir.glob('*_slc.fits')):
        hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        with fits.open(path) as h:
            d=h[1].data
            valid=np.isfinite(d['TIME']) & np.isfinite(d['PDCSAP_FLUX']) & np.isfinite(d['PDCSAP_FLUX_ERR']) & (d['PDCSAP_FLUX_ERR']>0)
            times.extend(d['TIME'][valid]); clean.extend(d['SAP_QUALITY'][valid]==0)
            headers.append(dict(file=path.name,bjdrefi=h[1].header.get('BJDREFI'),
                                bjdreff=h[1].header.get('BJDREFF'),timesys=h[1].header.get('TIMESYS')))
    times=np.asarray(times); clean=np.asarray(clean,dtype=bool)
    if len(times)==0 or len(np.unique(times))!=len(times):
        raise ValueError('Empty or duplicate raw times')
    half=a.duration_hours/48
    rows=[]
    for index,independent,regularized,*_ in timings:
        for label,tt in [('independent',independent),('regularized',regularized)]:
            offset=times-tt
            inside=np.abs(offset)<=half
            if not inside.any():
                continue
            ingress=(offset>=-half)&(offset<-.6*half)
            egress=(offset<=half)&(offset>.6*half)
            rows.append(dict(transit_index=int(index),timing=label,center=tt,
                in_transit=int(inside.sum()),unflagged_in_transit=int((inside&clean).sum()),
                unflagged_ingress_proxy=int((ingress&clean).sum()),
                unflagged_egress_proxy=int((egress&clean).sum()),
                unflagged_left_baseline=int(((offset>=-3*half)&(offset<-half)&clean).sum()),
                unflagged_right_baseline=int(((offset<=3*half)&(offset>half)&clean).sum())))
    a.output.mkdir(parents=True,exist_ok=True)
    frame=pd.DataFrame(rows)
    frame.to_csv(a.output/'transit_coverage.csv',index=False)
    scope=dict(raw_hashes=hashes,timing_sha256=hashlib.sha256(a.timings.read_bytes()).hexdigest(),
        duration_hours=a.duration_hours,raw_headers=headers,total_finite=len(times),
        quality_zero=int(clean.sum()),
        caveats='Saved timings omit rejected events. Quality==0 is a conservative diagnostic, not ALDERAAN mask. Outer 20 percent windows are geometry proxies, not fitted contacts. Assumes timing-file BKJD matches FITS reference; headers recorded.')
    (a.output/'scope.json').write_text(json.dumps(scope,indent=2))
    print(frame.to_string(index=False))


if __name__=='__main__':
    main()
