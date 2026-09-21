"""Retrieve public MAST cadence files with hashes and FITS identity checks."""
import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import urljoin
import requests
from astropy.io import fits
import numpy as np


class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]
    def handle_starttag(self,tag,attrs):
        if tag=='a':
            self.links.extend(v for k,v in attrs if k=='href')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kepid',type=int,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--download',action='store_true')
    p.add_argument('--cadence',choices=['short','long'],default='short')
    p.add_argument('--filenames',nargs='+',help='Restrict to exact filenames present in the MAST listing')
    a=p.parse_args()
    kid=f'{a.kepid:09d}'
    url=f'https://archive.stsci.edu/pub/kepler/lightcurves/{kid[:4]}/{kid}/'
    response=requests.get(url,timeout=45); response.raise_for_status()
    links=Links(); links.feed(response.text)
    suffix='_slc.fits' if a.cadence=='short' else '_llc.fits'
    names=sorted(set(x for x in links.links if x.startswith('kplr'+kid+'-') and x.endswith(suffix) and '/' not in x))
    if a.filenames:
        missing=set(a.filenames)-set(names)
        if missing:
            raise ValueError(f'Requested files absent from cadence listing: {sorted(missing)}')
        names=sorted(set(a.filenames))
    a.output.mkdir(parents=True,exist_ok=True)
    rows=[]
    for name in names:
        row=dict(filename=name,url=urljoin(url,name))
        if a.download:
            path=a.output/name
            if not path.exists():
                data=requests.get(row['url'],timeout=90); data.raise_for_status()
                temporary=path.with_suffix('.fits.partial')
                temporary.write_bytes(data.content)
                with fits.open(temporary) as h:
                    if int(h[0].header['KEPLERID'])!=a.kepid:
                        raise ValueError('KIC identity mismatch')
                temporary.rename(path)
            with fits.open(path) as h:
                if int(h[0].header['KEPLERID'])!=a.kepid:
                    raise ValueError('Cached KIC identity mismatch')
                d=h[1].data
                valid=np.isfinite(d['TIME']) & np.isfinite(d['PDCSAP_FLUX']) & np.isfinite(d['PDCSAP_FLUX_ERR']) & (d['PDCSAP_FLUX_ERR']>0)
                if valid.sum()<2:
                    raise ValueError(f'Insufficient finite PDCSAP rows: {name}')
                if a.cadence not in str(h[0].header.get('OBSMODE','')).lower():
                    raise ValueError(f'Unexpected observing mode: {name}')
                row.update(quarter=int(h[0].header['QUARTER']),obsmode=h[0].header.get('OBSMODE'),
                    finite_pdcsap_rows=int(valid.sum()),
                    median_spacing_seconds=float(np.median(np.diff(d['TIME'][valid]))*86400),
                    time_min=float(d['TIME'][valid].min()),time_max=float(d['TIME'][valid].max()),
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        rows.append(row)
        print(row,flush=True)
    manifest='mast_sc_manifest.json' if a.cadence=='short' else 'mast_lc_manifest.json'
    (a.output/manifest).write_text(json.dumps(dict(kepid=a.kepid,cadence=a.cadence,listing=url,
        listing_sha256=hashlib.sha256(response.content).hexdigest(),files=rows,
        caveat='Finite PDCSAP availability, not quality-filtered transit coverage.'),indent=2))


if __name__=='__main__':
    main()
