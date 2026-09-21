"""Measure out-of-transit residual scales without changing noise assumptions."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from audit_transit_timing_fold import nearest_phase


def robust_scale(x):
    return float(1.4826*np.median(abs(x-np.median(x))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    timing = args.directory/'K02712_00_quick.ttvs'
    centers = np.loadtxt(timing)[:, 2]
    rows, hashes = [], {timing.name: hashlib.sha256(timing.read_bytes()).hexdigest()}
    for stage in ['detrended', 'filtered']:
        path = args.directory/f'K02712_lc_{stage}.fits'
        hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        with fits.open(path) as hdul:
            time, flux, error, quarter = [np.asarray(hdul[key].data) for key in ['TIME','FLUX','ERROR','QUARTER']]
        finite = np.isfinite(time+flux+error) & (error > 0)
        out = finite & (abs(nearest_phase(time, centers)) > 8)
        for q in np.unique(quarter[out]):
            selected = out & (quarter == q)
            f, t, err = flux[selected], time[selected], error[selected]
            order = np.argsort(t)
            f, t, err = f[order], t[order], err[order]
            if len(f) < 10:
                continue
            close = (np.diff(t) > 0) & (np.diff(t) < .03)
            differences = np.diff(f)[close]
            normalized = differences/np.sqrt(err[:-1][close]**2+err[1:][close]**2)
            if len(differences) < 5:
                continue
            rows.append(dict(stage=stage, quarter=int(q), n=len(f), adjacent_pairs=len(differences),
                median_error_ppm=float(np.median(err)*1e6), scatter_ppm=robust_scale(f)*1e6,
                adjacent_difference_scatter_ppm=robust_scale(differences)*1e6/np.sqrt(2),
                normalized_adjacent_scale=robust_scale(normalized)))
    args.output.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(args.output/'quarter_noise.csv', index=False)
    (args.output/'provenance.json').write_text(json.dumps(dict(input_sha256=hashes,
        scope='Outside 8 hours of saved centers. Adjacent differences suppress slowly correlated noise; not a GP validation or revised likelihood.'), indent=2)+'\n')
    print(result.to_string(index=False))
    print(result.groupby('stage')[['median_error_ppm','scatter_ppm','normalized_adjacent_scale']].median().to_string())


if __name__ == '__main__':
    main()
