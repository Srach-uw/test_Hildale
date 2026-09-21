"""Align raw and processed photometry exactly for a preprocessing control."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits


def match_unique_times(raw, processed):
    for time in (raw, processed):
        if not np.isfinite(time).all() or len(np.unique(time)) != len(time):
            raise ValueError('Times must be finite and unique')
    common, ri, pi = np.intersect1d(raw, processed, return_indices=True)
    if not len(common):
        raise ValueError('No exact timestamp intersection')
    return common, ri, pi


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-dir', type=Path, required=True)
    p.add_argument('--processed', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new diagnostic output directory')
    arrays, hashes = [], {}
    for path in sorted(a.raw_dir.glob('*_llc.fits')):
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        with fits.open(path) as h:
            d = h[1].data
            valid = (np.isfinite(d['TIME']) & np.isfinite(d['PDCSAP_FLUX']) &
                     np.isfinite(d['PDCSAP_FLUX_ERR']) & (d['PDCSAP_FLUX_ERR'] > 0) &
                     (d['SAP_QUALITY'] == 0))
            q = np.full(valid.sum(), h[0].header['QUARTER'])
            arrays.append(np.column_stack([d['TIME'][valid], d['PDCSAP_FLUX'][valid],
                                           d['PDCSAP_FLUX_ERR'][valid], q]))
    raw = np.concatenate(arrays)
    hashes[str(a.processed)] = hashlib.sha256(a.processed.read_bytes()).hexdigest()
    with fits.open(a.processed) as h:
        processed = np.column_stack([np.asarray(h[k].data).ravel() for k in ['TIME', 'FLUX', 'ERROR']])
    processed = processed[np.isfinite(processed).all(axis=1) & (processed[:, 2] > 0)]
    time, ri, pi = match_unique_times(raw[:, 0], processed[:, 0])
    raw, processed = raw[ri], processed[pi]
    rows = []
    for quarter in np.unique(raw[:, 3]):
        use = raw[:, 3] == quarter
        # Normalizing within each shared quarter.
        for data in (raw, processed):
            scale = np.median(data[use, 1])
            if scale <= 0:
                raise ValueError('Nonpositive flux normalization')
            data[use, 1:3] /= scale
        ratio = processed[use, 2] / raw[use, 2]
        rows.append(dict(quarter=int(quarter), points=int(use.sum()),
                         error_ratio_p16=float(np.quantile(ratio, .16)),
                         error_ratio_median=float(np.median(ratio)),
                         error_ratio_p84=float(np.quantile(ratio, .84))))
    a.output.mkdir(parents=True)
    for name, data in [('raw', raw), ('detrended', processed)]:
        hdus = [fits.PrimaryHDU()]
        for key, values in [('TIME', time), ('FLUX', data[:, 1]), ('ERROR', data[:, 2])]:
            hdus.append(fits.ImageHDU(values, name=key))
        fits.HDUList(hdus).writeto(a.output / f'{name}.fits')
    pd.DataFrame(rows).to_csv(a.output / 'uncertainty_comparison.csv', index=False)
    (a.output / 'provenance.json').write_text(json.dumps(dict(
        input_hashes=hashes, exact_shared_points=len(time),
        timestamp_sha256=hashlib.sha256(time.astype('<f8').tobytes()).hexdigest(),
        scope='Exact timestamp intersection after each input quality mask; quarterly median normalization. Diagnostic only; not posterior data.'), indent=2))
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()
