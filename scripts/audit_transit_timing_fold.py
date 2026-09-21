"""Inspect timing sensitivity using identical saved detrended photometry."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


def nearest_phase(time, centers):
    centers = np.sort(centers)
    right = np.searchsorted(centers, time).clip(0, len(centers)-1)
    left = (right-1).clip(0, len(centers)-1)
    index = np.where(abs(time-centers[left]) < abs(time-centers[right]), left, right)
    return (time-centers[index])*24


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    lcpath = args.directory/'K02712_lc_detrended.fits'
    ttpath = args.directory/'K02712_00_quick.ttvs'
    with fits.open(lcpath) as f:
        time, flux, error = [np.asarray(f[key].data, float) for key in ['TIME','FLUX','ERROR']]
    good = np.isfinite(time+flux+error) & (error > 0)
    time, flux, error = time[good], flux[good], error[good]
    tt = np.loadtxt(ttpath)
    slope, intercept = np.polyfit(tt[:, 0], tt[:, 2], 1)
    linear = intercept+slope*tt[:, 0]
    phases = {'Constant period': nearest_phase(time, linear),
              'Regularized transit times': nearest_phase(time, tt[:, 2])}
    edges = np.linspace(-8, 8, 49)
    fig, axes = plt.subplots(2, 1, figsize=(9, 7))
    table = []
    for label, phase in phases.items():
        x, y, se = [], [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            keep = (phase >= lo) & (phase < hi)
            if not keep.any():
                continue
            w = 1/error[keep]**2
            mean = np.average(flux[keep], weights=w)
            uncertainty = np.sqrt(1/w.sum())
            x.append((lo+hi)/2); y.append((mean-1)*1e6); se.append(uncertainty*1e6)
            table.append([label, (lo+hi)/2, mean, uncertainty, int(keep.sum())])
        axes[0].errorbar(x, y, yerr=se, label=label, marker='.', linewidth=1)
    axes[0].set(xlabel='Hours from transit center', ylabel='Relative flux (ppm)', title='KOI-2712: same detrending, different timing folds')
    axes[0].legend()
    axes[1].plot(tt[:, 0], (tt[:, 1]-linear)*24, '.', alpha=.5, label='Individual estimates')
    axes[1].plot(tt[:, 0], (tt[:, 2]-linear)*24, '.', label='Regularized estimates')
    axes[1].set(xlabel='Transit index', ylabel='Timing offset (hours)')
    axes[1].legend()
    fig.tight_layout()
    args.output.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output/'timing_fold.png', dpi=150)
    import pandas as pd
    pd.DataFrame(table, columns=['fold','hours','flux','formal_standard_error','n']).to_csv(args.output/'fold_bins.csv',index=False)
    result = {'n_photometry': len(time), 'n_transit_times': len(tt), 'linear_period_days': slope,
              'regularized_timing_rms_hours': float(np.std((tt[:,2]-linear)*24)),
              'individual_timing_rms_hours': float(np.std((tt[:,1]-linear)*24)),
              'regularized_timing_range_hours': float(np.ptp(tt[:,2]-linear)*24),
              'inputs': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [lcpath,ttpath]},
              'scope': 'Same detrended data; not independent preprocessing. Errors ignore correlated noise. No hypothesis test or exclusion.'}
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
