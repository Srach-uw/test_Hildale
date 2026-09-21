"""Compare constant-speed chord and circular-orbit transit models."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dependency-path', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.dependency_path))
    import batman
    from batman import _quadratic_ld
    rows = []
    period = 7.507866572210437
    for hours in [5.2, 6.05]:
        duration = hours / 24
        for impact in [0., .5, .8, .96]:
            radius = .014
            p = batman.TransitParams()
            p.t0 = 0.; p.per = period; p.rp = radius
            p.a = np.sqrt(((1+radius)**2-impact**2)/np.sin(np.pi*duration/period)**2+impact**2)
            p.inc = np.degrees(np.arccos(impact/p.a))
            p.ecc = 0.; p.w = 90.; p.u = [.481, .149]; p.limb_dark = 'quadratic'
            time = np.linspace(-duration, duration, 2001)
            for cadence, exposure, oversample in [('short', 58.85/86400, 3), ('long', 1765.5/86400, 31)]:
                orbit = batman.TransitModel(p, time, exp_time=exposure, supersample_factor=oversample).light_curve(p)
                expanded = (time[:, None]+np.linspace(-exposure/2, exposure/2, oversample)).ravel()
                separation = np.sqrt(impact**2+4*((1+radius)**2-impact**2)*(expanded/duration)**2)
                separation[abs(expanded)>duration] = 100.
                chord = _quadratic_ld._quadratic_ld(separation, radius, *p.u, 1).reshape(-1, oversample).mean(axis=1)
                rows.append(dict(duration_hours=hours, impact=impact, cadence=cadence,
                                 max_flux_difference_ppm=float(np.max(abs(orbit-chord))*1e6),
                                 rms_flux_difference_ppm=float(np.sqrt(np.mean((orbit-chord)**2))*1e6)))
    args.output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output/'comparison.csv', index=False)
    (args.output/'scope.json').write_text(json.dumps(dict(period_days=period, radius_ratio=.014,
        limb_darkening=[.481,.149], batman_version=batman.__version__,
        source='https://raw.githubusercontent.com/gjgilbert/batman/master/c_src/_rsky.c',
        scope='Published current fork formula, not attested historical runtime. Fixed-parameter flux differences, not refitted posterior shifts.'), indent=2))
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()
