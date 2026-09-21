"""Compare direct density integration with duration-posterior reweighting.

Geometry is fixed and the duration likelihood is Gaussian. These are numerical
integrals, independent of ALDERAAN sampling and real-data catalog choices.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import roots_legendre

from common import trapezoid
from extract_eccentricity_posteriors_direct import DAY_S, G_SI, RHO_SUN_KG_M3


def duration_from_density(rho, e, w, period=12.0, ror=0.05, b=0.5):
    a2 = (G_SI * RHO_SUN_KG_M3 * rho * (period*DAY_S)**2 / (3*np.pi))**(2/3)
    g = (1+e*np.sin(w))/np.sqrt(1-e**2)
    return period/np.pi * np.arcsin(np.sqrt(((1+ror)**2-b*b)/(a2-b*b))) / g


def density_and_derivative(t, e, w, period=12.0, ror=0.05, b=0.5):
    g = (1+e*np.sin(w))/np.sqrt(1-e**2)
    k = np.pi/period*g
    angle = t*k
    chord = (1+ror)**2-b*b
    a2 = chord/np.sin(angle)**2+b*b
    c = 3*np.pi/(G_SI*(period*DAY_S)**2*RHO_SUN_KG_M3)
    rho = c*a2**1.5
    derivative = 3*c*np.sqrt(a2)*chord*np.cos(angle)/np.sin(angle)**3*k
    return rho, derivative


def calculate_case(density_error, duration_fractional_error, order=384):
    egrid = np.linspace(0, 0.8, 161)
    omega = np.linspace(0, 2*np.pi, 120, endpoint=False)[:, None]
    nodes, weights = roots_legendre(order)
    lo, hi = 0.02, 1+9*density_error
    rho = (lo+hi)/2 + nodes*(hi-lo)/2
    rho_weights = weights*(hi-lo)/2
    central_t = float(duration_from_density(1.0, 0, 0))
    terr = central_t*duration_fractional_error
    curves = {key: [] for key in ["direct", "log_duration_uncorrected", "flat_duration_uncorrected", "measure_corrected"]}
    for e in egrid:
        t = duration_from_density(rho[None, :], e, omega)
        direct = np.exp(-0.5*((t-central_t)/terr)**2-0.5*((rho-1)/density_error)**2)
        curves["direct"].append(float(np.mean(direct @ rho_weights)))
        tlo = duration_from_density(hi, e, omega)
        thi = duration_from_density(lo, e, omega)
        tnodes = (tlo+thi)/2 + nodes*(thi-tlo)/2
        tw = weights*(thi-tlo)/2
        implied, derivative = density_and_derivative(tnodes, e, omega)
        product = np.exp(-0.5*((tnodes-central_t)/terr)**2-0.5*((implied-1)/density_error)**2)
        curves["log_duration_uncorrected"].append(float(np.mean(np.sum(product/tnodes*tw, axis=1))))
        curves["flat_duration_uncorrected"].append(float(np.mean(np.sum(product*tw, axis=1))))
        curves["measure_corrected"].append(float(np.mean(np.sum(product*derivative*tw, axis=1))))
    curves = {key: np.asarray(value)/trapezoid(value, egrid) for key, value in curves.items()}
    rows = []
    for key, pdf in curves.items():
        cdf = np.r_[0, np.cumsum((pdf[1:]+pdf[:-1])/2*np.diff(egrid))]
        rows.append({"density_error_fraction": density_error,
                     "duration_error_fraction": duration_fractional_error,
                     "method": key, "e50": float(np.interp(.5, cdf, egrid)),
                     "mean_e": float(trapezoid(pdf*egrid, egrid)),
                     "total_variation_from_direct": float(.5*trapezoid(abs(pdf-curves['direct']), egrid))})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for density_error in [.1, .5]:
        for duration_error in [.03, .2]:
            rows.extend(calculate_case(density_error, duration_error))
            print(f"Finished density={density_error}, duration={duration_error}", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
