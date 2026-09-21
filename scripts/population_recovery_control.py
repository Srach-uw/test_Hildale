"""Recover known populations from a noisy photoeccentric density observable.

This isolates hierarchy and transit selection. It does not validate detrending,
transit fitting, or the transformation of ALDERAAN posterior samples.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from hierarchical_rayleigh import fit_from_mass_matrix


def draw_transiting_population(rng, n, mean, emax=0.6):
    accepted_e, accepted_w = [], []
    sigma = mean / np.sqrt(np.pi / 2.0)
    while sum(len(v) for v in accepted_e) < n:
        e = rng.rayleigh(sigma, size=n * 3)
        w = rng.uniform(0.0, 2.0 * np.pi, size=len(e))
        support = e < emax
        e, w = e[support], w[support]
        transit = (1.0 + e * np.sin(w)) / (1.0 - e**2)
        keep = rng.uniform(size=len(e)) < transit * (1.0 - emax)
        accepted_e.append(e[keep])
        accepted_w.append(w[keep])
    return np.concatenate(accepted_e)[:n], np.concatenate(accepted_w)[:n]


def density_observable(e, omega):
    return 3.0 * np.log10((1.0 + e * np.sin(omega)) / np.sqrt(1.0 - e**2))


def likelihood_masses(observed, error, egrid, omega):
    predicted = density_observable(egrid[:, None], omega[None, :])
    transit = (1.0 + egrid[:, None] * np.sin(omega)) / (1.0 - egrid[:, None]**2)
    widths = np.gradient(egrid)
    widths[[0, -1]] *= 0.5
    matrices = {mode: [] for mode in ["none", "legacy_forward_norm", "manuscript_reciprocal"]}
    for y in observed:
        loglike = -0.5 * ((y - predicted) / error)**2
        likelihood = np.exp(loglike - loglike.max())
        posterior = likelihood * widths[:, None]
        posterior /= posterior.sum()
        matrices["none"].append(posterior.sum(axis=1))
        matrices["legacy_forward_norm"].append((posterior * transit).sum(axis=1))
        matrices["manuscript_reciprocal"].append((posterior / transit).sum(axis=1))
    return {mode: np.asarray(rows) for mode, rows in matrices.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--planets", type=int, default=600)
    parser.add_argument("--seeds", type=int, default=2)
    args = parser.parse_args()
    egrid = np.linspace(0, 0.6, 401)
    omega = np.linspace(0, 2*np.pi, 180, endpoint=False)
    sigmas = np.linspace(0.002, 0.3, 350)
    rows = []
    for seed in range(args.seeds):
        for truth in [0.022, 0.066, 0.2]:
            for error in [0.03, 0.3]:
                rng = np.random.default_rng(20260915 + 1000*seed + int(truth*10000) + int(error*100))
                e, w = draw_transiting_population(rng, args.planets, truth)
                observed = density_observable(e, w) + rng.normal(0, error, len(e))
                matrices = likelihood_masses(observed, error, egrid, omega)
                for mode, matrix in matrices.items():
                    fit = fit_from_mass_matrix(matrix, egrid, sigmas, mode != "none", mode)
                    rows.append({"seed": seed, "n": args.planets, "true_mean": truth,
                                 "density_error_dex": error, "selection_mode": mode,
                                 "observed_true_mean": float(e.mean()),
                                 "recovered_mean": fit["expected_e"],
                                 "lower": fit["expected_e_lo"], "upper": fit["expected_e_hi"],
                                 "boundary": fit["boundary_flag"]})
                print(f"Completed seed={seed}, mean={truth}, error={error}", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = pd.DataFrame(rows)
    results.to_csv(args.output_dir / "recovery.csv", index=False)
    results["covers_truth"] = (results.lower <= results.true_mean) & (results.upper >= results.true_mean)
    summary = results.groupby(["true_mean", "density_error_dex", "selection_mode"]).agg(
        median_recovered=("recovered_mean", "median"),
        minimum=("recovered_mean", "min"), maximum=("recovered_mean", "max"),
        coverage=("covers_truth", "mean"),
    ).reset_index()
    summary.to_csv(args.output_dir / "summary.csv", index=False)
    (args.output_dir / "scope.json").write_text(json.dumps({
        "observable": "3 log10((1+e sin(omega))/sqrt(1-e^2)) plus Gaussian measurement error",
        "selection": "forward geometric transit probability, rejection sampled",
        "prior": "uniform e and omega for individual likelihood grids",
        "limits": "Observable-level experiment. No ALDERAAN, light curves, or post-fit exclusions.",
        "e_max": 0.6, "sigma_min": float(sigmas.min()), "sigma_max": float(sigmas.max()),
        "n_e": len(egrid), "n_omega": len(omega), "n_sigma": len(sigmas)
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
