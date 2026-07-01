from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit Rayleigh population models to e posterior grids.")
    parser.add_argument("--summary", default=None)
    parser.add_argument("--sigma-min", type=float, default=0.005)
    parser.add_argument("--sigma-max", type=float, default=0.25)
    parser.add_argument("--n-sigma", type=int, default=400)
    parser.add_argument(
        "--ignore-transit-selection",
        action="store_true",
        help="Diagnostic only: omit the reciprocal geometric transit-probability factor.",
    )
    args = parser.parse_args()

    summary_path = Path(args.summary) if args.summary else output_dir() / "eccentricity_posterior_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(f"Posterior summary not found: {summary_path}")
    summary = pd.read_csv(summary_path)
    sigmas = np.linspace(args.sigma_min, args.sigma_max, args.n_sigma)

    rows = []
    for disk, system, label in [
        ("thick", "single", "thick_singles"),
        ("thin", "single", "thin_singles"),
        ("thick", "multi", "thick_multis"),
        ("thin", "multi", "thin_multis"),
    ]:
        sub = summary[(summary["disk"] == disk) & (summary["system"] == system)]
        if len(sub) < 5:
            rows.append({"population": label, "n": len(sub), "status": "skipped_n_lt_5"})
            continue
        fit = fit_rayleigh(sub["posterior_file"].tolist(), sigmas, apply_transit_selection=not args.ignore_transit_selection)
        rows.append({"population": label, "n": len(sub), "status": "ok", **fit})

    out = pd.DataFrame(rows)
    suffix = "no_transit_selection" if args.ignore_transit_selection else "transit_selection"
    out_path = output_dir() / f"rayleigh_population_fit_{suffix}.csv"
    out.to_csv(output_dir() / "rayleigh_population_fit.csv", index=False)
    out.to_csv(out_path, index=False)
    print(out.to_string(index=False))
    print(f"\nWrote: {out_path}")


def fit_rayleigh(files: list[str], sigmas: np.ndarray, apply_transit_selection: bool = True) -> dict:
    e_masses = []
    e_grid = None
    for file in files:
        data = np.load(file)
        e_grid = data["e_grid"]
        omega_grid = data["omega_grid"]
        posterior = np.asarray(data["posterior"], dtype=float)
        posterior = np.clip(posterior, 0.0, np.inf)
        total = posterior.sum()
        if total <= 0:
            posterior = np.full_like(posterior, 1.0 / posterior.size)
        else:
            posterior = posterior / total
        if apply_transit_selection:
            selection = transit_selection_correction(e_grid, omega_grid)
            e_masses.append(np.sum(posterior * selection, axis=1))
        else:
            e_masses.append(np.sum(posterior, axis=1))
    assert e_grid is not None
    mass_matrix = np.vstack(e_masses)
    lls = []
    for sigma in sigmas:
        ray = (e_grid / sigma**2) * np.exp(-(e_grid**2) / (2.0 * sigma**2))
        ray = ray / np.trapz(ray, e_grid)
        terms = mass_matrix @ ray
        lls.append(np.sum(np.log(np.clip(terms, 1e-300, None))))
    lls = np.array(lls)
    best = int(np.argmax(lls))
    e_mean = sigmas[best] * np.sqrt(np.pi / 2.0)
    ok = lls > lls[best] - 0.5
    sig_ci = sigmas[ok]
    e_lo = sig_ci[0] * np.sqrt(np.pi / 2.0)
    e_hi = sig_ci[-1] * np.sqrt(np.pi / 2.0)
    return {
        "sigma_rayleigh": sigmas[best],
        "expected_e": e_mean,
        "expected_e_lo": e_lo,
        "expected_e_hi": e_hi,
        "ll_max": lls[best],
        "transit_selection_applied": bool(apply_transit_selection),
    }


def transit_selection_correction(e_grid: np.ndarray, omega_grid: np.ndarray) -> np.ndarray:
    e, omega = np.meshgrid(e_grid, omega_grid, indexing="ij")
    correction = (1.0 - e**2) / (1.0 + e * np.sin(omega))
    return np.clip(correction, 1e-12, np.inf)


if __name__ == "__main__":
    main()
