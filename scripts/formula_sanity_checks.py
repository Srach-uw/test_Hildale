from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from common import output_dir
from extract_eccentricity_posteriors import posterior_grid_from_zeta
from hierarchical_rayleigh import fit_rayleigh, transit_selection_correction


def main() -> None:
    out = output_dir() / "formula_sanity_checks"
    out.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, object]] = []

    e_grid = np.linspace(0.001, 0.95, 240)
    omega_grid = np.linspace(0.0, 2.0 * np.pi, 180)

    checks.extend(check_transit_selection(e_grid, omega_grid))
    checks.extend(check_zeta_posteriors(e_grid, omega_grid))
    checks.extend(check_rayleigh_recovery(out, e_grid, omega_grid))

    status = "pass" if all(row["passed"] for row in checks) else "fail"
    payload = {"status": status, "checks": checks}
    path = out / "formula_sanity_checks.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print(f"\nWrote: {path}")
    if status != "pass":
        raise SystemExit(1)


def check_transit_selection(e_grid: np.ndarray, omega_grid: np.ndarray) -> list[dict[str, object]]:
    corr = transit_selection_correction(np.array([0.001, 0.5]), np.array([np.pi / 2, 3 * np.pi / 2]))
    return [
        {
            "name": "transit_selection_finite_positive",
            "passed": bool(np.all(np.isfinite(corr)) and np.all(corr > 0)),
            "details": corr.tolist(),
        },
        {
            "name": "transit_selection_periastron_less_than_apastron",
            "passed": bool(corr[1, 0] < corr[1, 1]),
            "details": {"e_0p5_periastron": float(corr[1, 0]), "e_0p5_apastron": float(corr[1, 1])},
        },
    ]


def check_zeta_posteriors(e_grid: np.ndarray, omega_grid: np.ndarray) -> list[dict[str, object]]:
    rng = np.random.default_rng(123)
    zeta_circular = rng.normal(1.0, 0.03, size=1200)
    zeta_short = rng.normal(0.6, 0.03, size=1200)
    post_circular = posterior_grid_from_zeta(zeta_circular, e_grid, omega_grid, include_transit_prior=False)
    post_short = posterior_grid_from_zeta(zeta_short, e_grid, omega_grid, include_transit_prior=False)
    med_circular = e_median(post_circular, e_grid)
    med_short = e_median(post_short, e_grid)
    return [
        {
            "name": "zeta_one_prefers_lower_e_than_short_duration",
            "passed": bool(med_circular < med_short),
            "details": {"median_e_zeta_1": med_circular, "median_e_zeta_0p6": med_short},
        },
        {
            "name": "posterior_grids_normalized",
            "passed": bool(np.isclose(post_circular.sum(), 1.0) and np.isclose(post_short.sum(), 1.0)),
            "details": {"sum_zeta_1": float(post_circular.sum()), "sum_zeta_0p6": float(post_short.sum())},
        },
    ]


def check_rayleigh_recovery(out: Path, e_grid: np.ndarray, omega_grid: np.ndarray) -> list[dict[str, object]]:
    files = []
    e0 = 0.035
    sigma_e = 0.008
    e_profile = np.exp(-0.5 * ((e_grid - e0) / sigma_e) ** 2)
    omega_profile = np.ones_like(omega_grid)
    posterior = e_profile[:, None] * omega_profile[None, :]
    posterior = posterior / posterior.sum()
    e_pdf = posterior.sum(axis=1)
    for i in range(12):
        path = out / f"synthetic_low_e_{i:02d}.npz"
        np.savez_compressed(
            path,
            e_grid=e_grid,
            omega_grid=omega_grid,
            posterior=posterior,
            e_pdf=e_pdf,
            kepid=1000000 + i,
            kepoi_name=f"KTEST{i:02d}.01",
            koi_target=f"KTEST{i:02d}",
            include_transit_prior=False,
        )
        files.append(str(path))
    sigmas = np.linspace(0.005, 0.12, 200)
    fit = fit_rayleigh(files, sigmas, apply_transit_selection=True)
    return [
        {
            "name": "rayleigh_fit_low_e_synthetic_reasonable",
            "passed": bool(0.015 <= fit["expected_e"] <= 0.08),
            "details": fit,
        }
    ]


def e_median(posterior: np.ndarray, e_grid: np.ndarray) -> float:
    e_mass = posterior.sum(axis=1)
    cdf = np.cumsum(e_mass)
    cdf = cdf / cdf[-1]
    return float(np.interp(0.5, cdf, e_grid))


if __name__ == "__main__":
    main()
