from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from common import circular_duration_days, output_dir, trapezoid
from extract_eccentricity_posteriors import posterior_grid_from_zeta
from hierarchical_rayleigh import fit_from_mass_matrix, rayleigh_grid, transit_probability_weight


def main() -> None:
    out = output_dir() / "formula_sanity_checks"
    out.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, object]] = []

    e_grid = np.linspace(0.001, 0.95, 240)
    omega_grid = np.linspace(0.0, 2.0 * np.pi, 180, endpoint=False)

    checks.extend(check_transit_selection(e_grid, omega_grid))
    checks.extend(check_exact_duration())
    checks.extend(check_zeta_posteriors(e_grid, omega_grid))
    checks.extend(check_rayleigh_recovery(e_grid, omega_grid))

    status = "pass" if all(row["passed"] for row in checks) else "fail"
    payload = {"status": status, "checks": checks}
    path = out / "formula_sanity_checks.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print(f"\nWrote: {path}")
    if status != "pass":
        raise SystemExit(1)


def check_transit_selection(e_grid: np.ndarray, omega_grid: np.ndarray) -> list[dict[str, object]]:
    weight = transit_probability_weight(np.array([0.001, 0.5]), np.array([np.pi / 2, 3 * np.pi / 2]))
    rays, normalizers = rayleigh_grid(
        e_grid,
        np.array([0.02, 0.05, 0.1]),
        apply_transit_selection=True,
        selection_mode="legacy_forward_norm",
    )
    uniform_posterior_mass = transit_probability_weight(e_grid, omega_grid).mean(axis=1) / len(e_grid)
    uninformative = np.vstack([uniform_posterior_mass] * 12)
    terms = (uninformative @ rays) / normalizers
    flat_after_normalization = np.nanmax(np.abs(terms / terms[:, :1] - 1.0))
    recovery = transit_selected_recovery(e_grid, omega_grid)
    return [
        {
            "name": "transit_probability_weight_finite_positive",
            "passed": bool(np.all(np.isfinite(weight)) and np.all(weight > 0)),
            "details": weight.tolist(),
        },
        {
            "name": "transit_probability_periastron_greater_than_apastron",
            "passed": bool(weight[1, 0] > weight[1, 1]),
            "details": {"e_0p5_periastron": float(weight[1, 0]), "e_0p5_apastron": float(weight[1, 1])},
        },
        {
            "name": "population_transit_normalization_flattens_uninformative_likelihood",
            "passed": bool(flat_after_normalization < 1e-2),
            "details": {"max_fractional_deviation": float(flat_after_normalization)},
        },
        {
            "name": "forward_selected_likelihood_recovers_intrinsic_rayleigh",
            "passed": bool(abs(recovery["forward_sigma_map"] - recovery["true_sigma"]) < 0.02),
            "details": recovery,
        },
        {
            "name": "manuscript_reciprocal_is_not_generative_selection_correction",
            "passed": bool(
                abs(recovery["reciprocal_sigma_map"] - recovery["true_sigma"])
                > abs(recovery["forward_sigma_map"] - recovery["true_sigma"])
            ),
            "details": recovery,
        },
    ]


def transit_selected_recovery(e_grid: np.ndarray, omega_grid: np.ndarray) -> dict[str, float]:
    """Demonstrate the distinction between paper replication and generative selection."""
    rng = np.random.default_rng(20260710)
    true_sigma = 0.30
    n_proposal = 250_000
    e = rng.rayleigh(true_sigma, size=n_proposal)
    omega = rng.uniform(0.0, 2.0 * np.pi, size=n_proposal)
    physical = e < 0.95
    e = e[physical]
    omega = omega[physical]
    ptrans = (1.0 + e * np.sin(omega)) / (1.0 - e**2)
    draw_prob = ptrans / ptrans.sum()
    chosen = rng.choice(np.arange(len(e)), size=3000, replace=True, p=draw_prob)
    e_obs = e[chosen]
    omega_obs = omega[chosen]
    e_idx = np.clip(np.searchsorted(e_grid, e_obs), 1, len(e_grid) - 1)
    left = e_grid[e_idx - 1]
    e_idx -= (np.abs(e_obs - left) < np.abs(e_obs - e_grid[e_idx])).astype(int)
    omega_idx = np.floor((omega_obs % (2.0 * np.pi)) / (2.0 * np.pi) * len(omega_grid)).astype(int)
    obs_weight = (1.0 + e_grid[e_idx] * np.sin(omega_grid[omega_idx])) / (1.0 - e_grid[e_idx] ** 2)

    forward_mass = np.zeros((len(chosen), len(e_grid)))
    reciprocal_mass = np.zeros_like(forward_mass)
    row = np.arange(len(chosen))
    forward_mass[row, e_idx] = obs_weight
    reciprocal_mass[row, e_idx] = 1.0 / obs_weight
    sigmas = np.linspace(0.05, 0.60, 800)
    forward = fit_from_mass_matrix(
        forward_mass,
        e_grid,
        sigmas,
        apply_transit_selection=True,
        selection_mode="legacy_forward_norm",
    )
    reciprocal = fit_from_mass_matrix(
        reciprocal_mass,
        e_grid,
        sigmas,
        apply_transit_selection=True,
        selection_mode="manuscript_reciprocal",
    )
    return {
        "true_sigma": true_sigma,
        "forward_sigma_map": float(forward["sigma_map"]),
        "reciprocal_sigma_map": float(reciprocal["sigma_map"]),
    }


def check_exact_duration() -> list[dict[str, object]]:
    period = np.array([10.0])
    ror = np.array([0.05])
    impact = np.array([0.0])
    a_over_r = np.array([20.0])
    got = circular_duration_days(period, impact, ror, a_over_r)[0]
    expected = (period[0] / np.pi) * np.arcsin((1.0 + ror[0]) / a_over_r[0])
    grazing = circular_duration_days(period, np.array([1.2]), ror, a_over_r)[0]
    return [
        {
            "name": "exact_circular_duration_matches_central_formula",
            "passed": bool(np.isclose(got, expected, rtol=0, atol=1e-14)),
            "details": {"got_days": float(got), "expected_days": float(expected)},
        },
        {
            "name": "exact_circular_duration_invalid_grazing_is_nan",
            "passed": bool(np.isnan(grazing)),
            "details": {"grazing_duration": None if np.isnan(grazing) else float(grazing)},
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


def check_rayleigh_recovery(e_grid: np.ndarray, omega_grid: np.ndarray) -> list[dict[str, object]]:
    e0 = 0.035
    sigma_e = 0.008
    e_profile = np.exp(-0.5 * ((e_grid - e0) / sigma_e) ** 2)
    posterior = e_profile[:, None] * np.ones_like(omega_grid)[None, :]
    posterior = posterior / posterior.sum()
    weight = transit_probability_weight(e_grid, omega_grid)
    mass = np.sum(posterior * weight, axis=1)
    mass_matrix = np.vstack([mass] * 12)
    sigmas = np.linspace(1e-4, 0.12, 300)
    fit = fit_from_mass_matrix(
        mass_matrix,
        e_grid,
        sigmas,
        apply_transit_selection=True,
        selection_mode="legacy_forward_norm",
    )
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
