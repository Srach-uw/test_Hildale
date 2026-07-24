from __future__ import annotations

import numpy as np
import pytest

from hierarchical_rayleigh import fit_from_mass_matrix, leave_one_planet_influence, rayleigh_grid


def test_zero_outlier_floor_preserves_plain_rayleigh() -> None:
    e_grid = np.linspace(0.0, 0.95, 300)
    sigmas = np.array([0.03, 0.1, 0.3])
    plain, _ = rayleigh_grid(e_grid, sigmas, False, "none")
    explicit_zero, _ = rayleigh_grid(e_grid, sigmas, False, "none", 0.0)
    np.testing.assert_allclose(explicit_zero, plain, rtol=0.0, atol=0.0)


def test_positive_outlier_floor_adds_finite_tail_and_renormalizes() -> None:
    e_grid = np.linspace(0.0, 0.95, 300)
    sigmas = np.array([0.01, 0.03])
    rays, _ = rayleigh_grid(e_grid, sigmas, False, "none", 1e-6)
    assert np.all(rays[-1] > 0.0)
    np.testing.assert_allclose(np.trapezoid(rays, e_grid, axis=0), 1.0, rtol=1e-10)


def test_outlier_floor_rejects_invalid_values() -> None:
    e_grid = np.linspace(0.0, 0.95, 20)
    sigmas = np.array([0.1])
    for value in (-1e-6, np.nan, 1.1):
        with pytest.raises(ValueError, match="outlier_floor"):
            rayleigh_grid(e_grid, sigmas, False, "none", value)


def test_outlier_floor_reduces_single_extreme_planet_leverage() -> None:
    e_grid = np.linspace(0.0, 0.95, 300)
    sigmas = np.linspace(1e-4, 0.5, 1000)
    low_e = np.exp(-0.5 * ((e_grid - 0.02) / 0.01) ** 2)
    high_e = np.exp(-0.5 * ((e_grid - 0.85) / 0.01) ** 2)
    masses = np.vstack([low_e] * 40 + [high_e])

    plain = fit_from_mass_matrix(masses, e_grid, sigmas, False, "none", 0.0)
    robust = fit_from_mass_matrix(masses, e_grid, sigmas, False, "none", 1e-6)

    assert robust["expected_e"] < plain["expected_e"]
    assert robust["outlier_floor"] == pytest.approx(1e-6)


def test_modified_model_mean_and_low_sigma_boundary_are_explicit() -> None:
    e_grid = np.linspace(0.0, 0.95, 300)
    sigmas = np.linspace(1e-4, 0.1, 2000)
    low_e = np.exp(-0.5 * ((e_grid - 0.003) / 0.002) ** 2)
    masses = np.vstack([low_e] * 100)

    fit = fit_from_mass_matrix(masses, e_grid, sigmas, False, "none", 1e-6)

    assert fit["expected_e_model"] == pytest.approx(fit["expected_e_truncated"])
    assert fit["expected_e_model"] > fit["expected_e"]
    assert fit["sigma_near_grid_lower_edge"]
    assert fit["boundary_flag"]


def test_leave_one_planet_influence_has_the_correct_direction() -> None:
    import pandas as pd

    e_grid = np.linspace(0.0, 0.95, 300)
    sigmas = np.linspace(1e-4, 0.5, 1000)
    low_e = np.exp(-0.5 * ((e_grid - 0.02) / 0.01) ** 2)
    high_e = np.exp(-0.5 * ((e_grid - 0.5) / 0.03) ** 2)
    masses = np.vstack([low_e] * 20 + [high_e])
    summary = pd.DataFrame(
        {
            "kepoi_name": [f"K{i:05d}.01" for i in range(len(masses))],
            "disk": ["thin"] * len(masses),
            "system": ["single"] * len(masses),
            "e50": [0.02] * 20 + [0.5],
        }
    )

    influence = leave_one_planet_influence(
        summary,
        masses,
        e_grid,
        sigmas,
        False,
        "thin_singles",
        "none",
    )

    high_row = influence.iloc[-1]
    assert high_row["fractional_shift_when_removed"] < 0.0
