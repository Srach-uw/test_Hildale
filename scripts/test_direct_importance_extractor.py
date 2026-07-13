from __future__ import annotations

import numpy as np

from extract_eccentricity_posteriors_direct import (
    DAY_S,
    G_SI,
    RHO_SUN_KG_M3,
    direct_importance_posterior,
    macdougall_rho_star_samp,
    weighted_posterior_grid,
    weighted_quantile,
)


def exact_duration_days(period_days: float, rho_solar: float, ror: float, impact: float) -> float:
    period_s = period_days * DAY_S
    a_over_r = (G_SI * rho_solar * RHO_SUN_KG_M3 * period_s**2 / (3.0 * np.pi)) ** (1.0 / 3.0)
    argument = np.sqrt(((1.0 + ror) ** 2 - impact**2) / (a_over_r**2 - impact**2))
    return period_days * np.arcsin(argument) / np.pi


def test_exact_macdougall_equation_sanity() -> None:
    period_days = 12.3
    rho_true = 1.27
    ror = 0.041
    impact = 0.37
    duration_days = exact_duration_days(period_days, rho_true, ror, impact)
    recovered = macdougall_rho_star_samp(
        period_days * DAY_S,
        np.array([duration_days * DAY_S]),
        np.array([ror]),
        np.array([impact]),
        np.array([0.0]),
        np.array([1.234]),
    )[0]
    assert np.isclose(recovered, rho_true, rtol=2e-13, atol=0.0)
    assert np.isnan(
        macdougall_rho_star_samp(
            period_days * DAY_S,
            np.array([duration_days * DAY_S]),
            np.array([ror]),
            np.array([1.2]),
            np.array([0.0]),
            np.array([0.0]),
        )[0]
    )


def test_weighted_quantiles() -> None:
    values = np.array([0.0, 1.0, 2.0, 3.0])
    weights = np.array([1.0, 1.0, 6.0, 2.0])
    got = weighted_quantile(values, weights, [0.25, 0.5, 0.75])
    assert np.allclose(got, [9.0 / 7.0, 2.0, 2.625])


def test_weighted_grid_normalization() -> None:
    e_grid = np.linspace(0.0, 0.95, 20)
    omega_grid = np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False)
    posterior = weighted_posterior_grid(
        np.array([0.01, 0.22, 0.22, 0.91]),
        np.array([-0.1, 0.2, 2.0 * np.pi + 0.2, 4.0]),
        np.array([1.0, 2.0, 3.0, 4.0]),
        e_grid,
        omega_grid,
    )
    assert posterior.shape == (20, 24)
    assert np.all(np.isfinite(posterior))
    assert np.all(posterior >= 0.0)
    assert np.isclose(posterior.sum(), 1.0)
    assert np.count_nonzero(posterior) == 3


def synthetic_result(seed: int = 991, impact_center: float = 0.2) -> dict[str, object]:
    period_days = 10.0
    rho_true = 1.0
    ror = np.full(250, 0.05)
    impact = np.linspace(impact_center - 0.025, impact_center + 0.025, len(ror))
    duration = np.array([exact_duration_days(period_days, rho_true, ror_i, impact_i) for ror_i, impact_i in zip(ror, impact)])
    nested_weights = np.linspace(1.0, 2.0, len(ror))
    nested_weights /= nested_weights.sum()
    e_grid = np.linspace(0.0, 0.95, 80)
    omega_grid = np.linspace(0.0, 2.0 * np.pi, 72, endpoint=False)
    return direct_importance_posterior(
        ror,
        impact,
        duration,
        nested_weights,
        period_days,
        rho_true,
        0.02,
        0.02,
        n_proposals=120_000,
        e_max=0.95,
        e_grid=e_grid,
        omega_grid=omega_grid,
        density_error_mode="symmetric-average",
        seed=seed,
    )


def test_deterministic_output() -> None:
    first = synthetic_result()
    second = synthetic_result()
    assert np.array_equal(first["posterior"], second["posterior"])
    assert np.array_equal(first["e_quantiles"], second["e_quantiles"])
    assert first["importance_ess"] == second["importance_ess"]


def test_synthetic_circular_recovery_across_impact() -> None:
    for index, impact_center in enumerate((0.0, 0.5, 0.8)):
        result = synthetic_result(seed=12345 + index, impact_center=impact_center)
        posterior = np.asarray(result["posterior"])
        e_mass = posterior.sum(axis=1)
        e_grid = np.linspace(0.0, 0.95, 80)
        map_e = e_grid[np.argmax(e_mass)]
        e16, e50, _ = np.asarray(result["e_quantiles"])
        assert np.isclose(posterior.sum(), 1.0)
        assert map_e <= 0.06
        assert e16 < 0.08
        assert e50 < 0.22
        assert result["importance_ess"] > 100.0


def test_split_density_likelihood_includes_scale_normalization() -> None:
    from extract_eccentricity_posteriors_direct import density_log_likelihood

    got = density_log_likelihood(np.array([0.9, 1.2]), 1.0, 0.2, 0.1, "split")
    expected = np.array([-0.5 - np.log(0.1), -0.5 - np.log(0.2)])
    assert np.allclose(got, expected)
