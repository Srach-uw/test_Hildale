from __future__ import annotations

import numpy as np
import pandas as pd
from astropy.io import fits

from extract_eccentricity_posteriors_direct import (
    DAY_S,
    G_SI,
    RHO_SUN_KG_M3,
    adjust_stellar_density,
    direct_importance_posterior,
    macdougall_rho_star_samp,
    nested_sample_weights,
    paired_period_samples,
    periastron_eccentricity_cut,
    resampled_posterior_grid,
    weighted_posterior_grid,
    weighted_quantile,
)


def test_periastron_eccentricity_cut_matches_gilbert_condition() -> None:
    period_days = np.array([10.0, 10.1, 9.9])
    rho_solar = 1.2
    period_s = period_days * DAY_S
    a_over_r = (
        G_SI * rho_solar * RHO_SUN_KG_M3 * period_s**2 / (3.0 * np.pi)
    ) ** (1.0 / 3.0)
    expected = 1.0 - 1.0 / np.mean(a_over_r)
    got = periastron_eccentricity_cut(period_days, rho_solar, ceiling=0.99)
    assert np.isclose(got, expected)
    assert periastron_eccentricity_cut(period_days, rho_solar, ceiling=0.5) == 0.5


def test_periastron_eccentricity_cut_rejects_invalid_inputs() -> None:
    with np.testing.assert_raises(ValueError):
        periastron_eccentricity_cut(np.array([np.nan]), 1.0, ceiling=0.95)
    with np.testing.assert_raises(ValueError):
        periastron_eccentricity_cut(10.0, -1.0, ceiling=0.95)


def test_paired_period_samples_match_alderaan_polyfit() -> None:
    index = np.arange(12, dtype=float)
    model = 100.0 + 9.75 * index
    c0 = np.array([-0.02, 0.0, 0.03])
    c1 = np.array([-0.04, 0.01, 0.05])
    samples = pd.DataFrame({"C0_0": c0, "C1_0": c1})
    columns = [
        fits.Column(name="INDEX", format="D", array=index),
        fits.Column(name="MODEL", format="D", array=model),
    ]
    hdul = fits.HDUList(
        [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="TTIMES_00")]
    )
    got = paired_period_samples(samples, hdul, 0)
    centered = index - index[-1] // 2
    legx = centered / (index[-1] / 2.0)
    expected = []
    for row_c0, row_c1 in zip(c0, c1):
        ephem = model + row_c0 + row_c1 * legx
        expected.append(np.polynomial.polynomial.polyfit(index, ephem, 1)[1])
    np.testing.assert_allclose(got, expected, rtol=0.0, atol=1e-13)


def test_berger2018_density_source_requires_explicit_columns(tmp_path) -> None:
    from extract_eccentricity_posteriors_direct import load_density_source

    bad = tmp_path / "bad.csv"
    bad.write_text("kepid,rho_log\n1,0\n", encoding="utf-8")
    with np.testing.assert_raises(ValueError):
        load_density_source({"paths": {}}, "berger2018_kg", bad)


def test_berger2018_density_source_reads_requested_csv(tmp_path) -> None:
    from extract_eccentricity_posteriors_direct import load_density_source

    good = tmp_path / "density.csv"
    good.write_text("kepid,rho_log,rho_log_upper,rho_log_lower\n1,0.1,0.2,0.3\n", encoding="utf-8")
    got = load_density_source({"paths": {}}, "berger2018_kg", good)
    assert got.to_dict("records") == [{"kepid": 1, "rho_log": 0.1, "rho_log_upper": 0.2, "rho_log_lower": 0.3}]


def test_custom_density_source_reads_requested_csv(tmp_path) -> None:
    from extract_eccentricity_posteriors_direct import load_density_source

    good = tmp_path / "custom.csv"
    good.write_text(
        "kepid,rho_log,rho_log_upper,rho_log_lower\n1,-0.1,-0.8,-0.9\n",
        encoding="utf-8",
    )
    got = load_density_source({}, "custom_csv", good)
    assert got.to_dict("records") == [
        {"kepid": 1, "rho_log": -0.1, "rho_log_upper": -0.8, "rho_log_lower": -0.9}
    ]


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


def test_unweighted_resampling_is_deterministic_and_normalized() -> None:
    e_grid = np.linspace(0.0, 0.95, 40)
    omega_grid = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    kwargs = dict(
        eccentricity=np.array([0.01, 0.25, 0.8]),
        omega=np.array([0.1, 1.2, 4.8]),
        weights=np.array([0.1, 0.2, 0.7]),
        e_grid=e_grid,
        omega_grid=omega_grid,
        n_draws=20_000,
    )
    first = resampled_posterior_grid(rng=np.random.default_rng(12), **kwargs)
    second = resampled_posterior_grid(rng=np.random.default_rng(12), **kwargs)
    np.testing.assert_array_equal(first[0], second[0])
    np.testing.assert_array_equal(first[1], second[1])
    np.testing.assert_array_equal(first[2], second[2])
    assert np.isclose(first[0].sum(), 1.0)
    assert np.isclose(np.mean(first[1] == 0.8), 0.7, atol=0.02)


def synthetic_result(
    seed: int = 991,
    impact_center: float = 0.2,
    density_sampling_mode: str = "fixed_central",
) -> dict[str, object]:
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
        density_sampling_mode=density_sampling_mode,
        seed=seed,
    )


def test_deterministic_output() -> None:
    first = synthetic_result()
    second = synthetic_result()
    assert np.array_equal(first["posterior"], second["posterior"])
    assert np.array_equal(first["e_quantiles"], second["e_quantiles"])
    assert first["importance_ess"] == second["importance_ess"]


def test_gilbert_two_stage_sampling_is_deterministic_and_uses_requested_size() -> None:
    kwargs = dict(
        ror=np.full(12, 0.05),
        impact=np.linspace(0.1, 0.4, 12),
        dur14_days=np.full(12, 0.18),
        nested_weights=np.full(12, 1.0 / 12.0),
        period_days=10.0,
        rho_true=1.0,
        err_hi=0.1,
        err_lo=0.1,
        n_proposals=17,
        e_max=0.95,
        e_grid=np.linspace(0.0, 0.95, 30),
        omega_grid=np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False),
        density_error_mode="symmetric-average",
        seed=123,
        proposal_sampling_mode="gilbert_two_stage",
        transit_resample_draws=9,
        proposals_per_transit_draw=7,
    )
    first = direct_importance_posterior(**kwargs)
    second = direct_importance_posterior(**kwargs)
    assert first["proposal_count"] == 63
    assert first["proposal_sampling_mode"] == "gilbert_two_stage"
    assert first["transit_resample_draws"] == 9
    assert first["proposals_per_transit_draw"] == 7
    np.testing.assert_array_equal(first["posterior"], second["posterior"])


def test_gilbert_inverse_jacobian_is_explicit_and_deterministic() -> None:
    base = synthetic_result(seed=20260810)
    kwargs = dict(
        ror=np.full(250, 0.05),
        impact=np.linspace(0.175, 0.225, 250),
        dur14_days=np.full(250, 0.2),
        nested_weights=np.full(250, 1.0 / 250.0),
        period_days=10.0,
        rho_true=1.0,
        err_hi=0.1,
        err_lo=0.1,
        n_proposals=20_000,
        e_max=0.95,
        e_grid=np.linspace(0.0, 0.95, 80),
        omega_grid=np.linspace(0.0, 2.0 * np.pi, 72, endpoint=False),
        density_error_mode="symmetric-average",
        jacobian_mode="gilbert_inverse",
        seed=20260810,
    )
    first = direct_importance_posterior(**kwargs)
    second = direct_importance_posterior(**kwargs)
    assert first["jacobian_mode"] == "gilbert_inverse"
    np.testing.assert_array_equal(first["posterior"], second["posterior"])
    assert np.isclose(np.asarray(first["posterior"]).sum(), 1.0)
    assert not np.array_equal(first["posterior"], base["posterior"])


def test_source_style_sir_uses_explicit_final_draw_count() -> None:
    kwargs = dict(
        ror=np.full(80, 0.05),
        impact=np.linspace(0.1, 0.3, 80),
        dur14_days=np.full(80, 0.2),
        nested_weights=np.full(80, 1.0 / 80.0),
        period_days=10.0,
        rho_true=1.0,
        err_hi=0.1,
        err_lo=0.1,
        n_proposals=10_000,
        e_max=0.95,
        e_grid=np.linspace(0.0, 0.95, 40),
        omega_grid=np.linspace(0.0, 2.0 * np.pi, 36, endpoint=False),
        density_error_mode="symmetric-rms",
        posterior_sampling_mode="unweighted_resample",
        posterior_resample_draws=1_000,
        seed=20260807,
    )
    first = direct_importance_posterior(**kwargs)
    second = direct_importance_posterior(**kwargs)
    assert first["posterior_draw_count"] == 1_000
    np.testing.assert_array_equal(first["posterior"], second["posterior"])
    np.testing.assert_array_equal(first["e_quantiles"], second["e_quantiles"])
    assert np.isclose(np.asarray(first["posterior"]).sum(), 1.0)


def test_sir_and_weighted_grid_share_proposals_for_the_same_seed() -> None:
    kwargs = dict(
        ror=np.full(80, 0.05),
        impact=np.linspace(0.1, 0.3, 80),
        dur14_days=np.full(80, 0.2),
        nested_weights=np.full(80, 1.0 / 80.0),
        period_days=10.0,
        rho_true=1.0,
        err_hi=0.1,
        err_lo=0.1,
        n_proposals=30_000,
        e_max=0.95,
        e_grid=np.linspace(0.0, 0.95, 40),
        omega_grid=np.linspace(0.0, 2.0 * np.pi, 36, endpoint=False),
        density_error_mode="symmetric-average",
        seed=20260810,
    )
    weighted = direct_importance_posterior(
        **kwargs,
        posterior_sampling_mode="weighted_grid",
    )
    sir = direct_importance_posterior(
        **kwargs,
        posterior_sampling_mode="unweighted_resample",
        posterior_resample_draws=30_000,
    )
    np.testing.assert_allclose(
        np.asarray(weighted["e_pdf"]),
        np.asarray(sir["e_pdf"]),
        atol=0.02,
    )


def test_sir_rejects_nonpositive_final_draw_count() -> None:
    with np.testing.assert_raises_regex(
        ValueError, "posterior_resample_draws must be positive"
    ):
        direct_importance_posterior(
            ror=np.full(20, 0.05),
            impact=np.full(20, 0.2),
            dur14_days=np.full(20, 0.2),
            nested_weights=np.full(20, 1.0 / 20.0),
            period_days=10.0,
            rho_true=1.0,
            err_hi=0.1,
            err_lo=0.1,
            n_proposals=1_000,
            e_max=0.95,
            e_grid=np.linspace(0.0, 0.95, 20),
            omega_grid=np.linspace(0.0, 2.0 * np.pi, 18, endpoint=False),
            density_error_mode="symmetric-rms",
            posterior_sampling_mode="unweighted_resample",
            posterior_resample_draws=0,
            seed=20260808,
        )


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


def test_symmetric_rms_density_likelihood_matches_gilbert_error_rule() -> None:
    from extract_eccentricity_posteriors_direct import density_log_likelihood

    rho = np.array([0.9, 1.2])
    sigma = np.sqrt(0.2**2 + 0.1**2) / np.sqrt(2.0)
    got = density_log_likelihood(rho, 1.0, 0.2, 0.1, "symmetric-rms")
    expected = -0.5 * ((rho - 1.0) / sigma) ** 2 - np.log(sigma)
    np.testing.assert_allclose(got, expected)


def test_nested_sample_weight_modes_are_explicit() -> None:
    log_weights = np.log(np.array([0.1, 0.2, 0.7]))
    log_likelihood = np.log(np.array([0.6, 0.3, 0.1]))
    np.testing.assert_allclose(nested_sample_weights(log_weights, "dynesty"), [0.1, 0.2, 0.7])
    np.testing.assert_allclose(nested_sample_weights(log_weights, "equal"), np.full(3, 1.0 / 3.0))
    np.testing.assert_allclose(
        nested_sample_weights(log_weights, "likelihood", log_likelihood),
        [0.6, 0.3, 0.1],
    )


def test_likelihood_nested_weight_mode_requires_likelihood_column() -> None:
    with np.testing.assert_raises_regex(ValueError, "LN_LIKE"):
        nested_sample_weights(np.zeros(3), "likelihood")


def test_density_sensitivity_transform_preserves_fractional_errors() -> None:
    rho, err_hi, err_lo = adjust_stellar_density(
        1.0,
        0.2,
        0.1,
        offset_dex=np.log10(2.0),
        error_scale=1.5,
    )
    np.testing.assert_allclose([rho, err_hi, err_lo], [2.0, 0.6, 0.3])


def test_density_draw_mode_is_reproducible_and_finite() -> None:
    kwargs = dict(
        ror=np.full(120, 0.05),
        impact=np.linspace(0.15, 0.25, 120),
        dur14_days=np.full(120, 0.22),
        nested_weights=np.full(120, 1.0 / 120.0),
        period_days=10.0,
        rho_true=1.0,
        err_hi=0.1,
        err_lo=0.1,
        n_proposals=10_000,
        e_max=0.95,
        e_grid=np.linspace(0.0, 0.95, 40),
        omega_grid=np.linspace(0.0, 2.0 * np.pi, 36, endpoint=False),
        density_error_mode="symmetric-average",
        density_sampling_mode="draw_gaussian",
        seed=20260723,
    )
    first = direct_importance_posterior(**kwargs)
    second = direct_importance_posterior(**kwargs)
    np.testing.assert_array_equal(first["posterior"], second["posterior"])
    assert np.all(np.isfinite(first["posterior"]))
    assert np.all(np.isfinite(first["e_quantiles"]))


def test_split_normal_density_draw_is_reproducible_and_finite() -> None:
    kwargs = dict(
        ror=np.full(120, 0.05),
        impact=np.linspace(0.15, 0.25, 120),
        dur14_days=np.full(120, 0.22),
        nested_weights=np.full(120, 1.0 / 120.0),
        period_days=10.0,
        rho_true=1.0,
        err_hi=0.2,
        err_lo=0.1,
        n_proposals=10_000,
        e_max=0.95,
        e_grid=np.linspace(0.0, 0.95, 40),
        omega_grid=np.linspace(0.0, 2.0 * np.pi, 36, endpoint=False),
        density_error_mode="symmetric-average",
        density_sampling_mode="draw_split_normal",
        seed=20260724,
    )
    first = direct_importance_posterior(**kwargs)
    second = direct_importance_posterior(**kwargs)
    np.testing.assert_array_equal(first["posterior"], second["posterior"])
    assert np.all(np.isfinite(first["posterior"]))
    assert np.all(np.isfinite(first["e_quantiles"]))


def test_density_sampling_mode_rejects_unknown_value() -> None:
    with np.testing.assert_raises(ValueError):
        synthetic_result(density_sampling_mode="unknown")
