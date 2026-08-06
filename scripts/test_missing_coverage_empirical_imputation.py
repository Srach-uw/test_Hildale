from __future__ import annotations

import numpy as np

from hierarchical_rayleigh import posterior_weights_from_ll, weighted_quantile
from missing_coverage_empirical_imputation import population_log_terms


def test_duplicate_log_likelihood_moves_toward_analog() -> None:
    sigmas = np.linspace(0.001, 0.5, 500)
    expected = sigmas * np.sqrt(np.pi / 2.0)
    full_ll = -0.5 * np.square((sigmas - 0.2) / 0.03)
    low_e_analog = -0.5 * np.square((sigmas - 0.05) / 0.02)
    baseline = weighted_quantile(
        expected, posterior_weights_from_ll(sigmas, full_ll), [0.5]
    )[0]
    imputed = weighted_quantile(
        expected,
        posterior_weights_from_ll(sigmas, full_ll + 10 * low_e_analog),
        [0.5],
    )[0]
    assert imputed < baseline


def test_forward_mode_divides_population_normalizer() -> None:
    masses = np.array([[0.25, 0.75]])
    densities = np.array([[1.0, 2.0], [3.0, 4.0]])
    normalizers = np.array([2.0, 5.0])
    expected = np.log((masses @ densities) / normalizers)
    actual = population_log_terms(
        masses,
        densities,
        normalizers,
        "legacy_forward_norm",
    )
    assert np.allclose(actual, expected)


def test_literal_reciprocal_mode_does_not_divide_normalizer() -> None:
    masses = np.array([[0.25, 0.75]])
    densities = np.array([[1.0], [3.0]])
    actual = population_log_terms(
        masses,
        densities,
        np.array([99.0]),
        "manuscript_reciprocal",
    )
    assert np.allclose(actual, np.log(masses @ densities))
