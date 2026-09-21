"""Check density distributions against their integral definitions."""

import numpy as np
from scipy.integrate import quad

from extract_eccentricity_posteriors_direct import density_log_likelihood


def test_continuous_split_density_has_unit_integral_and_expected_side_mass():
    def pdf(x):
        return float(np.exp(density_log_likelihood(np.array(x), 1.0, 0.2, 0.1, "split-continuous")))

    lower = quad(pdf, -np.inf, 1.0)[0]
    upper = quad(pdf, 1.0, np.inf)[0]
    np.testing.assert_allclose([lower, upper], [1.0 / 3.0, 2.0 / 3.0], atol=1e-10)
    np.testing.assert_allclose(lower + upper, 1.0, atol=1e-10)


def test_continuous_split_density_has_no_jump_at_the_mode():
    values = density_log_likelihood(np.array([1.0 - 1e-9, 1.0 + 1e-9]), 1.0, 0.2, 0.1, "split-continuous")
    np.testing.assert_allclose(values[0], values[1], atol=1e-12)


def test_continuous_split_density_reduces_to_gaussian():
    x = np.linspace(0.5, 1.5, 101)
    actual = density_log_likelihood(x, 1.0, 0.1, 0.1, "split-continuous")
    expected = -0.5 * ((x - 1.0) / 0.1) ** 2 - np.log(0.1 * np.sqrt(2.0 * np.pi))
    np.testing.assert_allclose(actual, expected)


def test_legacy_split_remains_explicitly_reproducible():
    x = np.array([1.0 - 1e-9, 1.0 + 1e-9])
    legacy = density_log_likelihood(x, 1.0, 0.2, 0.1, "split")
    np.testing.assert_allclose(np.exp(legacy[1] - legacy[0]), 0.5)
