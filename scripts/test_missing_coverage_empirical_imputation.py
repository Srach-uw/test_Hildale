from __future__ import annotations

import numpy as np

from hierarchical_rayleigh import posterior_weights_from_ll, weighted_quantile


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
