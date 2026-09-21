"""Verify the simulation's selection law independently of hierarchy fitting."""

import numpy as np
from scipy.integrate import quad

from population_recovery_control import density_observable, draw_transiting_population


def test_circular_density_observable_is_zero():
    np.testing.assert_array_equal(density_observable(0.0, np.linspace(0, 2*np.pi, 90)), 0.0)


def test_transit_rejection_sampler_matches_analytic_selected_mean():
    mean = 0.2
    sigma = mean / np.sqrt(np.pi / 2)
    def selected(e):
        return e / sigma**2 * np.exp(-e**2 / (2*sigma**2)) / (1-e**2)
    norm = quad(selected, 0, 0.6)[0]
    expected = quad(lambda e: e*selected(e), 0, 0.6)[0] / norm
    e, w = draw_transiting_population(np.random.default_rng(302), 100000, mean)
    assert np.all((e >= 0) & (e < 0.6))
    assert abs(e.mean() - expected) < 0.0015
    # Conditional transit selection gives E[sin(omega) | e] = e/2.
    assert abs(np.mean(np.sin(w) - e/2)) < 0.008
