from __future__ import annotations

import numpy as np

from missing_coverage_extreme_bound import synthetic_near_circular_mass


def test_near_circular_mass_is_finite_and_low_e() -> None:
    e_grid = np.linspace(0.0, 0.95, 300)
    omega_grid = np.linspace(-np.pi, np.pi, 240, endpoint=False)
    mass = synthetic_near_circular_mass(
        e_grid,
        omega_grid,
        scale=0.005,
        selection_mode="manuscript_reciprocal",
    )
    assert mass.shape == e_grid.shape
    assert np.all(np.isfinite(mass))
    assert np.all(mass >= 0)
    assert np.argmax(mass) <= 2
