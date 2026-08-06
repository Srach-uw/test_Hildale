from __future__ import annotations

import numpy as np

from hierarchical_table3_order_diagnostic import population_normalizer


def test_forward_population_normalizer_matches_definition() -> None:
    e = np.linspace(0.0, 0.9, 1000)
    density = np.full_like(e, 1.0 / 0.9)
    expected = np.trapezoid(density / (1.0 - e**2), e)
    assert np.isclose(
        population_normalizer(density, e, "legacy_forward_norm"),
        expected,
    )
    assert expected > 1.0


def test_non_normalized_sensitivity_modes_return_unity() -> None:
    e = np.linspace(0.0, 0.9, 10)
    density = np.ones_like(e)
    assert population_normalizer(density, e, "none") == 1.0
    assert (
        population_normalizer(
            density,
            e,
            "manuscript_reciprocal",
        )
        == 1.0
    )
