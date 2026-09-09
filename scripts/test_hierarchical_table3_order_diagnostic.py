from hierarchical_table3_order_diagnostic import (
    PUBLISHED_BY_MODEL,
    PUBLISHED_ORDER,
    beta_density,
    half_gaussian_density,
)

import numpy as np

from common import trapezoid


def test_each_model_uses_its_own_literal_table3_values() -> None:
    assert tuple(PUBLISHED_BY_MODEL) == ("beta", "monotonic_beta", "half_gaussian")
    assert PUBLISHED_ORDER == (
        "thick_singles",
        "thick_multis",
        "thin_singles",
        "thin_multis",
    )
    assert PUBLISHED_BY_MODEL["beta"]["thick_singles"] == 0.058
    assert PUBLISHED_BY_MODEL["monotonic_beta"]["thick_singles"] == 0.041
    assert PUBLISHED_BY_MODEL["half_gaussian"]["thick_multis"] == 0.037


def test_outlier_floor_adds_tail_mass_and_preserves_normalization() -> None:
    e = np.linspace(0.0, 0.95, 200)
    for density, args in ((beta_density, (0.5, 8.0)), (half_gaussian_density, (0.05,))):
        plain = density(e, *args)
        floored = density(e, *args, outlier_floor=1e-6)
        assert np.isclose(trapezoid(plain, e), 1.0)
        assert np.isclose(trapezoid(floored, e), 1.0)
        assert floored[-1] > plain[-1]
