import numpy as np
import pandas as pd

from gilbert_postfit_qc_audit import (
    catalog_density_solar,
    radius_fractional_uncertainty,
    weighted_fraction,
)


def test_weighted_fraction_uses_nested_weights() -> None:
    mask = np.array([False, True, True])
    weights = np.array([0.8, 0.1, 0.1])
    assert np.isclose(weighted_fraction(mask, weights), 0.2)


def test_radius_fractional_uncertainty_combines_independent_terms() -> None:
    ror = np.array([0.09, 0.10, 0.11])
    weights = np.ones(3) / 3
    _, median, _, frac = radius_fractional_uncertainty(
        ror,
        weights,
        stellar_radius=1.0,
        radius_err_hi=0.1,
        radius_err_lo=0.1,
    )
    assert np.isclose(median, 0.10)
    assert 0.10 < frac < 0.15


def test_catalog_density_prefers_direct_solar_value() -> None:
    row = pd.Series({"rho_true_solar": 0.7, "rho_log": -2.0})
    assert np.isclose(catalog_density_solar(row), 0.7)


def test_catalog_density_falls_back_to_log_value() -> None:
    row = pd.Series({"rho_true_solar": np.nan, "rho_log": -0.3})
    assert np.isclose(catalog_density_solar(row), 10**-0.3)
