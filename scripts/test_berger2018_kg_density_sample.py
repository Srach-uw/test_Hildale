from pathlib import Path

import numpy as np
import pytest
import pandas as pd

from build_berger2018_kg_density_sample import (
    LOGG_SUN,
    density_with_asymmetric_errors,
    read_kg,
)
from extract_eccentricity_posteriors import absolute_density_error


def test_official_kg_radii_catalog_has_required_columns_and_host_coverage():
    path = Path(__file__).resolve().parents[1] / "data" / "berger2018_kg_radii_star_cat.fits"
    if not path.exists():
        pytest.skip("optional Berger 2018 radii catalog is not distributed with the repository")
    frame = read_kg(path)
    assert len(frame) == 177_911
    assert {"kepid", "berger2018_kg_logg", "berger2018_kg_radius"}.issubset(frame.columns)
    row = frame.iloc[0]
    rho = 10 ** (float(row["berger2018_kg_logg"]) - LOGG_SUN) / float(row["berger2018_kg_radius"])
    assert np.isfinite(rho) and rho > 0


def test_kg_density_identity_is_dimensionally_consistent():
    logg = 4.438
    radius = 1.0
    rho = 10 ** (logg - LOGG_SUN) / radius
    assert np.isclose(rho, 1.0)


def test_kg_uncertainties_use_extractor_absolute_error_contract():
    rho, err_hi, err_lo = density_with_asymmetric_errors(
        pd.Series([LOGG_SUN]),
        pd.Series([0.1]),
        pd.Series([1.0]),
        pd.Series([0.05]),
        pd.Series([0.04]),
    )
    assert np.isclose(rho.iloc[0], 1.0)
    assert err_hi.iloc[0] > 0
    assert err_lo.iloc[0] > 0
    stored_hi = np.log10(err_hi.iloc[0])
    stored_lo = np.log10(err_lo.iloc[0])
    assert np.isclose(absolute_density_error(stored_hi), err_hi.iloc[0])
    assert np.isclose(absolute_density_error(stored_lo), err_lo.iloc[0])
