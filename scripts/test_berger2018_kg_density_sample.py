from pathlib import Path

import numpy as np
import pytest

from build_berger2018_kg_density_sample import LOGG_SUN, read_kg


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
