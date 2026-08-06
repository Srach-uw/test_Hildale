from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from photoeccentric_density_audit import prepare_planets, summarize


def test_prepare_planets_computes_log_density_metrics() -> None:
    source = pd.DataFrame(
        {
            "kepoi_name": ["K00001.01"],
            "disk": ["thin"],
            "system": ["single"],
            "rho_circular_p16_solar": [0.5],
            "rho_circular_median_solar": [2.0],
            "rho_circular_p84_solar": [4.0],
            "rho_catalog_solar": [1.0],
            "e50": [0.2],
        }
    )
    result = prepare_planets(source).iloc[0]
    assert result["population"] == "thin_single"
    assert np.isclose(result["delta_log10_rho"], np.log10(2.0))
    assert np.isclose(result["rho_circular_width_dex"], np.log10(8.0))


def test_prepare_planets_rejects_invalid_quantile_order() -> None:
    source = pd.DataFrame(
        {
            "kepoi_name": ["K00001.01"],
            "disk": ["thin"],
            "system": ["single"],
            "rho_circular_p16_solar": [2.0],
            "rho_circular_median_solar": [1.0],
            "rho_circular_p84_solar": [4.0],
            "rho_catalog_solar": [1.0],
            "e50": [0.2],
        }
    )
    assert prepare_planets(source).empty


def test_summarize_reports_each_population() -> None:
    source = pd.DataFrame(
        {
            "kepoi_name": ["K00001.01", "K00002.01"],
            "disk": ["thin", "thick"],
            "system": ["single", "multi"],
            "rho_circular_p16_solar": [0.5, 0.8],
            "rho_circular_median_solar": [1.0, 1.2],
            "rho_circular_p84_solar": [2.0, 1.8],
            "rho_catalog_solar": [1.0, 1.0],
            "e50": [0.1, 0.2],
        }
    )
    result = summarize(prepare_planets(source))
    assert set(result["population"]) == {"thin_single", "thick_multi"}


def test_prepare_planets_requires_provenance_fields() -> None:
    with pytest.raises(ValueError, match="rho_catalog_solar"):
        prepare_planets(pd.DataFrame({"kepoi_name": ["K00001.01"]}))
