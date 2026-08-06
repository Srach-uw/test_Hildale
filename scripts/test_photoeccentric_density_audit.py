from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from photoeccentric_density_audit import (
    add_impact_widths,
    fit_impact_constrained_populations,
    normalized_log_weights,
    merge_qc_audit,
    prepare_planets,
    summarize,
    summarize_impact_width,
    sha256_file,
    weighted_quantile,
)


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


def test_log_weights_and_weighted_quantiles_are_stable() -> None:
    weights = normalized_log_weights(np.array([-1002.0, -1001.0, -1000.0]))
    assert np.isclose(weights.sum(), 1.0)
    assert np.all(np.isfinite(weights))
    quantile = weighted_quantile(
        np.array([0.1, 0.5, 0.9]),
        weights,
        [0.5],
    )[0]
    assert 0.5 < quantile < 0.9


def test_impact_summary_keeps_populations_and_width_classes() -> None:
    source = pd.DataFrame(
        {
            "kepoi_name": ["K1.01", "K2.01"],
            "disk": ["thin", "thin"],
            "system": ["single", "single"],
            "rho_circular_p16_solar": [0.5, 0.5],
            "rho_circular_median_solar": [1.0, 2.0],
            "rho_circular_p84_solar": [2.0, 4.0],
            "rho_catalog_solar": [1.0, 1.0],
            "e50": [0.1, 0.3],
        }
    )
    planets = prepare_planets(source)
    planets["impact_width"] = [0.2, 0.7]
    result = summarize_impact_width(planets, threshold=0.4)
    assert result["n"].sum() == 2
    assert set(result["impact_width_class"]) == {"constrained", "broad"}


def test_sha256_file(tmp_path) -> None:
    path = tmp_path / "input.csv"
    path.write_bytes(b"abc")
    assert (
        sha256_file(path)
        == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )


def test_merge_qc_audit_requires_the_same_unique_planets() -> None:
    summary = pd.DataFrame({"kepoi_name": ["K1.01"], "disk": ["thin"]})
    audit = pd.DataFrame(
        {"kepoi_name": ["K1.01"], "raw_result_file": ["K1-results.fits"]}
    )
    merged = merge_qc_audit(summary, audit)
    assert merged.loc[0, "raw_result_file"] == "K1-results.fits"

    duplicate = pd.concat([audit, audit], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        merge_qc_audit(summary, duplicate)


def test_add_impact_widths_uses_planet_index_and_nested_weights(tmp_path) -> None:
    result_file = tmp_path / "K00001-results.fits"
    samples = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="LN_WT", format="D", array=np.zeros(4)),
            fits.Column(
                name="IMPACT_0",
                format="D",
                array=np.array([0.1, 0.2, 0.3, 0.4]),
            ),
            fits.Column(
                name="IMPACT_1",
                format="D",
                array=np.array([0.5, 0.6, 0.7, 0.8]),
            ),
        ],
        name="SAMPLES",
    )
    fits.HDUList([fits.PrimaryHDU(), samples]).writeto(result_file)
    planets = pd.DataFrame(
        {
            "kepoi_name": ["K00001.01", "K00001.02"],
            "raw_result_file": [str(result_file), str(result_file)],
            "alderaan_planet_index": [0, 1],
        }
    )
    result = add_impact_widths(planets)
    assert result.loc[0, "impact_p50"] < result.loc[1, "impact_p50"]
    assert np.isclose(result.loc[0, "impact_width"], result.loc[1, "impact_width"])


def test_impact_constrained_fit_excludes_primary_qc_failures(tmp_path) -> None:
    e_grid = np.linspace(0.0, 0.95, 20)
    omega_grid = np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False)
    posterior = np.ones((len(e_grid), len(omega_grid)))
    rows = []
    for index in range(6):
        posterior_path = tmp_path / f"posterior_{index}.npz"
        np.savez(
            posterior_path,
            e_grid=e_grid,
            omega_grid=omega_grid,
            posterior=posterior,
            include_transit_prior=np.array(False),
        )
        rows.append(
            {
                "kepoi_name": f"K00001.{index + 1:02d}",
                "disk": "thin",
                "system": "single",
                "impact_width": 0.2,
                "qc_primary_exclude": index == 0,
                "posterior_file": str(posterior_path),
            }
        )
    result = fit_impact_constrained_populations(pd.DataFrame(rows))
    assert set(result["selection_mode"]) == {
        "legacy_forward_norm",
        "manuscript_reciprocal",
    }
    assert result["n"].eq(5).all()
    assert result["n_primary_qc_excluded_before_fit"].eq(1).all()


def test_impact_constrained_fit_requires_qc_manifest() -> None:
    with pytest.raises(ValueError, match="qc_primary_exclude"):
        fit_impact_constrained_populations(
            pd.DataFrame(
                {
                    "impact_width": [0.2],
                    "disk": ["thin"],
                    "system": ["single"],
                }
            )
        )
