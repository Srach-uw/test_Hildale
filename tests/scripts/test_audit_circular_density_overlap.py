from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_circular_density_overlap import (  # noqa: E402
    CircularDensityOverlapError,
    read_paired_samples_from_results,
    summarize_overlap,
    write_audit,
)


def duration_days(period_days: float, rho_solar: float, ror: float, impact: float) -> float:
    period_s = period_days * 86400.0
    a_over_r = (6.67430e-11 * rho_solar * 1408.0 * period_s**2 / (3.0 * np.pi)) ** (1.0 / 3.0)
    argument = np.sqrt(((1.0 + ror) ** 2 - impact**2) / (a_over_r**2 - impact**2))
    return period_days * np.arcsin(argument) / np.pi


def paired_samples() -> pd.DataFrame:
    period = np.array([10.0, 10.0])
    ror = np.array([0.03, 0.03])
    impact = np.array([0.2, 0.2])
    return pd.DataFrame(
        {
            "P": period,
            "T14": [duration_days(10.0, 0.8, 0.03, 0.2), duration_days(10.0, 1.4, 0.03, 0.2)],
            "ROR": ror,
            "IMPACT": impact,
            "LN_WT": [-40.0, 0.0],
        }
    )


def test_dynesty_weighting_controls_circular_density_median() -> None:
    summary = summarize_overlap(
        paired_samples(), stellar_density_samples=np.full(100, 1.4), grid_size=64
    )
    assert summary["rho_circular_solar"]["p50"] == pytest.approx(1.4, rel=1e-8)
    assert summary["delta_log10_rho_circular_minus_stellar_dex"]["p50"] == pytest.approx(0.0, abs=1e-8)
    assert summary["dynesty"]["effective_sample_size"] == pytest.approx(1.0, rel=1e-8)


def test_outputs_are_deterministic_and_record_provenance(tmp_path: Path) -> None:
    samples = tmp_path / "paired.csv"
    stellar = tmp_path / "stellar.csv"
    paired_samples().to_csv(samples, index=False)
    pd.DataFrame({"rho_star_solar": [1.2, 1.3, 1.4, 1.5]}).to_csv(stellar, index=False)

    first_csv, first_json, first = write_audit(samples, tmp_path / "first", stellar_density_samples_path=stellar, grid_size=64)
    second_csv, second_json, second = write_audit(samples, tmp_path / "second", stellar_density_samples_path=stellar, grid_size=64)

    for key in ("configuration", "input", "dynesty", "rho_circular_solar", "rho_stellar_solar", "delta_log10_rho_circular_minus_stellar_dex"):
        assert first[key] == second[key]
    assert first_csv.read_text(encoding="utf-8") == second_csv.read_text(encoding="utf-8")
    manifest = json.loads(first_json.read_text(encoding="utf-8"))
    assert manifest["diagnostic_label"] == "diagnostic_only_not_a_significance_or_population_claim"
    assert len(manifest["provenance"]["samples_sha256"]) == 64
    assert len(manifest["provenance"]["stellar_samples_sha256"]) == 64
    assert manifest["configuration"]["grid_size"] == 64
    assert manifest["dynesty"]["weight_column"] == "LN_WT"
    assert second_json.exists()


def test_quantile_stellar_input_and_overlap_interpretation() -> None:
    summary = summarize_overlap(
        paired_samples(),
        stellar_density_quantiles={"p16": 1.2, "p50": 1.4, "p84": 1.6},
        tolerance_dex=0.05,
        grid_size=64,
    )
    delta = summary["delta_log10_rho_circular_minus_stellar_dex"]
    assert summary["input"]["stellar_density_source"] == "declared_p16_p50_p84_quantiles"
    assert delta["p16"] <= 0.0 <= delta["p84"]
    assert delta["zero_in_68pct_credible_interval"]
    assert 0.0 <= delta["probability_abs_delta_within_tolerance"] <= 1.0


@pytest.mark.parametrize(
    "frame,match",
    [
        (pd.DataFrame({"P": [10.0]}), "missing columns"),
        (pd.DataFrame({"P": [-1.0], "T14": [0.1], "ROR": [0.03], "IMPACT": [0.2], "LN_WT": [0.0]}), "P must"),
    ],
)
def test_invalid_paired_inputs_are_rejected(frame: pd.DataFrame, match: str) -> None:
    with pytest.raises(CircularDensityOverlapError, match=match):
        summarize_overlap(frame, stellar_density_samples=[1.0])


def test_invalid_stellar_declaration_is_rejected() -> None:
    with pytest.raises(CircularDensityOverlapError, match="exactly one"):
        summarize_overlap(
            paired_samples(),
            stellar_density_samples=[1.0],
            stellar_density_quantiles={"p16": 0.9, "p50": 1.0, "p84": 1.1},
        )
    with pytest.raises(CircularDensityOverlapError, match="p16 <= p50 <= p84"):
        summarize_overlap(
            paired_samples(), stellar_density_quantiles={"p16": 1.2, "p50": 1.0, "p84": 1.4}
        )


def make_result_fits(path: Path) -> None:
    samples = np.zeros(
        2,
        dtype=[
            ("DUR14_0", "f8"),
            ("ROR_0", "f8"),
            ("IMPACT_0", "f8"),
            ("LN_WT", "f8"),
            ("C0_0", "f8"),
            ("C1_0", "f8"),
        ],
    )
    samples["DUR14_0"] = [duration_days(10.0, 0.8, 0.03, 0.2), duration_days(10.0, 1.4, 0.03, 0.2)]
    samples["ROR_0"] = 0.03
    samples["IMPACT_0"] = 0.2
    samples["LN_WT"] = [-40.0, 0.0]
    tt = np.zeros(3, dtype=[("INDEX", "f8"), ("MODEL", "f8")])
    tt["INDEX"] = [0.0, 1.0, 2.0]
    tt["MODEL"] = [100.0, 110.0, 120.0]
    primary = fits.PrimaryHDU()
    primary.header["NPL"] = 1
    fits.HDUList(
        [primary, fits.BinTableHDU(samples, name="SAMPLES"), fits.BinTableHDU(tt, name="TTIMES_00")]
    ).writeto(path)


def test_results_fits_mode_preserves_paired_draws_and_hashes(tmp_path: Path) -> None:
    result_path = tmp_path / "K00001-results.fits"
    make_result_fits(result_path)
    paired = read_paired_samples_from_results(result_path, 0)
    assert paired.columns.tolist() == ["P", "T14", "ROR", "IMPACT", "LN_WT"]
    assert paired["P"].tolist() == pytest.approx([10.0, 10.0])
    csv_path, manifest_path, _ = write_audit(
        None,
        tmp_path / "result_mode",
        results_fits_path=result_path,
        planet_index=0,
        stellar_density_quantiles={"p16": 1.3, "p50": 1.4, "p84": 1.5},
        grid_size=64,
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert csv_path.exists()
    assert manifest["provenance"]["samples_path"] is None
    assert manifest["provenance"]["result_fits"]["planet_index"] == 0
    assert len(manifest["provenance"]["result_fits"]["results_fits_sha256"]) == 64


def test_results_fits_mode_requires_index_and_valid_index(tmp_path: Path) -> None:
    result_path = tmp_path / "K00001-results.fits"
    make_result_fits(result_path)
    with pytest.raises(CircularDensityOverlapError, match="planet_index is required"):
        write_audit(
            None,
            tmp_path / "missing_index",
            results_fits_path=result_path,
            stellar_density_quantiles={"p16": 0.9, "p50": 1.0, "p84": 1.1},
        )
    with pytest.raises(CircularDensityOverlapError, match="outside"):
        read_paired_samples_from_results(result_path, 1)
