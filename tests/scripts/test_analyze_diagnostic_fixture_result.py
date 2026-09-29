from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_diagnostic_fixture_result import DiagnosticResultError, analyze_result  # noqa: E402
from diagnostic_transit_injection import circular_first_to_fourth_duration_days  # noqa: E402


def write_manifest(
    path: Path,
    *,
    period_days: float = 12.0,
    rho_star_solar: float = 1.1,
    prepared_fixture: bool = False,
) -> None:
    specification = {
        "period_days": period_days,
        "rho_star_solar": rho_star_solar,
        "impact_parameter": 0.2,
        "radius_ratio": 0.03,
        "seed": 7,
    }
    payload = {"injection": {"specification": specification}} if prepared_fixture else {
        "specification": specification
    }
    path.write_text(
        json.dumps(payload),
        encoding="utf-8",
    )


def write_result(path: Path, *, include_duration: bool = True) -> None:
    period_days = 12.0
    rho_star_solar = 1.1
    ror = np.full(9, 0.03)
    impact = np.full(9, 0.2)
    duration = circular_first_to_fourth_duration_days(
        period_days, rho_star_solar, ror[0], impact[0]
    )
    columns = [fits.Column(name="LN_WT", format="D", array=np.linspace(-8.0, 0.0, 9))]
    if include_duration:
        columns.append(fits.Column(name="DUR14_0", format="D", array=np.full(9, duration)))
    columns.extend(
        [
            fits.Column(name="ROR_0", format="D", array=ror),
            fits.Column(name="IMPACT_0", format="D", array=impact),
        ]
    )
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="SAMPLES")]).writeto(path)


def test_analyzer_recovers_declared_circular_density(tmp_path: Path) -> None:
    manifest = tmp_path / "injection.json"
    result = tmp_path / "result.fits"
    write_manifest(manifest)
    write_result(result)

    summary = analyze_result(result, manifest)

    assert summary["recovered_rho_circular_solar"]["p50"] == pytest.approx(1.1, rel=1e-10)
    assert summary["injected_within_recovered_68pct_interval"]
    assert summary["injected_within_recovered_95pct_interval"]
    assert not summary["median_bias_exceeds_warning_threshold"]
    assert summary["density_effective_sample_size"] > 1
    assert summary["bias_decomposition"]["injected_duration_days"] > 0
    assert summary["bias_decomposition"]["counterfactual_circular_density_solar"]["paired_posterior"]["p50"] == pytest.approx(1.1, rel=1e-10)


def test_analyzer_reads_prepared_fixture_manifest(tmp_path: Path) -> None:
    manifest = tmp_path / "fixture_manifest.json"
    result = tmp_path / "result.fits"
    write_manifest(manifest, prepared_fixture=True)
    write_result(result)

    summary = analyze_result(result, manifest)

    assert summary["injected"]["rho_star_solar"] == pytest.approx(1.1)


def test_analyzer_rejects_missing_transit_shape_column(tmp_path: Path) -> None:
    manifest = tmp_path / "injection.json"
    result = tmp_path / "result.fits"
    write_manifest(manifest)
    write_result(result, include_duration=False)

    with pytest.raises(DiagnosticResultError, match="DUR14_0"):
        analyze_result(result, manifest)
