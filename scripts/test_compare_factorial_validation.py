from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

import compare_factorial_validation as factorial

from compare_factorial_validation import (
    ARM_SPECS,
    ValidationAnalysisError,
    attach_repeatability_evidence,
    discover_validation_fits,
    read_target_sets,
    run_analysis,
    system_cluster_bootstrap,
)
from extract_eccentricity_posteriors_direct import DAY_S, G_SI, RHO_SUN_KG_M3


TARGET = "K00001"
KEPOI = "K00001.01"
KEPID = 1234567
PERIOD_DAYS = 10.0


def write_metadata(root: Path) -> tuple[Path, Path]:
    root.mkdir(parents=True, exist_ok=True)
    manifest = pd.DataFrame({"koi_target": [TARGET], "kepid": [KEPID]})
    for filename in (
        "targets_ld_reference_validation.csv",
        "targets_repeatability_validation.csv",
        "targets_short_cadence_validation.csv",
    ):
        manifest.to_csv(root / filename, index=False)
    inventory = pd.DataFrame(
        {
            "koi_target": [TARGET],
            "kepoi_name": [KEPOI],
            "kepid": [KEPID],
            "period_days": [PERIOD_DAYS],
            "included_in_alderaan_system": [True],
        }
    )
    inventory_path = root / "full_system_inventory.csv"
    inventory.to_csv(inventory_path, index=False)
    sample = pd.DataFrame(
        {
            "kepid": [KEPID],
            "rho_log": [0.0],
            "rho_log_upper": [0.08],
            "rho_log_lower": [-0.08],
            "disk": ["thin"],
            "system": ["single"],
        }
    )
    sample_path = root / "sample.csv"
    sample.to_csv(sample_path, index=False)
    return inventory_path, sample_path


def circular_duration_days(ror: float, impact: float) -> float:
    period_s = PERIOD_DAYS * DAY_S
    a_over_r = (G_SI * RHO_SUN_KG_M3 * period_s**2 / (3.0 * np.pi)) ** (1.0 / 3.0)
    argument = np.sqrt(((1.0 + ror) ** 2 - impact**2) / (a_over_r**2 - impact**2))
    return PERIOD_DAYS * np.arcsin(argument) / np.pi


def write_result(path: Path, duration_offset: float, ror_offset: float = 0.0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 80
    phase = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    ror = 0.05 + ror_offset + 0.0003 * np.sin(phase)
    impact = 0.2 + 0.01 * np.cos(phase)
    duration = np.array([circular_duration_days(r, b) for r, b in zip(ror, impact)])
    duration += duration_offset
    sample_columns = [
        fits.Column(name="LN_WT", format="D", array=np.zeros(count)),
        fits.Column(name="ROR_0", format="D", array=ror),
        fits.Column(name="IMPACT_0", format="D", array=impact),
        fits.Column(name="DUR14_0", format="D", array=duration),
    ]
    sample_hdu = fits.BinTableHDU.from_columns(sample_columns, name="SAMPLES")
    transit_columns = [
        fits.Column(name="INDEX", format="K", array=np.arange(5)),
        fits.Column(name="TTIME", format="D", array=np.arange(5) * PERIOD_DAYS),
        fits.Column(name="OUT_FLAG", format="K", array=np.zeros(5, dtype=int)),
    ]
    transit_hdu = fits.BinTableHDU.from_columns(transit_columns, name="TTIMES_00")
    primary = fits.PrimaryHDU()
    primary.header["NPL"] = 1
    fits.HDUList([primary, sample_hdu, transit_hdu]).writeto(path)


def write_all_arms(root: Path) -> None:
    offsets = {
        "original_lc": (0.0, 0.0),
        "reference_lc": (0.0010, 0.0004),
        "original_lc_repeat": (0.0001, 0.00005),
        "original_lcsc": (0.0005, 0.0002),
        "reference_lcsc": (0.0015, 0.0006),
        "paper_priors_original_lc": (0.0018, 0.0008),
    }
    for arm, spec in ARM_SPECS.items():
        result = root / "projects" / arm / "Results" / spec.run_id / TARGET / f"{TARGET}-results.fits"
        write_result(result, *offsets[arm])


def write_ld_pair(root: Path) -> None:
    offsets = {
        "original_lc": (0.0, 0.0),
        "reference_lc": (0.0010, 0.0004),
    }
    for arm, values in offsets.items():
        spec = ARM_SPECS[arm]
        result = root / "projects" / arm / "Results" / spec.run_id / TARGET / f"{TARGET}-results.fits"
        write_result(result, *values)


def test_discovery_is_recursive_and_arm_specific(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata"
    write_metadata(metadata)
    write_all_arms(tmp_path)
    target_sets = read_target_sets(metadata)
    discovery = discover_validation_fits(tmp_path, target_sets)
    assert len(discovery) == len(ARM_SPECS)
    assert discovery["status"].eq("present").all()
    assert set(discovery["arm"]) == set(ARM_SPECS)


def test_discovery_rejects_duplicate_arm_target(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata"
    write_metadata(metadata)
    write_all_arms(tmp_path)
    duplicate = tmp_path / "copy" / "original_lc" / f"{TARGET}-results.fits"
    write_result(duplicate, 0.0)
    with pytest.raises(ValidationAnalysisError, match="Multiple FITS files"):
        discover_validation_fits(tmp_path, read_target_sets(metadata))


def test_system_cluster_bootstrap_is_deterministic_and_clustered() -> None:
    frame = pd.DataFrame(
        {
            "koi_target": ["K1", "K1", "K2"],
            "delta": [0.0, 10.0, 100.0],
        }
    )
    first = system_cluster_bootstrap(frame, "delta", np.median, n_bootstrap=200, seed=42)
    second = system_cluster_bootstrap(frame, "delta", np.median, n_bootstrap=200, seed=42)
    assert np.array_equal(first, second)
    assert set(np.unique(first)).issubset({5.0, 10.0, 100.0})
    assert {5.0, 10.0, 100.0}.issubset(set(first))


def test_practical_evidence_uses_repeatability_distribution() -> None:
    metrics = pd.DataFrame(
        {
            "comparison_id": ["sampler_repeatability", "sampler_repeatability", "ld_reference_lc"],
            "parameter": ["e", "e", "e"],
            "abs_delta": [0.01, 0.02, 0.05],
        }
    )
    labeled, thresholds = attach_repeatability_evidence(metrics)
    e_threshold = thresholds.loc[thresholds["parameter"].eq("e"), "repeatability_abs_delta_p95"].iloc[0]
    assert np.isclose(e_threshold, 0.0195)
    assert labeled.iloc[2]["practical_evidence"] == "exceeds_observed_repeatability_p95_exploratory"
    assert labeled.iloc[0]["practical_evidence"] == "repeatability_reference_distribution"


def test_missing_fits_fails_after_writing_discovery_audit(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata"
    inventory, sample = write_metadata(metadata)
    output = tmp_path / "output"
    with pytest.raises(ValidationAnalysisError, match="FITS are incomplete"):
        run_analysis(
            validation_root=tmp_path,
            metadata_root=metadata,
            sample_path=sample,
            inventory_path=inventory,
            output_dir=output,
            n_proposals=100,
            e_max=0.95,
            density_error_mode="symmetric-average",
            period_tol=0.01,
            min_importance_ess=1.0,
            e_grid_size=10,
            omega_grid_size=12,
            n_bootstrap=20,
            seed=7,
        )
    audit = pd.read_csv(output / "factorial_validation_discovery.csv")
    assert len(audit) == len(ARM_SPECS)
    assert audit["status"].eq("missing").all()


def test_partial_snapshot_analyzes_only_complete_pairs(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata"
    inventory, sample = write_metadata(metadata)
    write_ld_pair(tmp_path)
    output = tmp_path / "output"

    paths = run_analysis(
        validation_root=tmp_path,
        metadata_root=metadata,
        sample_path=sample,
        inventory_path=inventory,
        output_dir=output,
        n_proposals=4_000,
        e_max=0.95,
        density_error_mode="symmetric-average",
        period_tol=0.01,
        min_importance_ess=1.0,
        e_grid_size=30,
        omega_grid_size=24,
        n_bootstrap=100,
        seed=7,
        allow_incomplete=True,
    )

    discovery = pd.read_csv(paths["discovery"])
    metrics = pd.read_csv(paths["paired_metrics"])
    thresholds = pd.read_csv(paths["repeatability"])
    report = Path(paths["report"]).read_text(encoding="utf-8")
    assert discovery["status"].eq("present").sum() == 2
    assert set(metrics["comparison_id"]) == {"ld_reference_lc"}
    assert thresholds["repeatability_abs_delta_p95"].isna().all()
    assert metrics["practical_evidence"].eq("not_assessable_without_repeatability").all()
    assert "Analysis mode: partial snapshot" in report


def test_direct_exclusions_fail_clearly_and_are_written(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = tmp_path / "metadata"
    inventory, sample = write_metadata(metadata)
    write_all_arms(tmp_path)

    def excluded_target(*args: object, **kwargs: object) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
        results_file = Path(args[0])
        planets = args[1]
        return [], [
            {
                "koi_target": TARGET,
                "kepid": KEPID,
                "kepoi_name": KEPOI,
                "koi_period": PERIOD_DAYS,
                "results_file": str(results_file),
                "stage": "importance",
                "reason": "synthetic direct failure",
            }
            for _ in range(len(planets))
        ]

    monkeypatch.setattr(factorial, "extract_direct_target", excluded_target)
    output = tmp_path / "output"
    with pytest.raises(ValidationAnalysisError, match="Direct eccentricity outputs could not be produced"):
        run_analysis(
            validation_root=tmp_path,
            metadata_root=metadata,
            sample_path=sample,
            inventory_path=inventory,
            output_dir=output,
            n_proposals=100,
            e_max=0.95,
            density_error_mode="symmetric-average",
            period_tol=0.01,
            min_importance_ess=1.0,
            e_grid_size=10,
            omega_grid_size=12,
            n_bootstrap=20,
            seed=7,
        )
    exclusions = pd.read_csv(output / "factorial_validation_direct_exclusions.csv")
    assert len(exclusions) == len(ARM_SPECS)
    assert exclusions["reason"].eq("synthetic direct failure").all()


def test_end_to_end_synthetic_analysis_is_deterministic(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata"
    inventory, sample = write_metadata(metadata)
    write_all_arms(tmp_path)

    common = dict(
        validation_root=tmp_path,
        metadata_root=metadata,
        sample_path=sample,
        inventory_path=inventory,
        n_proposals=4_000,
        e_max=0.95,
        density_error_mode="symmetric-average",
        period_tol=0.01,
        min_importance_ess=1.0,
        e_grid_size=30,
        omega_grid_size=24,
        n_bootstrap=200,
        seed=1234,
    )
    first_paths = run_analysis(output_dir=tmp_path / "out_first", **common)
    second_paths = run_analysis(output_dir=tmp_path / "out_second", **common)

    arm_planets = pd.read_csv(first_paths["arm_planets"])
    metrics = pd.read_csv(first_paths["paired_metrics"])
    summaries_first = pd.read_csv(first_paths["arm_summaries"])
    summaries_second = pd.read_csv(second_paths["arm_summaries"])
    assert len(arm_planets) == len(ARM_SPECS)
    assert len(metrics) == 5 * 5
    assert len(summaries_first) == 5 * 5
    assert np.isfinite(
        arm_planets[["e50", "zeta_median", "t14_hr_p50", "impact_p50", "rp_over_rs_p50"]]
    ).all().all()
    pd.testing.assert_frame_equal(summaries_first, summaries_second)
    assert Path(first_paths["report"]).read_text(encoding="utf-8").find("not p-values") >= 0
