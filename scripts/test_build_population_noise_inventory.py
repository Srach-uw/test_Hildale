import os
from pathlib import Path

import pandas as pd
import pytest
from astropy.io import fits

from build_population_noise_inventory import build_inventory, sha256


POPULATIONS = [
    ("thin", "single"),
    ("thick", "single"),
    ("thin", "multi"),
    ("thick", "multi"),
]
SOURCE_ROOT = "inputs/alderaan_archive_raw_20260612/ALDERAAN_posteriors"


def _result(path, target, npl, marker="selected"):
    header = fits.Header({"TARGET": target, "NPL": npl, "MARKER": marker})
    fits.PrimaryHDU(header=header).writeto(path)


def _fixture(
    tmp_path,
    *,
    first_npl=1,
    first_required_npl=None,
    first_index=0,
    omit_second_timing=False,
):
    source_dir = tmp_path / "source"
    staging_dir = tmp_path / "staging"
    companion_root = tmp_path / "companions"
    leverage_dir = tmp_path / "leverage"
    output_dir = tmp_path / "output"
    for path in (source_dir, staging_dir, companion_root, leverage_dir):
        path.mkdir()

    summary_rows = []
    manifest_rows = []
    for index, (disk, system) in enumerate(POPULATIONS, start=1):
        target = f"K{index:05d}"
        planet = f"{target}.01"
        npl = first_npl if index == 1 else 1
        required_npl = (first_required_npl if first_required_npl is not None else npl) if index == 1 else 1
        source = source_dir / f"{target}-results.fits"
        staging = staging_dir / f"{target}-results.fits"
        _result(source, target, npl)
        os.link(source, staging)
        summary_rows.append(
            {"kepoi_name": planet, "koi_target": target, "alderaan_planet_index": first_index if index == 1 else 0, "disk": disk, "system": system, "qc_primary_exclude": False}
        )
        manifest_rows.append(
            {
                "koi_target": target,
                "selected_root": SOURCE_ROOT,
                "selected_file": source,
                "required_npl": required_npl,
                "selected_planet_count": required_npl,
                "selected_full_system_match": True,
                "selected_sha256": sha256(source),
                "candidate_count": 1,
                "candidate_files": source,
                "candidate_planet_counts": f"{source}:{npl}",
                "staged_file": staging,
                "staging_method": "hardlink",
            }
        )
        target_dir = companion_root / target
        target_dir.mkdir()
        for suffix in ("lc_filtered.fits", "transit_parameters.csv"):
            (target_dir / f"{target}_{suffix}").write_bytes(b"fixture")
        for planet_index in range(npl):
            if index == 1 and omit_second_timing and planet_index == 1:
                continue
            (target_dir / f"{target}_{planet_index:02d}_quick.ttvs").write_text("0 0 0\n", encoding="utf-8")

        pd.DataFrame(
            [{
                "kepoi_name": planet,
                "koi_target": target,
                "disk": disk,
                "system": system,
                "posterior_source": "alderaan_direct_importance",
                "rho_true_solar": 1.0,
                "e50": 0.1,
                "log_contrast_fit_over_paper": float(index),
            }]
        ).to_csv(leverage_dir / f"{disk}_{system}_leverage.csv", index=False)

    summary_path = tmp_path / "summary.csv"
    manifest_path = tmp_path / "manifest.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)
    return {
        "summary": summary_path,
        "manifest": manifest_path,
        "leverage": leverage_dir,
        "output": output_dir,
        "companions": [companion_root],
    }


def _build(paths, **kwargs):
    return build_inventory(
        paths["summary"],
        paths["manifest"],
        paths["leverage"],
        paths["output"],
        companion_roots=paths["companions"],
        expected_summary_sha256=sha256(paths["summary"]),
        expected_planets=4,
        expected_targets=4,
        expected_exclusions=set(),
        **kwargs,
    )


def test_builds_inventory_but_does_not_infer_replay_readiness(tmp_path):
    paths = _fixture(tmp_path)

    inventory, sources, pilots = _build(paths)

    assert len(inventory) == 4
    assert set(inventory["source_result_state"]) == {"verified"}
    assert set(inventory["lc_filtered_state"]) == {"present_unverified_source"}
    assert set(inventory["companion_source_match"]) == {"unverified_no_companion_hash_manifest"}
    assert set(inventory["replay_readiness"]) == {"blocked_unverified_companion_provenance"}
    assert int(sources.iloc[0]["verified_results"]) == 4
    assert len(pilots) == 4
    assert set(pilots["authorization_state"]) == {"candidate_not_authorized"}


def test_checks_timing_for_full_joint_system(tmp_path):
    paths = _fixture(tmp_path, first_npl=2, first_required_npl=1, omit_second_timing=True)

    inventory, _, _ = _build(paths)

    first = inventory.loc[inventory["koi_target"] == "K00001"].iloc[0]
    assert first["selected_planet_count_qc"] == 1
    assert first["manifest_required_planet_count"] == 1
    assert first["full_joint_planet_count"] == 2
    assert first["joint_fit_membership_state"] == "fits_has_additional_joint_planets"
    assert not first["all_full_system_timing_present"]
    assert str(first["missing_timing_indices"]) == "1"
    assert first["replay_readiness"] == "blocked_missing_companions"


def test_result_hash_mismatch_stops_and_reports_conflict(tmp_path):
    paths = _fixture(tmp_path)
    manifest = pd.read_csv(paths["manifest"])
    manifest.loc[0, "selected_sha256"] = "0" * 64
    manifest.to_csv(paths["manifest"], index=False)

    with pytest.raises(RuntimeError, match="conflicts"):
        _build(paths)

    conflicts = pd.read_csv(paths["output"] / "conflicts.csv")
    assert "source_hash_mismatch" in set(conflicts["conflict_type"])


def test_duplicate_source_manifest_target_is_rejected(tmp_path):
    paths = _fixture(tmp_path)
    manifest = pd.read_csv(paths["manifest"])
    manifest = pd.concat([manifest, manifest.iloc[[0]]], ignore_index=True)
    manifest.to_csv(paths["manifest"], index=False)

    with pytest.raises(ValueError, match="duplicate targets"):
        _build(paths)


def test_conflicting_full_system_alternative_stops(tmp_path):
    paths = _fixture(tmp_path)
    manifest = pd.read_csv(paths["manifest"])
    selected = Path(manifest.loc[0, "selected_file"])
    alternative = selected.with_name("K00001-alternative-results.fits")
    _result(alternative, "K00001", 1, marker="alternative")
    manifest.loc[0, "candidate_count"] = 2
    manifest.loc[0, "candidate_files"] = f"{selected}|{alternative}"
    manifest.loc[0, "candidate_planet_counts"] = f"{selected}:1|{alternative}:1"
    manifest.to_csv(paths["manifest"], index=False)

    with pytest.raises(RuntimeError, match="conflicts"):
        _build(paths)

    conflicts = pd.read_csv(paths["output"] / "conflicts.csv")
    assert "conflicting_full_system_fit" in set(conflicts["conflict_type"])


def test_incompatible_selected_planet_index_stops(tmp_path):
    paths = _fixture(tmp_path, first_npl=2, first_required_npl=1, first_index=2)

    with pytest.raises(RuntimeError, match="conflicts"):
        _build(paths)

    conflicts = pd.read_csv(paths["output"] / "conflicts.csv")
    assert "selected_planet_index_incompatible" in set(conflicts["conflict_type"])
