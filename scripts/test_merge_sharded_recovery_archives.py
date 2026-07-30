from __future__ import annotations

import tarfile
from pathlib import Path

import pandas as pd
import pytest

from merge_sharded_recovery_archives import merge_archives, safe_extract, write_tar


def make_archive(
    root: Path,
    name: str,
    targets: list[str],
    run_id: str = "sagear_published_inventory_missing",
) -> Path:
    payload = root / f"{name}_payload"
    for target in targets:
        result = (
            payload
            / "Results"
            / run_id
            / target
            / f"{target}-results.fits"
        )
        result.parent.mkdir(parents=True, exist_ok=True)
        result.write_bytes((target + "\n").encode("ascii"))
    archive = root / f"{name}.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        handle.add(payload / "Results", arcname="Results")
    return archive


def test_merge_union_and_roundtrip_archive(tmp_path: Path) -> None:
    first = make_archive(tmp_path, "first", ["K00001"])
    second = make_archive(tmp_path, "second", ["K00002"])
    expected = pd.DataFrame({"koi_target": ["K00001", "K00002", "K00003"]})
    staging = tmp_path / "staging"
    manifest = merge_archives(
        [first, second],
        expected_targets=expected,
        staging=staging,
        run_id="sagear_published_inventory_missing",
    )
    assert manifest["merge_status"].tolist() == [
        "result_present",
        "result_present",
        "missing_result",
    ]
    merged = tmp_path / "merged.tar.gz"
    write_tar(staging, merged)
    extracted = tmp_path / "roundtrip"
    safe_extract(merged, extracted)
    assert len(list(extracted.rglob("*-results.fits"))) == 2


def test_duplicate_target_across_archives_is_fatal(tmp_path: Path) -> None:
    first = make_archive(tmp_path, "first", ["K00001"])
    second = make_archive(tmp_path, "second", ["K00001"])
    expected = pd.DataFrame({"koi_target": ["K00001"]})
    with pytest.raises(ValueError, match="duplicate target"):
        merge_archives(
            [first, second],
            expected_targets=expected,
            staging=tmp_path / "staging",
            run_id="sagear_published_inventory_missing",
        )


def test_target_table_must_be_unique_and_nonblank(tmp_path: Path) -> None:
    archive = make_archive(tmp_path, "one", ["K00001"])
    with pytest.raises(ValueError, match="duplicate koi_target"):
        merge_archives(
            [archive],
            expected_targets=pd.DataFrame({"koi_target": ["K00001", "K00001"]}),
            staging=tmp_path / "duplicates",
            run_id="sagear_published_inventory_missing",
        )
    with pytest.raises(ValueError, match="blank koi_target"):
        merge_archives(
            [archive],
            expected_targets=pd.DataFrame({"koi_target": ["K00001", " "]}),
            staging=tmp_path / "blank",
            run_id="sagear_published_inventory_missing",
        )


def test_result_from_a_different_run_is_fatal(tmp_path: Path) -> None:
    archive = make_archive(tmp_path, "wrong_run", ["K00001"], run_id="other_run")
    with pytest.raises(ValueError, match="not under Results/sagear_published_inventory_missing/K00001"):
        merge_archives(
            [archive],
            expected_targets=pd.DataFrame({"koi_target": ["K00001"]}),
            staging=tmp_path / "wrong_run",
            run_id="sagear_published_inventory_missing",
        )
