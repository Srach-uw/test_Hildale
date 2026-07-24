from __future__ import annotations

import tarfile
from pathlib import Path

import pandas as pd
import pytest

from merge_sharded_recovery_archives import merge_archives, safe_extract, write_tar


def make_archive(root: Path, name: str, targets: list[str]) -> Path:
    payload = root / f"{name}_payload"
    for target in targets:
        result = (
            payload
            / "Results"
            / "sagear_published_inventory_missing"
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
