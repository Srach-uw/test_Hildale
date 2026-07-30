from pathlib import Path

import pandas as pd
import pytest

from postprocess_published_inventory_recovery import (
    coverage_status,
    find_results_dir,
    project_root_for_script,
    safe_extract,
)


def test_find_results_dir_requires_one_results_root(tmp_path: Path) -> None:
    result = tmp_path / "Results" / "run" / "K00001" / "K00001-results.fits"
    result.parent.mkdir(parents=True)
    result.touch()
    root, count = find_results_dir(tmp_path)
    assert root == tmp_path / "Results"
    assert count == 1


def test_coverage_gate_counts_exact_missing(tmp_path: Path) -> None:
    sample = pd.DataFrame(
        [
            {"kepoi_name": "K00001.01", "disk": "thin", "system": "single"},
            {"kepoi_name": "K00002.01", "disk": "thick", "system": "multi"},
        ]
    )
    summary = pd.DataFrame([{"kepoi_name": "K00001.01"}])
    queue = pd.DataFrame([{"koi_target": "K00002"}])
    sample_path = tmp_path / "sample.csv"
    summary_path = tmp_path / "summary.csv"
    queue_path = tmp_path / "queue.csv"
    sample.to_csv(sample_path, index=False)
    summary.to_csv(summary_path, index=False)
    queue.to_csv(queue_path, index=False)
    table, missing = coverage_status(sample_path, summary_path, queue_path)
    assert missing == 1
    assert int(table["missing_planets"].sum()) == 1


def test_safe_extract_rejects_path_traversal(tmp_path: Path) -> None:
    import io
    import tarfile

    archive = tmp_path / "bad.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        member = tarfile.TarInfo("../escape.txt")
        payload = b"bad"
        member.size = len(payload)
        handle.addfile(member, io.BytesIO(payload))
    with pytest.raises(ValueError, match="escapes destination"):
        safe_extract(archive, tmp_path / "out")


def test_failure_taxonomy_uses_discovered_results_dir() -> None:
    source = Path(__file__).with_name("postprocess_published_inventory_recovery.py")
    text = source.read_text(encoding="utf-8")
    assert "str(results_dir)," in text
    assert "str(results_root)," not in text


def test_project_root_supports_repository_and_root_level_layouts(tmp_path: Path) -> None:
    repository_script = tmp_path / "repo" / "scripts" / "postprocess.py"
    repository_script.parent.mkdir(parents=True)
    repository_script.touch()
    assert project_root_for_script(repository_script) == tmp_path / "repo"

    root_level_script = tmp_path / "research" / "postprocess.py"
    root_level_script.parent.mkdir(parents=True)
    root_level_script.touch()
    assert project_root_for_script(root_level_script) == tmp_path / "research"
