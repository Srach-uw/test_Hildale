from pathlib import Path


def test_repository_collector_requires_stopped_workers_and_uses_repo_paths() -> None:
    source = (
        Path(__file__).resolve().parent.parent
        / "cloud"
        / "collect_sharded_recovery.ps1"
    ).read_text(encoding="utf-8")
    assert 'if ($state -ne "TERMINATED")' in source
    assert "Refusing to collect shard" in source
    assert '$merger = Join-Path $projectRoot "scripts\\merge_sharded_recovery_archives.py"' in source
    assert '$targets = Join-Path $PSScriptRoot "recovery_142\\targets_missing_launchable.csv"' in source


def test_collector_has_no_stale_script_directory_merge_paths() -> None:
    source = (
        Path(__file__).resolve().parent.parent
        / "cloud"
        / "collect_sharded_recovery.ps1"
    ).read_text(encoding="utf-8")
    assert 'Join-Path $repo "merge_sharded_recovery_archives.py"' not in source
    assert "cloud_published_inventory_missing_batch\\targets_missing_launchable.csv" not in source
