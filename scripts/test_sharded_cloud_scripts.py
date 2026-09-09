from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "cloud"


def test_launch_uses_stop_and_distinct_regions() -> None:
    script = (ROOT / "launch_sharded_recovery.ps1").read_text(encoding="ascii")
    assert "--instance-termination-action=STOP" in script
    assert "--max-run-duration=24h" in script
    assert "--machine-type=c3d-highcpu-16" in script
    assert "--provisioning-model=STANDARD" in script
    assert "DELETE" not in script
    assert "[string]$LinuxUser" in script
    assert 'LinuxUser = "' + "shreshth" + '_rach1"' not in script
    assert "Assert-NativeSuccess" in script
    assert "--zone=$worker.Zone" not in script
    zones = [
        "us-central1-a",
        "us-east1-b",
    ]
    assert all(zone in script for zone in zones)
    assert len({zone.rsplit("-", 1)[0] for zone in zones}) == 2


def test_collection_verifies_every_hash_before_merge() -> None:
    script = (ROOT / "collect_sharded_recovery.ps1").read_text(encoding="ascii")
    assert script.index("Get-FileHash -Algorithm SHA256") < script.index(
        "merge_sharded_recovery_archives.py"
    )
    assert 'if ($archives.Count -ne 2)' in script
    assert "compute instances stop" in script
    assert "compute instances delete" not in script
