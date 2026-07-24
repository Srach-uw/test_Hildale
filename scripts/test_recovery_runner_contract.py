from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER_DIRS = [
    ROOT / "cloud" / "ld_validation",
    ROOT / "cloud" / "recovery_142",
]


def test_two_pass_runner_preserves_canonical_then_retry_order() -> None:
    for runner_dir in RUNNER_DIRS:
        script = (runner_dir / "run_recovery_two_pass.sh").read_text(encoding="ascii")
        canonical = script.index("SEED_OFFSET=0 bash run_batch.sh")
        retry = script.index("SEED_OFFSET=1 bash run_batch.sh")
        package = script.index("bash pack_results.sh")
        assert canonical < retry < package


def test_pack_results_survives_incomplete_summary() -> None:
    for runner_dir in RUNNER_DIRS:
        script = (runner_dir / "pack_results.sh").read_text(encoding="ascii")
        assert "set +e\nbash summarize_progress.sh" in script
        assert "SUMMARY_RC=$?" in script
        assert script.index('tar -czf "$OUT"') < script.index(
            'if [ "$SUMMARY_RC" -ne 0 ]'
        )


def test_deployment_guide_uses_two_pass_v4_bundle() -> None:
    for runner_dir in RUNNER_DIRS:
        guide = (runner_dir / "DEPLOY_ON_EXISTING_VM.md").read_text(encoding="ascii")
        assert "Missing_142_20260724_v4.zip" in guide
        assert "nohup bash run_recovery_two_pass.sh" in guide
        assert "run_recovery_two_pass|run_batch" in guide
        assert "16-hour" not in guide


def test_shard_runner_preserves_one_target_partition() -> None:
    for runner_dir in RUNNER_DIRS:
        script = (runner_dir / "run_shard.sh").read_text(encoding="ascii")
        assert 'TARGET_CSV="target_shards/targets_shard_${SHARD_PADDED}.csv"' in script
        assert 'export RUN_ID="sagear_published_inventory_missing"' in script
        assert "bash run_recovery_two_pass.sh" in script


def test_packaged_audit_includes_shard_assignments() -> None:
    for runner_dir in RUNNER_DIRS:
        script = (runner_dir / "pack_results.sh").read_text(encoding="ascii")
        assert "target_shards" in script


def test_shard_bootstrap_preserves_disk_and_predictable_archive() -> None:
    for runner_dir in RUNNER_DIRS:
        script = (runner_dir / "bootstrap_shard.sh").read_text(encoding="ascii")
        assert 'FINAL_ARCHIVE="$HOME/alderaan_shard_${SHARD_PADDED}_results.tar.gz"' in script
        assert 'sha256sum "$FINAL_ARCHIVE"' in script
        assert "sudo shutdown -h now" in script
        assert "delete" not in script.lower()
