from __future__ import annotations

import argparse
import hashlib
import shutil
import tarfile
import tempfile
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_extract(archive: Path, destination: Path) -> None:
    destination = destination.resolve()
    with tarfile.open(archive, "r:*") as handle:
        for member in handle.getmembers():
            target = (destination / member.name).resolve()
            if destination != target and destination not in target.parents:
                raise ValueError(f"archive member escapes destination: {member.name}")
        handle.extractall(destination)


def result_target(path: Path) -> str:
    suffix = "-results.fits"
    if not path.name.endswith(suffix):
        raise ValueError(path)
    return path.name[: -len(suffix)]


def expected_target_set(expected_targets: pd.DataFrame) -> set[str]:
    if "koi_target" not in expected_targets:
        raise ValueError("target table is missing required koi_target column")
    targets = expected_targets["koi_target"].astype(str).str.strip()
    if targets.eq("").any():
        raise ValueError("target table contains blank koi_target values")
    if targets.duplicated().any():
        duplicates = sorted(targets[targets.duplicated(keep=False)].unique())
        raise ValueError(f"target table contains duplicate koi_target values: {duplicates}")
    return set(targets)


def validate_result_layout(result: Path, run_id: str, target: str) -> None:
    expected_tail = ("Results", run_id, target, result.name)
    if result.parts[-4:] != expected_tail:
        raise ValueError(
            f"result is not under Results/{run_id}/{target}: {result}"
        )


def merge_archives(
    archives: list[Path],
    expected_targets: pd.DataFrame,
    staging: Path,
    run_id: str,
) -> pd.DataFrame:
    expected = expected_target_set(expected_targets)
    seen: dict[str, dict[str, object]] = {}
    staging.mkdir(parents=True, exist_ok=True)
    results_out = staging / "Results" / run_id
    status_out = staging / "status"
    logs_out = staging / "logs"
    provenance_out = staging / "provenance"
    for path in (results_out, status_out, logs_out, provenance_out):
        path.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="alderaan_shard_merge_") as tmp:
        tmp_root = Path(tmp)
        for archive_index, archive in enumerate(archives):
            extracted = tmp_root / f"archive_{archive_index:02d}"
            extracted.mkdir()
            safe_extract(archive, extracted)
            archive_hash = sha256(archive)

            for result in extracted.rglob("*-results.fits"):
                if result.stat().st_size == 0:
                    continue
                target = result_target(result)
                validate_result_layout(result, run_id, target)
                if target not in expected:
                    raise ValueError(f"unexpected result target {target} in {archive}")
                if target in seen:
                    raise ValueError(
                        f"duplicate target {target} occurs in both "
                        f"{seen[target]['archive']} and {archive}"
                    )
                target_dir = results_out / target
                target_dir.mkdir(parents=True, exist_ok=True)
                destination = target_dir / result.name
                shutil.copy2(result, destination)
                seen[target] = {
                    "koi_target": target,
                    "archive": str(archive),
                    "archive_sha256": archive_hash,
                    "result_sha256": sha256(destination),
                    "result_bytes": destination.stat().st_size,
                }

            for status in extracted.rglob("*.status"):
                target = status.stem
                if target in expected and not (status_out / status.name).exists():
                    shutil.copy2(status, status_out / status.name)

            for source_name, destination_root in (
                ("logs", logs_out),
                ("provenance", provenance_out),
            ):
                roots = [p for p in extracted.rglob(source_name) if p.is_dir()]
                for root in roots:
                    shard_root = destination_root / f"archive_{archive_index:02d}"
                    shutil.copytree(root, shard_root, dirs_exist_ok=True)

    rows = []
    for target in sorted(expected):
        if target in seen:
            rows.append({**seen[target], "merge_status": "result_present"})
        else:
            rows.append(
                {
                    "koi_target": target,
                    "archive": "",
                    "archive_sha256": "",
                    "result_sha256": "",
                    "result_bytes": 0,
                    "merge_status": "missing_result",
                }
            )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(staging / "shard_merge_manifest.csv", index=False)
    expected_targets.to_csv(staging / "targets_missing_launchable.csv", index=False)
    return manifest


def write_tar(staging: Path, output: Path) -> None:
    with tarfile.open(output, "w:gz") as handle:
        for path in sorted(staging.rglob("*")):
            if path.is_file():
                handle.add(path, arcname=path.relative_to(staging).as_posix())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("archives", nargs="+", type=Path)
    parser.add_argument(
        "--targets",
        default="cloud_published_inventory_missing_batch/targets_missing_launchable.csv",
        type=Path,
    )
    parser.add_argument(
        "--run-id",
        default="sagear_published_inventory_missing",
    )
    parser.add_argument(
        "--staging",
        default="tmp/merged_sharded_recovery",
        type=Path,
    )
    parser.add_argument(
        "--output",
        default="outputs/alderaan_results_published_inventory_missing_merged.tar.gz",
        type=Path,
    )
    args = parser.parse_args()

    archives = [path.resolve() for path in args.archives]
    if len(archives) != len(set(archives)):
        raise ValueError("the same archive was supplied more than once")
    for archive in archives:
        if not archive.is_file():
            raise FileNotFoundError(archive)

    staging = args.staging.resolve()
    if staging.exists() and any(staging.iterdir()):
        raise ValueError(f"staging directory is not empty: {staging}")
    expected = pd.read_csv(args.targets)
    manifest = merge_archives(
        archives,
        expected_targets=expected,
        staging=staging,
        run_id=args.run_id,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_tar(staging, args.output)
    present = int(manifest["merge_status"].eq("result_present").sum())
    missing = int(manifest["merge_status"].eq("missing_result").sum())
    print(f"MERGED: {present} results; {missing} missing; {len(expected)} expected")
    print(f"ARCHIVE: {args.output.resolve()}")
    print(f"SHA256: {sha256(args.output.resolve())}")


if __name__ == "__main__":
    main()
