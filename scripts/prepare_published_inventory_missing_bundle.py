from __future__ import annotations

import argparse
import hashlib
import shutil
from pathlib import Path

import pandas as pd

from alderaan_batch import build_alderaan_catalog
from build_recovery_shards import add_runtime_proxy, assign_lpt_shards, validate_partition
from build_limb_darkening_validation_batch import expand_to_full_koi_systems
from common import load_config
from prepare_gcp_missing_alderaan_bundle import build_targets


REQUIRED_RUNNER_FILES = [
    "bootstrap_shard.sh",
    "DEPLOY_ON_EXISTING_VM.md",
    "run_batch.sh",
    "run_one_target.sh",
    "setup_vm.sh",
    "summarize_progress.sh",
    "validate_bundle.py",
    "pack_results.sh",
    "run_recovery_two_pass.sh",
    "run_shard.sh",
    "patch_alderaan_repro.py",
    "patch_alderaan_paper_priors.py",
]


def classify_retry_targets(
    targets: pd.DataFrame, previously_launched: set[str]
) -> pd.DataFrame:
    out = targets.copy()
    out["prior_launch_status"] = out["koi_target"].map(
        lambda target: (
            "previously_launched_no_usable_result"
            if str(target) in previously_launched
            else "not_in_previous_592_target_launch"
        )
    )
    return out


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_readme(
    path: Path,
    targets: pd.DataFrame,
    catalog: pd.DataFrame,
    inventory: pd.DataFrame,
) -> None:
    counts = (
        targets.groupby(["disk", "system", "prior_launch_status"])
        .size()
        .rename("targets")
        .reset_index()
    )
    lines = [
        "# Exact published-inventory ALDERAAN completion batch",
        "",
        "This bundle covers every published-host system that still lacks a usable",
        "ALDERAAN result in the local archive plus the completed cloud run.",
        "",
        f"- target systems: {len(targets)}",
        f"- population planets missing posteriors: {int(targets['missing_planets'].sum())}",
        f"- full-system ALDERAAN catalog rows: {len(catalog)}",
        f"- full-system inventory rows including unseedable KOIs: {len(inventory)}",
        "- cadence mode: long plus available short cadence",
        "- ALDERAAN commit: 7443dff16b7f9092e14a6f0cc1f8948d457c9e0b",
        "- runner: repaired resumable target-specific logging and deterministic seed",
        "",
        "Multiplicity and population labels come from the exact published host",
        "inventory. The transit fit includes every non-false-positive, seedable",
        "KOI companion in each selected system before the 1-100 day population cut.",
        "",
        "## Target breakdown",
        "",
        "| disk | system | prior launch status | targets |",
        "|---|---|---|---:|",
    ]
    for _, row in counts.iterrows():
        lines.append(
            f"| {row['disk']} | {row['system']} | "
            f"{row['prior_launch_status']} | {int(row['targets'])} |"
        )
    lines.extend(
        [
            "",
            "## Run",
            "",
            "After setup and environment activation:",
            "",
            "```bash",
            "JOBS=20 nohup bash run_recovery_two_pass.sh > published_inventory_missing.log 2>&1 &",
            "```",
            "",
            "For a multi-VM run, assign exactly one shard to each VM:",
            "",
            "```bash",
            "nohup bash run_shard.sh 00 28 > shard_00.log 2>&1 &",
            "```",
            "",
            "The target shards are deterministic, mutually exclusive, and their",
            "union is validated against the complete 142-target master list.",
            "",
            "This runs the canonical seed, retries only unresolved targets with",
            "`SEED_OFFSET=1`, and packages all results and failure evidence.",
            "The run is resumable. Existing nonempty result FITS are skipped.",
            "Failures are recorded with their active stage (`download`, `detrend`,",
            "`noise`, or `transit_fit`) and are never silently replaced by a",
            "different noise or transit model.",
            "",
            "A successful retry records its seed offset in the target provenance.",
            "Do not alter GP jitter, priors, or the transit model in the canonical",
            "branch. Any such numerical rescue belongs in a separate sensitivity",
            "run.",
            "",
            "For a manual resume or extra snapshot, package the current results",
            "at any time with:",
            "",
            "```bash",
            "RUN_ID=sagear_published_inventory_missing bash pack_results.sh",
            "```",
            "",
            "On Windows, process the downloaded archive from the canonical",
            "`sagear_reproduction` directory with:",
            "",
            "```powershell",
            'python .\\postprocess_published_inventory_recovery.py --archive "C:\\path\\to\\alderaan_results.tar.gz"',
            "```",
            "",
            "The postprocessor withholds the final hierarchy until the exact",
            "published-inventory coverage gate is complete.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare the exact published-host missing-posterior batch."
    )
    parser.add_argument("--config", default=None)
    parser.add_argument(
        "--sample",
        default="outputs/sagear2026_planet_inventory_pre_visual_qc.csv",
    )
    parser.add_argument(
        "--missing",
        default="outputs/sagear2026_planets_missing_alderaan_target_results.csv",
    )
    parser.add_argument(
        "--previous-targets",
        default="cloud_missing_batch/targets_missing_launchable.csv",
    )
    parser.add_argument(
        "--runner-template",
        default="cloud_ld_validation_batch",
    )
    parser.add_argument(
        "--out",
        default="cloud_published_inventory_missing_batch",
    )
    parser.add_argument("--shards", type=int, default=5)
    args = parser.parse_args()

    cfg = load_config(args.config)
    sample = pd.read_csv(args.sample)
    missing = pd.read_csv(args.missing)
    previous = pd.read_csv(args.previous_targets)
    targets = build_targets(missing, sample)
    targets = classify_retry_targets(
        targets, set(previous["koi_target"].astype(str))
    )

    selected_rows = sample.loc[
        sample["koi_target"].isin(targets["koi_target"])
    ].copy()
    base_catalog = build_alderaan_catalog(selected_rows, cfg)
    full_catalog, full_inventory = expand_to_full_koi_systems(
        base_catalog, targets, cfg
    )

    target_set = set(targets["koi_target"])
    if set(full_catalog["koi_id"]) != target_set:
        missing_catalog = sorted(target_set - set(full_catalog["koi_id"]))
        raise ValueError(
            f"full-system catalog is missing targets: {missing_catalog[:10]}"
        )
    declared = full_catalog.groupby("koi_id")["npl"].first()
    actual = full_catalog.groupby("koi_id").size()
    if not declared.equals(actual):
        raise ValueError("ALDERAAN npl declarations do not match catalog rows")

    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    template = Path(args.runner_template).resolve()
    for name in REQUIRED_RUNNER_FILES:
        source = template / name
        if not source.exists():
            raise FileNotFoundError(source)
        shutil.copy2(source, out / name)

    targets.to_csv(out / "targets_missing_launchable.csv", index=False)
    full_catalog.to_csv(out / "sagear_missing_catalog.csv")
    selected_rows.to_csv(
        out / "published_inventory_missing_population_rows.csv", index=False
    )
    full_inventory.to_csv(
        out / "published_inventory_missing_full_system_inventory.csv",
        index=False,
    )
    sharded = assign_lpt_shards(
        add_runtime_proxy(targets, full_catalog), args.shards
    )
    validate_partition(targets, sharded)
    shard_dir = out / "target_shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    for old in shard_dir.glob("targets_shard_*.csv"):
        old.unlink()
    for shard_id, rows in sharded.groupby("shard_id", sort=True):
        rows.to_csv(
            shard_dir / f"targets_shard_{int(shard_id):02d}.csv",
            index=False,
        )
    (
        sharded.groupby("shard_id", as_index=False)
        .agg(
            targets=("koi_target", "size"),
            planets_missing=("missing_planets", "sum"),
            catalog_planets=("catalog_planets", "sum"),
            runtime_proxy=("runtime_proxy", "sum"),
        )
        .to_csv(shard_dir / "shard_summary.csv", index=False)
    )
    sharded.to_csv(shard_dir / "shard_assignment.csv", index=False)
    write_readme(out / "README.md", targets, full_catalog, full_inventory)

    manifest_rows = []
    manifest_name = "bundle_manifest_sha256.csv"
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != manifest_name:
            manifest_rows.append(
                {
                    "file": path.relative_to(out).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    pd.DataFrame(manifest_rows).to_csv(
        out / manifest_name, index=False
    )

    print(
        targets.groupby(["prior_launch_status", "disk", "system"])
        .size()
        .rename("targets")
        .to_string()
    )
    print(f"Targets: {len(targets)}")
    print(f"Missing population planets: {int(targets['missing_planets'].sum())}")
    print(f"Full-system catalog rows: {len(full_catalog)}")
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
