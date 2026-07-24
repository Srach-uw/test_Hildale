from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def add_runtime_proxy(targets: pd.DataFrame, catalog: pd.DataFrame) -> pd.DataFrame:
    """Add a queue-balancing proxy without changing target eligibility."""
    grouped = (
        catalog.groupby("koi_id", as_index=False)
        .agg(
            catalog_planets=("npl", "first"),
            mean_depth_ppm=("depth", "mean"),
        )
        .rename(columns={"koi_id": "koi_target"})
    )
    out = targets.merge(grouped, on="koi_target", how="left", validate="one_to_one")
    if out[["catalog_planets", "mean_depth_ppm"]].isna().any().any():
        missing = out.loc[out["catalog_planets"].isna(), "koi_target"].tolist()
        raise ValueError(f"targets missing from ALDERAAN catalog: {missing[:10]}")

    # In the completed 24-target validation, multiplicity was the dominant
    # observable runtime predictor (Pearson r about 0.84). Shallow systems were
    # also slower on average. This score is used only for LPT queue balancing.
    out["runtime_proxy"] = (
        1.0
        + 1.5 * (out["catalog_planets"].astype(float) - 1.0)
        + 0.35 * (out["mean_depth_ppm"].astype(float) < 500.0)
        + 0.20
        * out["prior_launch_status"]
        .eq("previously_launched_no_usable_result")
        .astype(float)
    )
    return out


def assign_lpt_shards(targets: pd.DataFrame, shard_count: int) -> pd.DataFrame:
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    if targets["koi_target"].duplicated().any():
        raise ValueError("koi_target must be unique before sharding")

    loads = [0.0] * shard_count
    counts = [0] * shard_count
    assignments: dict[str, int] = {}
    ordered = targets.sort_values(
        ["runtime_proxy", "catalog_planets", "koi_target"],
        ascending=[False, False, True],
        kind="stable",
    )
    for row in ordered.itertuples(index=False):
        shard = min(range(shard_count), key=lambda i: (loads[i], counts[i], i))
        assignments[str(row.koi_target)] = shard
        loads[shard] += float(row.runtime_proxy)
        counts[shard] += 1

    out = targets.copy()
    out["shard_id"] = out["koi_target"].map(assignments).astype(int)
    return out.sort_values(["shard_id", "runtime_proxy"], ascending=[True, False])


def validate_partition(master: pd.DataFrame, partitioned: pd.DataFrame) -> None:
    if len(master) != len(partitioned):
        raise ValueError("partition row count differs from master target count")
    if partitioned["koi_target"].duplicated().any():
        raise ValueError("a target occurs in more than one shard")
    if set(master["koi_target"]) != set(partitioned["koi_target"]):
        raise ValueError("shard union does not equal the master target set")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--targets",
        default="cloud_published_inventory_missing_batch/targets_missing_launchable.csv",
    )
    parser.add_argument(
        "--catalog",
        default="cloud_published_inventory_missing_batch/sagear_missing_catalog.csv",
    )
    parser.add_argument("--shards", type=int, default=5)
    parser.add_argument(
        "--out-dir",
        default="cloud_published_inventory_missing_batch/target_shards",
    )
    args = parser.parse_args()

    targets = pd.read_csv(args.targets)
    catalog = pd.read_csv(args.catalog)
    enriched = add_runtime_proxy(targets, catalog)
    partitioned = assign_lpt_shards(enriched, args.shards)
    validate_partition(targets, partitioned)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("targets_shard_*.csv"):
        old.unlink()
    for shard_id, rows in partitioned.groupby("shard_id", sort=True):
        rows.to_csv(out / f"targets_shard_{int(shard_id):02d}.csv", index=False)

    summary = (
        partitioned.groupby("shard_id", as_index=False)
        .agg(
            targets=("koi_target", "size"),
            planets_missing=("missing_planets", "sum"),
            catalog_planets=("catalog_planets", "sum"),
            runtime_proxy=("runtime_proxy", "sum"),
            prior_failures=(
                "prior_launch_status",
                lambda s: int(
                    s.eq("previously_launched_no_usable_result").sum()
                ),
            ),
        )
    )
    summary.to_csv(out / "shard_summary.csv", index=False)
    partitioned.to_csv(out / "shard_assignment.csv", index=False)
    print(summary.to_string(index=False))
    print(f"VALID: {len(partitioned)} targets assigned exactly once")


if __name__ == "__main__":
    main()
