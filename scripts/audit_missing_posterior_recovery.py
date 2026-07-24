"""Audit the exact published-inventory ALDERAAN recovery queue.

This deliberately audits the fixed 2,491-planet inventory rather than the
older broad 592-system queue.  A system is the unit of an ALDERAAN run, while
the population result is counted in planets, so both levels are reported.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", default="outputs/eccentricity_posterior_manifest_equal_nested_published_inventory_pre_visual_qc.csv")
    parser.add_argument("--batch", default="cloud_published_inventory_missing_batch")
    parser.add_argument("--output", default="outputs/missing_posterior_recovery_audit.csv")
    args = parser.parse_args()

    root = args.root.resolve()
    manifest = pd.read_csv(root / args.manifest)
    batch = root / args.batch
    queue = pd.read_csv(batch / "targets_missing_launchable.csv")
    missing = manifest.loc[manifest["posterior_status"] == "missing_after_uniform_assembly"].copy()

    expected_missing = set(missing["koi_target"])
    queued = set(queue["koi_target"])
    if expected_missing != queued:
        raise SystemExit(
            "queue mismatch: "
            f"missing_not_queued={len(expected_missing - queued)}, "
            f"queued_not_missing={len(queued - expected_missing)}"
        )

    group = (
        missing.groupby("koi_target", dropna=False)
        .agg(
            missing_population_planets=("koi_target", "size"),
            paper_category_target=("paper_category_target", "first"),
        )
        .reset_index()
    )
    system = queue.merge(group, on=["koi_target"], how="left", validate="one_to_one")
    system["queue_target"] = True
    system["missing_planets"] = system["missing_planets"].astype(int)
    system["prior_retry"] = system["prior_launch_status"].eq(
        "previously_launched_no_usable_result"
    )
    system["new_system"] = system["prior_launch_status"].eq(
        "not_in_previous_592_target_launch"
    )

    if system["missing_population_planets"].isna().any():
        raise SystemExit("some queued systems did not join to the missing population rows")

    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    system.sort_values(["disk", "system", "koi_target"]).to_csv(output, index=False)

    summary = {
        "missing_population_planets": len(missing),
        "missing_systems": missing["koi_target"].nunique(),
        "queued_systems": len(queue),
        "queued_catalog_rows": int(pd.read_csv(batch / "sagear_missing_catalog.csv").shape[0]),
        "new_systems": int(system["new_system"].sum()),
        "retry_systems": int(system["prior_retry"].sum()),
    }
    print("MISSING POSTERIOR RECOVERY AUDIT")
    for key, value in summary.items():
        print(f"{key},{value}")
    print("\nBY POPULATION GROUP")
    print(
        missing.groupby(["disk", "system", "paper_category_target"], dropna=False)
        .size()
        .rename("missing_planets")
        .to_string()
    )
    print(f"\nmanifest,{output}")


if __name__ == "__main__":
    main()
