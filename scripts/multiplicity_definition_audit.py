from __future__ import annotations

"""Audit the single/multi definition used by the published Sagear inventory.

The KOI catalog's ``koi_count`` is a host-level catalog field.  It is not
equivalent to recounting only the rows that survive the planet-level cuts.
This audit keeps those definitions separate and records which one reproduces
the published pre-visual category counts.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import load_config, output_dir, read_koi, read_sagear2026_kinematic_hosts


PUBLISHED = {
    ("thin", "single"): 1121,
    ("thick", "single"): 275,
    ("thin", "multi"): 883,
    ("thick", "multi"): 212,
}


def published_pre_visual_rows(koi: pd.DataFrame, hosts: pd.DataFrame) -> pd.DataFrame:
    host_ids = set(pd.to_numeric(hosts["kepid"], errors="coerce").dropna().astype(int))
    out = koi[
        koi["kepid"].isin(host_ids)
        & koi["koi_disposition"].astype(str).str.upper().ne("FALSE POSITIVE")
        & pd.to_numeric(koi["koi_period"], errors="coerce").between(1.0, 100.0)
    ].copy()
    out = out.merge(hosts[["kepid", "disk_published"]], on="kepid", validate="many_to_one")
    out["disk"] = out["disk_published"]
    return out


def definition_counts(rows: pd.DataFrame, host_counts: pd.Series) -> dict[tuple[str, str], int]:
    system = rows["kepid"].map(host_counts).eq(1).map({True: "single", False: "multi"})
    return rows.assign(system=system).groupby(["disk", "system"]).size().to_dict()


def build_audit(koi: pd.DataFrame, hosts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = published_pre_visual_rows(koi, hosts)
    host_ids = set(hosts["kepid"].astype(int))
    scoped = koi[koi["kepid"].isin(host_ids)].copy()
    scoped["koi_count"] = pd.to_numeric(scoped["koi_count"], errors="coerce")
    definitions = {
        "catalog_koi_count": scoped.groupby("kepid")["koi_count"].first(),
        "all_catalog_rows": scoped.groupby("kepid")["kepoi_name"].size(),
        "nonfp_all_periods": scoped[
            scoped["koi_disposition"].astype(str).str.upper().ne("FALSE POSITIVE")
        ].groupby("kepid")["kepoi_name"].size(),
        "nonfp_period_1_to_100": rows.groupby("kepid")["kepoi_name"].size(),
    }
    count_rows: list[dict[str, object]] = []
    for name, counts in definitions.items():
        counts = counts.reindex(rows["kepid"].unique())
        observed = definition_counts(rows, counts)
        count_rows.append(
            {
                "definition": name,
                "rows": int(len(rows)),
                "hosts": int(rows["kepid"].nunique()),
                "matches_published_pre_visual_counts": observed == PUBLISHED,
                **{f"{disk}_{system}": int(observed.get((disk, system), 0)) for disk, system in PUBLISHED},
            }
        )

    comparison = pd.DataFrame(
        {
            "kepid": rows["kepid"].drop_duplicates().astype(int),
        }
    ).set_index("kepid")
    for name, counts in definitions.items():
        comparison[name] = counts.reindex(comparison.index)
    comparison["koi_count_vs_nonfp_all_periods_diff"] = (
        comparison["catalog_koi_count"] != comparison["nonfp_all_periods"]
    )
    comparison["koi_count_vs_nonfp_period_1_to_100_diff"] = (
        comparison["catalog_koi_count"] != comparison["nonfp_period_1_to_100"]
    )
    return pd.DataFrame(count_rows), comparison.reset_index()


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit Sagear published-inventory multiplicity definitions.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--counts-out", default=None)
    parser.add_argument("--hosts-out", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    counts, hosts = build_audit(read_koi(cfg), read_sagear2026_kinematic_hosts(cfg))
    if not bool(counts.loc[counts["definition"] == "catalog_koi_count", "matches_published_pre_visual_counts"].iloc[0]):
        raise AssertionError("catalog koi_count no longer reproduces Sagear's pre-visual category counts")
    out = output_dir()
    counts_path = Path(args.counts_out) if args.counts_out else out / "multiplicity_definition_counts.csv"
    hosts_path = Path(args.hosts_out) if args.hosts_out else out / "multiplicity_definition_host_comparison.csv"
    counts.to_csv(counts_path, index=False)
    hosts.to_csv(hosts_path, index=False)
    print(counts.to_string(index=False))
    print(f"\nHosts with catalog/non-FP multiplicity disagreement: {int(hosts['koi_count_vs_nonfp_all_periods_diff'].sum())}")
    print(f"Wrote: {counts_path}")
    print(f"Wrote: {hosts_path}")


if __name__ == "__main__":
    main()
