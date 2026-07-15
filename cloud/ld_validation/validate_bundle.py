#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", default="targets_missing_launchable.csv")
    parser.add_argument("--catalog", default="sagear_missing_catalog.csv")
    args = parser.parse_args()

    targets = pd.read_csv(args.targets)
    catalog = pd.read_csv(args.catalog, index_col=0)
    failures = []
    if targets["koi_target"].duplicated().any():
        failures.append("duplicate koi_target rows in targets")
    missing_catalog_targets = set(targets["koi_target"]) - set(catalog["koi_id"])
    if missing_catalog_targets:
        failures.append(f"{len(missing_catalog_targets)} target(s) missing from catalog")
    for col in ["koi_id", "kic_id", "npl", "period", "epoch", "depth", "duration", "impact"]:
        if col not in catalog.columns:
            failures.append(f"catalog missing column {col}")
        elif catalog[col].isna().any():
            failures.append(f"catalog column {col} has NaNs")
    for col in ["period", "depth", "duration"]:
        if col in catalog.columns and (pd.to_numeric(catalog[col], errors="coerce") <= 0).any():
            failures.append(f"catalog column {col} has non-positive values")
    npl = catalog.groupby("koi_id").size().rename("actual").reset_index()
    decl = catalog.groupby("koi_id")["npl"].first().rename("declared").reset_index()
    bad = npl.merge(decl, on="koi_id")
    bad = bad[bad["actual"] != bad["declared"]]
    if len(bad):
        failures.append(f"{len(bad)} target(s) have npl mismatch")
    for col in ["limbdark_1", "limbdark_2"]:
        if col in catalog.columns:
            nunique = catalog.groupby("koi_id")[col].nunique(dropna=False)
            inconsistent = nunique[nunique > 1]
            if len(inconsistent):
                failures.append(
                    f"{len(inconsistent)} target(s) have inconsistent system-level {col}"
                )
    if failures:
        print("VALIDATION FAILED")
        for failure in failures:
            print("-", failure)
        sys.exit(2)
    print(f"VALIDATION OK: {len(targets)} targets, {len(catalog)} catalog rows")


if __name__ == "__main__":
    main()
