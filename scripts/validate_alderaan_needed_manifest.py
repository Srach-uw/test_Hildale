from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate the final pre-ALDERAAN needed-target manifest.")
    parser.add_argument("--copy-to-external-output", default=None)
    args = parser.parse_args()

    out = output_dir()
    sample = pd.read_csv(out / "canonical_sample_old_astropy_rawcc.csv")
    all_status = pd.read_csv(out / "alderaan_all_planets_status_best.csv")
    needed = pd.read_csv(out / "alderaan_needed_planets_best.csv")
    runnable = pd.read_csv(out / "alderaan_runnable_needed_planets_best.csv")
    unseeded = pd.read_csv(out / "alderaan_unseeded_needed_planets_best.csv")
    targets = pd.read_csv(out / "alderaan_needed_targets_best.csv")
    expanded = pd.read_csv(out / "alderaan_needed_catalog_rows_best.csv")
    catalog = pd.read_csv(out / "alderaan_needed_catalog_best.csv", index_col=0)

    needed_with_catalog = flag_needed_catalog_presence(runnable, catalog)
    missing_from_catalog = needed_with_catalog[~needed_with_catalog["present_in_alderaan_catalog"]].copy()

    checks = build_checks(sample, all_status, needed, runnable, unseeded, targets, expanded, catalog, missing_from_catalog)
    checks_path = out / "alderaan_needed_validation_checks.csv"
    missing_path = out / "alderaan_needed_planets_missing_from_catalog.csv"
    md_path = out / "alderaan_needed_validation.md"

    checks.to_csv(checks_path, index=False)
    missing_from_catalog.to_csv(missing_path, index=False)
    write_markdown(md_path, checks, missing_from_catalog, needed, targets, catalog)

    copy_root = Path(args.copy_to_external_output) if args.copy_to_external_output else None
    if copy_root:
        copy_root.mkdir(parents=True, exist_ok=True)
        for path in [checks_path, missing_path, md_path]:
            (copy_root / path.name).write_bytes(path.read_bytes())

    print("=== ALDERAAN needed manifest validation ===")
    print(checks.to_string(index=False))
    print(f"\nWrote: {checks_path}")
    print(f"Wrote: {missing_path}")
    print(f"Wrote: {md_path}")


def flag_needed_catalog_presence(needed: pd.DataFrame, catalog: pd.DataFrame) -> pd.DataFrame:
    left = needed.copy()
    left["_period_key"] = pd.to_numeric(left["koi_period"], errors="coerce").round(8)
    right = catalog.copy()
    right["_period_key"] = pd.to_numeric(right["period"], errors="coerce").round(8)
    keys = right[["koi_id", "kic_id", "_period_key"]].drop_duplicates()
    merged = left.merge(
        keys.assign(present_in_alderaan_catalog=True),
        left_on=["koi_target", "kepid", "_period_key"],
        right_on=["koi_id", "kic_id", "_period_key"],
        how="left",
    )
    merged["present_in_alderaan_catalog"] = merged["present_in_alderaan_catalog"].fillna(False)
    return merged.drop(columns=[c for c in ["koi_id", "kic_id", "_period_key"] if c in merged.columns])


def build_checks(
    sample: pd.DataFrame,
    all_status: pd.DataFrame,
    needed: pd.DataFrame,
    runnable: pd.DataFrame,
    unseeded: pd.DataFrame,
    targets: pd.DataFrame,
    expanded: pd.DataFrame,
    catalog: pd.DataFrame,
    missing_from_catalog: pd.DataFrame,
) -> pd.DataFrame:
    target_systems = set(targets["koi_target"])
    runnable_systems = set(runnable["koi_target"])
    expanded_systems = set(expanded["koi_target"])
    catalog_systems = set(catalog["koi_id"])

    required_catalog_cols = ["koi_id", "kic_id", "npl", "period", "epoch", "depth", "duration", "impact"]
    required_bad = 0
    for col in required_catalog_cols:
        required_bad += int(catalog[col].isna().sum())
    required_bad += int((pd.to_numeric(catalog["period"], errors="coerce") <= 0).sum())
    required_bad += int((pd.to_numeric(catalog["depth"], errors="coerce") <= 0).sum())
    required_bad += int((pd.to_numeric(catalog["duration"], errors="coerce") <= 0).sum())

    npl_check = catalog.groupby("koi_id").size().rename("actual_npl").reset_index()
    npl_declared = catalog.groupby("koi_id")["npl"].first().rename("declared_npl").reset_index()
    npl_bad = npl_check.merge(npl_declared, on="koi_id")
    npl_bad = npl_bad[npl_bad["actual_npl"] != npl_bad["declared_npl"]]

    rows = [
        check("sample_rows", len(sample), True, "Canonical best sample rows."),
        check("sample_unique_kepoi", sample["kepoi_name"].nunique(), sample["kepoi_name"].nunique() == len(sample), "No duplicate planet IDs expected."),
        check("all_status_rows", len(all_status), len(all_status) == len(sample), "Status table should mirror the sample."),
        check("all_status_unique_kepoi", all_status["kepoi_name"].nunique(), all_status["kepoi_name"].nunique() == len(all_status), "No duplicate status planet IDs."),
        check("needed_planets_rows", len(needed), int(all_status["needs_alderaan"].sum()) == len(needed), "Needed rows must equal status needs_alderaan count."),
        check("runnable_needed_planets_rows", len(runnable), int((needed["can_seed_alderaan"]).sum()) == len(runnable), "Runnable needed rows must have valid launch seeds."),
        check("unseeded_needed_planets_rows", len(unseeded), int((~needed["can_seed_alderaan"]).sum()) == len(unseeded), "Unseeded needed rows should be tracked separately."),
        check("needed_unique_kepoi", needed["kepoi_name"].nunique(), needed["kepoi_name"].nunique() == len(needed), "No duplicate needed planet IDs."),
        check("needed_targets_rows", len(targets), targets["koi_target"].nunique() == len(targets), "One row per KOI target."),
        check("runnable_planets_missing_target_row", len(runnable_systems - target_systems), len(runnable_systems - target_systems) == 0, "Every runnable needed planet target must be in target manifest."),
        check("target_rows_without_runnable_needed_planet", len(target_systems - runnable_systems), len(target_systems - runnable_systems) == 0, "Every target row should be justified by at least one runnable needed planet."),
        check("expanded_catalog_rows", len(expanded), expanded_systems == target_systems, "Expanded rows should include all sample planets in needed systems."),
        check("alderaan_catalog_rows", len(catalog), catalog_systems == target_systems, "ALDERAAN catalog should include every target system."),
        check("runnable_planets_missing_from_catalog", len(missing_from_catalog), len(missing_from_catalog) == 0, "Every runnable flagged planet should have valid ALDERAAN catalog parameters."),
        check("catalog_required_field_failures", required_bad, required_bad == 0, "Required catalog fields must be finite and positive where appropriate."),
        check("catalog_npl_mismatches", len(npl_bad), len(npl_bad) == 0, "Catalog npl should match rows per target."),
        check("system_expansion_extra_sibling_planets", len(catalog) - len(needed), len(catalog) >= len(needed), "Non-needed planets added because ALDERAAN fits whole systems."),
    ]
    return pd.DataFrame(rows)


def check(name: str, value: int, passed: bool, note: str) -> dict:
    return {"check": name, "value": int(value), "passed": bool(passed), "note": note}


def write_markdown(
    path: Path,
    checks: pd.DataFrame,
    missing_from_catalog: pd.DataFrame,
    needed: pd.DataFrame,
    targets: pd.DataFrame,
    catalog: pd.DataFrame,
) -> None:
    failures = checks[~checks["passed"]]
    tier_counts = (
        needed.groupby(["priority_tier", "recommended_action"])
        .size()
        .rename("planets")
        .reset_index()
        .sort_values(["priority_tier", "recommended_action"])
    )
    target_counts = (
        targets.groupby(["min_priority_tier", "target_run_type"])
        .size()
        .rename("targets")
        .reset_index()
        .sort_values(["min_priority_tier", "target_run_type"])
    )
    pop_counts = (
        needed.groupby(["disk", "system"])
        .size()
        .rename("needed_planets")
        .reset_index()
        .sort_values(["disk", "system"])
    )
    lines = [
        "# ALDERAAN Needed Manifest Validation",
        "",
        "## Pass/Fail Checks",
        "",
        checks.to_markdown(index=False),
        "",
        "## Summary",
        "",
        f"- Needed planet rows: {len(needed)}",
        f"- Runnable needed planet rows: {int(needed['can_seed_alderaan'].sum())}",
        f"- Unseeded needed planet rows: {int((~needed['can_seed_alderaan']).sum())}",
        f"- Unique KOI systems to run/refit: {len(targets)}",
        f"- ALDERAAN catalog rows after system expansion: {len(catalog)}",
        f"- Extra sibling planets included by system expansion: {len(catalog) - len(needed)}",
        f"- Failed validation checks: {len(failures)}",
        "",
        "## Needed Planets By Tier",
        "",
        tier_counts.to_markdown(index=False),
        "",
        "## Needed Targets By Tier",
        "",
        target_counts.to_markdown(index=False),
        "",
        "## Needed Planets By Population",
        "",
        pop_counts.to_markdown(index=False),
        "",
    ]
    if len(missing_from_catalog):
        lines += [
            "## Missing From ALDERAAN Catalog",
            "",
            missing_from_catalog[["kepoi_name", "koi_target", "kepid", "koi_period", "recommended_action"]]
            .head(50)
            .to_markdown(index=False),
            "",
        ]
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
