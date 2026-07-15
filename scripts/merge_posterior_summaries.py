from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge existing and newly extracted ALDERAAN e,omega posterior summaries.")
    parser.add_argument(
        "--base",
        default=None,
        help="Existing posterior summary CSV. Defaults to eccentricity_posterior_summary_old_astropy_rawcc.csv.",
    )
    parser.add_argument("--new", required=True, help="New posterior summary CSV, e.g. cloud extraction output.")
    parser.add_argument("--out", default=None, help="Merged output CSV path.")
    parser.add_argument("--coverage-out", default=None, help="Merged coverage output CSV path.")
    parser.add_argument("--sample", default=None, help="Canonical sample path for coverage.")
    parser.add_argument(
        "--replacements-out",
        default=None,
        help="Output CSV path listing kepoi_names present in both base and new (new supersedes old).",
    )
    parser.add_argument("--new-excluded", default=None, help="Optional extraction-exclusion CSV from the new run.")
    parser.add_argument("--manifest-out", default=None, help="Per-sample planet posterior/QC manifest CSV.")
    parser.add_argument("--qc-usable-out", default=None, help="Merged summary filtered to canonical primary-QC usable rows.")
    parser.add_argument("--veto-kois", default=None, help="Comma-separated kepoi_name values to mark as vetted QC vetoes.")
    args = parser.parse_args()

    out_dir = output_dir()
    base_path = Path(args.base) if args.base else out_dir / "eccentricity_posterior_summary_old_astropy_rawcc.csv"
    new_path = Path(args.new)
    out_path = Path(args.out) if args.out else out_dir / "eccentricity_posterior_summary_merged.csv"
    coverage_path = Path(args.coverage_out) if args.coverage_out else out_dir / "eccentricity_posterior_coverage_merged.csv"
    sample_path = Path(args.sample) if args.sample else out_dir / "canonical_sample_old_astropy_rawcc.csv"
    replacements_path = (
        Path(args.replacements_out) if args.replacements_out else out_dir / "eccentricity_posterior_merge_replacements.csv"
    )
    manifest_path = Path(args.manifest_out) if args.manifest_out else out_dir / "eccentricity_posterior_qc_manifest.csv"
    qc_usable_path = Path(args.qc_usable_out) if args.qc_usable_out else out_dir / "eccentricity_posterior_summary_qc_primary.csv"

    base = pd.read_csv(base_path)
    new = pd.read_csv(new_path)
    sample = pd.read_csv(sample_path)

    # kepoi_names present in both: the new ALDERAAN run refit an existing
    # archive planet (e.g. a joint fit for a missing-target system that
    # includes an already-modeled sibling). The dedup below keeps the new
    # row, i.e. treats it as a superseding refit -- make that explicit rather
    # than letting 30-ish silent row replacements hide inside a rowcount.
    replaced_ids = set(base["kepoi_name"]) & set(new["kepoi_name"])
    if replaced_ids:
        replacements = base[base["kepoi_name"].isin(replaced_ids)].merge(
            new[new["kepoi_name"].isin(replaced_ids)],
            on="kepoi_name",
            suffixes=("_old", "_new"),
        )
        replacements_path.parent.mkdir(parents=True, exist_ok=True)
        replacements.to_csv(replacements_path, index=False)

    base["posterior_source"] = "existing_archive"
    new["posterior_source"] = "new_alderaan"
    merged = pd.concat([base, new], ignore_index=True, sort=False)
    merged["_source_rank"] = merged["posterior_source"].map({"existing_archive": 0, "new_alderaan": 1}).fillna(0)
    merged = (
        merged.sort_values(["kepoi_name", "_source_rank"])
        .drop_duplicates("kepoi_name", keep="last")
        .drop(columns=["_source_rank"])
        .sort_values(["disk", "system", "koi_target", "koi_period"])
        .reset_index(drop=True)
    )

    # Boolean QC columns (e.g. incomplete_system) only exist on rows from
    # whichever extraction run introduced them; concat leaves the other
    # source's rows as NaN, which `.astype(bool)` downstream would silently
    # read as True. Default any such column to False for rows that lack it.
    bool_like_cols = [
        "incomplete_system",
        "zeta_median_outside_grid_support",
        "zeta_p16_outside_grid_support",
        "zeta_p84_outside_grid_support",
        "zeta_any_summary_outside_grid_support",
    ]
    for col in bool_like_cols:
        if col in merged.columns:
            merged[col] = merged[col].astype("boolean").fillna(False).astype(bool)

    veto_kois = {k.strip() for k in str(args.veto_kois or "").split(",") if k.strip()}
    merged = annotate_qc(merged, veto_kois=veto_kois)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    coverage_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)
    qc_usable = merged[~merged["qc_primary_exclude"]].reset_index(drop=True)
    qc_usable.to_csv(qc_usable_path, index=False)
    coverage = coverage_summary(sample, merged)
    coverage.to_csv(coverage_path, index=False)
    manifest = build_manifest(sample, merged, Path(args.new_excluded) if args.new_excluded else None)
    manifest.to_csv(manifest_path, index=False)

    print(f"Base rows: {len(base)}")
    print(f"New rows: {len(new)}")
    print(f"Merged unique posterior rows: {len(merged)}")
    print(f"Replacement rows (kepoi_name in both base and new, new kept): {len(replaced_ids)}")
    if replaced_ids:
        print(f"Wrote: {replacements_path}")
    print(f"Wrote: {out_path}")
    print(f"Wrote: {qc_usable_path}")
    print(f"Wrote: {coverage_path}")
    print(f"Wrote: {manifest_path}")
    print("\nCoverage:")
    print(coverage.to_string(index=False))


GRID_SUPPORT_FALLBACK_MIN = float(np.sqrt(1.0 - 0.95**2) / (1.0 + 0.95))
GRID_SUPPORT_FALLBACK_MAX = float(np.sqrt(1.0 - 0.95**2) / (1.0 - 0.95))


def annotate_qc(summary: pd.DataFrame, veto_kois: set[str] | None = None) -> pd.DataFrame:
    out = summary.copy()
    veto_kois = veto_kois or set()
    for col in ["zeta_median", "zeta_p16", "zeta_p84", "e84", "period_relative_difference"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    gmin = out["zeta_grid_min"] if "zeta_grid_min" in out.columns else pd.Series(np.nan, index=out.index)
    gmax = out["zeta_grid_max"] if "zeta_grid_max" in out.columns else pd.Series(np.nan, index=out.index)
    out["zeta_grid_min"] = pd.to_numeric(gmin, errors="coerce").fillna(GRID_SUPPORT_FALLBACK_MIN)
    out["zeta_grid_max"] = pd.to_numeric(gmax, errors="coerce").fillna(GRID_SUPPORT_FALLBACK_MAX)
    out["zeta_median_outside_grid_support"] = ~out["zeta_median"].between(out["zeta_grid_min"], out["zeta_grid_max"])
    out["zeta_p16_outside_grid_support"] = ~out["zeta_p16"].between(out["zeta_grid_min"], out["zeta_grid_max"])
    out["zeta_p84_outside_grid_support"] = ~out["zeta_p84"].between(out["zeta_grid_min"], out["zeta_grid_max"])
    out["zeta_any_summary_outside_grid_support"] = (
        out["zeta_median_outside_grid_support"]
        | out["zeta_p16_outside_grid_support"]
        | out["zeta_p84_outside_grid_support"]
    )
    out["qc_flag_zeta_tail"] = (out["zeta_median"] < 0.7) | (out["zeta_median"] > 1.3)
    out["qc_flag_e84_gt_0p9"] = out["e84"] > 0.9 if "e84" in out.columns else False
    out["qc_flag_period_mismatch"] = (
        out["period_relative_difference"] > 1e-3 if "period_relative_difference" in out.columns else False
    )
    out["qc_flag_vetted_veto"] = out["kepoi_name"].astype(str).isin(veto_kois)
    if "incomplete_system" not in out.columns:
        out["incomplete_system"] = pd.NA
    incomplete = out["incomplete_system"].astype("boolean")
    out["incomplete_system_unknown"] = incomplete.isna()
    # Unknown completeness is not evidence of a complete ALDERAAN system.
    out["incomplete_system"] = incomplete.fillna(True).astype(bool)
    out["qc_primary_exclude"] = (
        out["zeta_any_summary_outside_grid_support"]
        | out["incomplete_system"]
        | out["qc_flag_vetted_veto"]
    )
    reasons = []
    for _, row in out.iterrows():
        r = []
        if bool(row.get("zeta_any_summary_outside_grid_support", False)):
            r.append("zeta_summary_outside_e095_support")
        if bool(row.get("incomplete_system_unknown", False)):
            r.append("incomplete_system_unknown")
        elif bool(row.get("incomplete_system", False)):
            r.append("incomplete_system")
        if bool(row.get("qc_flag_vetted_veto", False)):
            r.append("vetted_qc_veto")
        if bool(row.get("qc_flag_zeta_tail", False)):
            r.append("zeta_tail_review")
        if bool(row.get("qc_flag_e84_gt_0p9", False)):
            r.append("e84_gt_0p9_review")
        if bool(row.get("qc_flag_period_mismatch", False)):
            r.append("period_mismatch_review")
        reasons.append(";".join(r))
    out["qc_reasons"] = reasons
    return out


def build_manifest(sample: pd.DataFrame, summary: pd.DataFrame, excluded_path: Path | None) -> pd.DataFrame:
    keep_cols = [
        "kepoi_name",
        "posterior_source",
        "posterior_file",
        "e16",
        "e50",
        "e84",
        "zeta_median",
        "zeta_p16",
        "zeta_p84",
        "impact_mode",
        "qc_primary_exclude",
        "qc_reasons",
    ]
    available = summary[[c for c in keep_cols if c in summary.columns]].copy()
    manifest = sample.merge(available, on="kepoi_name", how="left")
    manifest["posterior_status"] = np.where(manifest["posterior_file"].notna(), "posterior_available", "missing_after_merge")
    manifest["qc_primary_exclude"] = manifest["qc_primary_exclude"].astype("boolean").fillna(True).astype(bool)
    manifest["qc_reasons"] = manifest["qc_reasons"].fillna("missing_after_merge")
    if excluded_path is not None and excluded_path.exists():
        excluded = pd.read_csv(excluded_path)
        if {"kepoi_name", "reason"}.issubset(excluded.columns):
            reason_map = excluded.drop_duplicates("kepoi_name").set_index("kepoi_name")["reason"].to_dict()
            missing = manifest["posterior_status"].eq("missing_after_merge")
            manifest.loc[missing, "qc_reasons"] = (
                manifest.loc[missing, "kepoi_name"].map(reason_map).fillna(manifest.loc[missing, "qc_reasons"])
            )
    return manifest


def coverage_summary(sample: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    summary_ids = set(summary["kepoi_name"])
    for (disk, system), sub in sample.groupby(["disk", "system"]):
        posterior_planets = int(sub["kepoi_name"].isin(summary_ids).sum())
        rows.append(
            {
                "disk": disk,
                "system": system,
                "sample_planets": int(len(sub)),
                "posterior_planets": posterior_planets,
                "missing_planets": int(len(sub) - posterior_planets),
                "coverage_fraction": posterior_planets / len(sub) if len(sub) else 0.0,
            }
        )
    return pd.DataFrame(rows).sort_values(["disk", "system"])


if __name__ == "__main__":
    main()
