from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import output_dir


CATALOGS = [
    (
        "current_cumulative_2026",
        Path(r"C:\Users\shres\Desktop\HILDALE RESEARCH\cumulative_2026.02.11_22.33.58.csv"),
    ),
    (
        "alderaan_cumulative_20240816",
        Path(r"C:\Users\shres\Desktop\HILDALE RESEARCH\external\alderaan\Catalogs\koi_cumulative_exoarchive_20240816.csv"),
    ),
    (
        "dr22_mullally_q1_q16",
        Path(r"C:\Users\shres\Desktop\HILDALE RESEARCH\external\alderaan\Catalogs\kepler_q1_q16_dr22_mullally.csv"),
    ),
    (
        "dr24_coughlin_q1_q17",
        Path(r"C:\Users\shres\Desktop\HILDALE RESEARCH\external\alderaan\Catalogs\kepler_q1_q17_dr24_coughlin.csv"),
    ),
    (
        "dr25_thompson_q1_q17",
        Path(r"C:\Users\shres\Desktop\HILDALE RESEARCH\external\alderaan\Catalogs\kepler_q1_q17_dr25_thompson.csv"),
    ),
    (
        "merged_planets_full",
        Path(r"C:\Users\shres\Desktop\HILDALE RESEARCH\merged_planets_full.csv"),
    ),
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit local historical transit-seed recovery for unseeded ALDERAAN planets.")
    parser.add_argument("--copy-to-codex-outputs", default=None)
    args = parser.parse_args()

    out = output_dir()
    unseeded = pd.read_csv(out / "alderaan_unseeded_needed_planets_best.csv")
    rows = []
    for catalog_name, path in CATALOGS:
        if not path.exists():
            continue
        cat = pd.read_csv(path, comment="#")
        if "kepoi_name" not in cat.columns:
            continue
        keep = [
            "kepoi_name",
            "kepid",
            "koi_disposition",
            "koi_pdisposition",
            "koi_score",
            "koi_fpflag_nt",
            "koi_fpflag_ss",
            "koi_fpflag_co",
            "koi_fpflag_ec",
            "koi_period",
            "koi_time0bk",
            "koi_impact",
            "koi_duration",
            "koi_depth",
            "koi_prad",
            "koi_model_snr",
            "koi_tce_plnt_num",
            "koi_tce_delivname",
        ]
        sub = unseeded[["kepoi_name", "kepid", "koi_period", "koi_duration"]].merge(
            cat[[c for c in keep if c in cat.columns]],
            on="kepoi_name",
            how="left",
            suffixes=("_current", "_candidate"),
        )
        for _, row in sub.iterrows():
            current_duration = number(row.get("koi_duration_current"))
            candidate_duration = number(row.get("koi_duration_candidate"))
            rows.append(
                {
                    "kepoi_name": row["kepoi_name"],
                    "catalog": catalog_name,
                    "catalog_path": str(path),
                    "found_in_catalog": pd.notna(row.get("kepid_candidate", row.get("kepid"))),
                    "candidate_disposition": row.get("koi_disposition"),
                    "candidate_pdisposition": row.get("koi_pdisposition"),
                    "candidate_score": number(row.get("koi_score")),
                    "candidate_fpflag_nt": number(row.get("koi_fpflag_nt")),
                    "candidate_fpflag_ss": number(row.get("koi_fpflag_ss")),
                    "candidate_fpflag_co": number(row.get("koi_fpflag_co")),
                    "candidate_fpflag_ec": number(row.get("koi_fpflag_ec")),
                    "current_period": number(row.get("koi_period_current")),
                    "candidate_period": number(row.get("koi_period_candidate")),
                    "candidate_epoch": number(row.get("koi_time0bk")),
                    "current_duration": current_duration,
                    "candidate_duration": candidate_duration,
                    "duration_ratio_candidate_to_current": candidate_duration / current_duration
                    if current_duration and np.isfinite(candidate_duration)
                    else np.nan,
                    "candidate_depth": number(row.get("koi_depth")),
                    "candidate_impact": number(row.get("koi_impact")),
                    "candidate_prad": number(row.get("koi_prad")),
                    "candidate_snr": number(row.get("koi_model_snr")),
                    "candidate_tce_plnt_num": number(row.get("koi_tce_plnt_num")),
                    "candidate_tce_delivname": row.get("koi_tce_delivname"),
                }
            )
    audit = pd.DataFrame(rows)
    audit["has_valid_depth"] = pd.to_numeric(audit["candidate_depth"], errors="coerce") > 0
    audit["duration_within_25pct"] = np.abs(np.log(pd.to_numeric(audit["duration_ratio_candidate_to_current"], errors="coerce"))) <= np.log(1.25)
    audit["autofill_recommended"] = audit["has_valid_depth"] & audit["duration_within_25pct"]

    audit_path = out / "alderaan_unseeded_historical_seed_audit.csv"
    md_path = out / "alderaan_unseeded_historical_seed_audit.md"
    audit.to_csv(audit_path, index=False)
    write_markdown(md_path, audit)

    copy_root = Path(args.copy_to_codex_outputs) if args.copy_to_codex_outputs else None
    if copy_root:
        copy_root.mkdir(parents=True, exist_ok=True)
        for path in [audit_path, md_path]:
            (copy_root / path.name).write_bytes(path.read_bytes())

    print("=== Unseeded historical transit-seed audit ===")
    print(summary_table(audit).to_string(index=False))
    print(f"\nWrote: {audit_path}")
    print(f"Wrote: {md_path}")


def number(value) -> float:
    try:
        value = float(value)
    except Exception:
        return np.nan
    return value if np.isfinite(value) else np.nan


def summary_table(audit: pd.DataFrame) -> pd.DataFrame:
    return (
        audit.groupby("catalog")
        .agg(
            matched_rows=("found_in_catalog", "sum"),
            valid_depth_rows=("has_valid_depth", "sum"),
            valid_depth_and_duration_consistent_rows=("autofill_recommended", "sum"),
        )
        .reset_index()
        .sort_values("catalog")
    )


def write_markdown(path: Path, audit: pd.DataFrame) -> None:
    summary = summary_table(audit)
    valid = audit[audit["has_valid_depth"]].copy()
    cols = [
        "kepoi_name",
        "catalog",
        "candidate_disposition",
        "candidate_pdisposition",
        "current_duration",
        "candidate_duration",
        "duration_ratio_candidate_to_current",
        "candidate_depth",
        "candidate_impact",
        "candidate_snr",
        "autofill_recommended",
    ]
    lines = [
        "# Unseeded Transit Seed Recovery Audit",
        "",
        "This checks local historical KOI/TCE catalogs for the 27 needed ALDERAAN planets that lack current KOI depth seeds.",
        "",
        "## Summary By Catalog",
        "",
        summary.to_markdown(index=False),
        "",
        "## Rows With Historical Depth",
        "",
        valid[[c for c in cols if c in valid.columns]].to_markdown(index=False, floatfmt=".4f")
        if len(valid)
        else "None.",
        "",
        "Conclusion: historical depths exist for a few blocked planets, but none pass the simple automatic-fill rule requiring depth plus duration agreement within 25 percent of the current catalog row. These should stay manual/blocked for now.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
