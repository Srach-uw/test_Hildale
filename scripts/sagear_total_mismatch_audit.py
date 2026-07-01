from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common import load_config, output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit the 2474 vs 2465 Sagear total-planet mismatch.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--summary", default=None)
    parser.add_argument("--shape", default=None)
    parser.add_argument("--main-tex", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = output_dir()
    sample_path = Path(args.sample) if args.sample else out / "canonical_sample_old_astropy_rawcc.csv"
    summary_path = Path(args.summary) if args.summary else out / "eccentricity_posterior_summary_old_astropy_rawcc.csv"
    shape_path = Path(args.shape) if args.shape else out / "alderaan_shape_diagnostics.csv"
    tex_path = Path(args.main_tex) if args.main_tex else Path(cfg["_root"]) / "main.tex"

    macros = parse_macros(tex_path)
    sample = pd.read_csv(sample_path)
    summary = pd.read_csv(summary_path) if summary_path.exists() else pd.DataFrame()
    shape = pd.read_csv(shape_path) if shape_path.exists() else pd.DataFrame()

    macro_audit = build_macro_audit(macros, sample)
    category_audit = build_category_audit(macros, sample)
    implied_removed = build_implied_removed(macros)
    candidates = build_candidate_planets(sample, summary, shape)
    host_candidates = build_host_candidates(candidates)

    macro_path = out / "sagear_total_mismatch_macro_audit.csv"
    cat_path = out / "sagear_total_mismatch_category_audit.csv"
    implied_path = out / "sagear_total_mismatch_implied_removed.csv"
    cand_path = out / "sagear_total_mismatch_candidate_planets.csv"
    host_path = out / "sagear_total_mismatch_candidate_hosts.csv"
    md_path = out / "sagear_total_mismatch_audit.md"

    macro_audit.to_csv(macro_path, index=False)
    category_audit.to_csv(cat_path, index=False)
    implied_removed.to_csv(implied_path, index=False)
    candidates.to_csv(cand_path, index=False)
    host_candidates.to_csv(host_path, index=False)
    write_markdown(md_path, macro_audit, category_audit, implied_removed, host_candidates, cand_path)

    print("=== Macro total audit ===")
    print(macro_audit.to_string(index=False))
    print("\n=== Category audit ===")
    print(category_audit.to_string(index=False))
    print("\n=== Implied manuscript removed/missing bucket ===")
    print(implied_removed.to_string(index=False))
    print("\n=== Top candidate hosts to inspect, not confirmed Sagear IDs ===")
    print(host_candidates.head(20).to_string(index=False))
    print(f"\nWrote: {macro_path}")
    print(f"Wrote: {cat_path}")
    print(f"Wrote: {implied_path}")
    print(f"Wrote: {cand_path}")
    print(f"Wrote: {host_path}")
    print(f"Wrote: {md_path}")


def parse_macros(tex_path: Path) -> dict[str, int]:
    text = tex_path.read_text(errors="replace")
    return {
        match.group(1): int(match.group(2))
        for match in re.finditer(r"\\newcommand\\(all\w+)\{(\d+)\s*\}", text)
    }


def build_macro_audit(macros: dict[str, int], sample: pd.DataFrame) -> pd.DataFrame:
    thin = int(macros["allthinplanets"])
    thick = int(macros["allthickplanets"])
    subgroup = (
        int(macros["allthinsingles"])
        + int(macros["allthicksingles"])
        + int(macros["allthinmultiplanets"])
        + int(macros["allthickmultiplanets"])
    )
    rows = [
        {
            "quantity": "our_current_best_sample_rows",
            "value": int(len(sample)),
            "delta_vs_our_sample": 0,
            "interpretation": "Current best reconstructed sample size.",
        },
        {
            "quantity": "sagear_allplanets_macro",
            "value": int(macros["allplanets"]),
            "delta_vs_our_sample": int(len(sample) - macros["allplanets"]),
            "interpretation": "The total stated in the sample-selection text.",
        },
        {
            "quantity": "sagear_disk_total_macros",
            "value": thin + thick,
            "delta_vs_our_sample": int(len(sample) - (thin + thick)),
            "interpretation": "allthinplanets + allthickplanets. This equals our current total.",
        },
        {
            "quantity": "sagear_subgroup_total_macros",
            "value": subgroup,
            "delta_vs_our_sample": int(len(sample) - subgroup),
            "interpretation": "thin/thick singles + multis. This equals allplanets, not disk total.",
        },
    ]
    return pd.DataFrame(rows)


def build_category_audit(macros: dict[str, int], sample: pd.DataFrame) -> pd.DataFrame:
    targets = {
        "thin_singles": int(macros["allthinsingles"]),
        "thick_singles": int(macros["allthicksingles"]),
        "thin_multis": int(macros["allthinmultiplanets"]),
        "thick_multis": int(macros["allthickmultiplanets"]),
        "thin_total_planets": int(macros["allthinplanets"]),
        "thick_total_planets": int(macros["allthickplanets"]),
        "thin_multi_hosts": int(macros["allthinmultistars"]),
        "thick_multi_hosts": int(macros["allthickmultistars"]),
    }
    ours = {
        "thin_singles": int(((sample["disk"] == "thin") & (sample["system"] == "single")).sum()),
        "thick_singles": int(((sample["disk"] == "thick") & (sample["system"] == "single")).sum()),
        "thin_multis": int(((sample["disk"] == "thin") & (sample["system"] == "multi")).sum()),
        "thick_multis": int(((sample["disk"] == "thick") & (sample["system"] == "multi")).sum()),
        "thin_total_planets": int((sample["disk"] == "thin").sum()),
        "thick_total_planets": int((sample["disk"] == "thick").sum()),
        "thin_multi_hosts": int(sample.loc[(sample["disk"] == "thin") & (sample["system"] == "multi"), "kepid"].nunique()),
        "thick_multi_hosts": int(sample.loc[(sample["disk"] == "thick") & (sample["system"] == "multi"), "kepid"].nunique()),
    }
    rows = []
    for key, target in targets.items():
        rows.append({"metric": key, "ours": ours[key], "sagear_macro": target, "delta_ours_minus_sagear": ours[key] - target})
    rows.append(
        {
            "metric": "removal_only_can_match_subgroups",
            "ours": 0,
            "sagear_macro": 0,
            "delta_ours_minus_sagear": 0,
            "note": (
                "No. We are below Sagear in thin singles, thick singles, and thin multis, but above in thick multis. "
                "Deleting nine planets alone cannot recover the subgroup targets; labels/classification also differ."
            ),
        }
    )
    return pd.DataFrame(rows)


def build_implied_removed(macros: dict[str, int]) -> pd.DataFrame:
    rows = [
        {
            "bucket": "thick_planets_in_disk_total_but_not_subgroups",
            "count": int(macros["allthickplanets"] - macros["allthicksingles"] - macros["allthickmultiplanets"]),
            "calculation": "allthickplanets - allthicksingles - allthickmultiplanets",
            "interpretation": "Nine thick planets are present in the thick total macro but absent from final subgroup planet macros.",
        },
        {
            "bucket": "thick_hosts_in_disk_total_but_not_subgroups",
            "count": int(macros["allthickstars"] - macros["allthicksingles"] - macros["allthickmultistars"]),
            "calculation": "allthickstars - allthicksingles - allthickmultistars",
            "interpretation": "Five thick hosts are present in the thick host total macro but absent from subgroup host macros.",
        },
    ]
    return pd.DataFrame(rows)


def build_candidate_planets(sample: pd.DataFrame, summary: pd.DataFrame, shape: pd.DataFrame) -> pd.DataFrame:
    df = sample[sample["disk"] == "thick"].copy()
    if not summary.empty:
        has_post = set(summary["kepoi_name"].astype(str))
        df["has_existing_posterior"] = df["kepoi_name"].astype(str).isin(has_post)
        ecols = ["kepoi_name", "e50", "zeta_median"]
        df = df.merge(summary[[c for c in ecols if c in summary.columns]], on="kepoi_name", how="left")
    else:
        df["has_existing_posterior"] = False
    if not shape.empty:
        keep = [
            "kepoi_name",
            "alderaan_to_koi_duration_ratio",
            "alderaan_to_koi_ror_ratio",
            "alderaan_impact_med",
            "koi_model_snr",
        ]
        df = df.merge(shape[[c for c in keep if c in shape.columns]], on="kepoi_name", how="left")

    if "koi_model_snr_x" in df.columns and "koi_model_snr" not in df.columns:
        df["koi_model_snr"] = df["koi_model_snr_x"]

    df["missing_posterior_score"] = (~df["has_existing_posterior"]).astype(int) * 3
    df["low_snr_score"] = (numeric_series(df, "koi_model_snr") < 20).astype(int)
    df["duration_shift_score"] = (
        np.abs(np.log(numeric_series(df, "alderaan_to_koi_duration_ratio"))) > np.log(1.25)
    ).astype(int)
    df["ror_shift_score"] = (
        np.abs(np.log(numeric_series(df, "alderaan_to_koi_ror_ratio"))) > np.log(1.25)
    ).astype(int)
    df["impact_score"] = (numeric_series(df, "alderaan_impact_med") > 0.85).astype(int)
    df["candidate_score"] = (
        df["missing_posterior_score"]
        + df["low_snr_score"]
        + df["duration_shift_score"]
        + df["ror_shift_score"]
        + df["impact_score"]
    )

    keep_cols = [
        "kepid",
        "kepoi_name",
        "disk",
        "system",
        "koi_period",
        "koi_model_snr",
        "koi_prad",
        "P_thick",
        "has_existing_posterior",
        "e50",
        "zeta_median",
        "alderaan_to_koi_duration_ratio",
        "alderaan_to_koi_ror_ratio",
        "alderaan_impact_med",
        "candidate_score",
    ]
    return df[[c for c in keep_cols if c in df.columns]].sort_values(
        ["candidate_score", "system", "P_thick"], ascending=[False, False, False]
    )


def numeric_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.nan, index=df.index)
    return pd.to_numeric(df[col], errors="coerce")


def build_host_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        candidates.groupby("kepid")
        .agg(
            candidate_planets=("kepoi_name", "count"),
            systems=("system", lambda x: ",".join(sorted(set(map(str, x))))),
            planet_names=("kepoi_name", lambda x: ",".join(map(str, x))),
            missing_posterior_planets=("has_existing_posterior", lambda x: int((~x).sum())),
            max_candidate_score=("candidate_score", "max"),
            sum_candidate_score=("candidate_score", "sum"),
            median_p_thick=("P_thick", "median"),
            median_snr=("koi_model_snr", "median"),
        )
        .reset_index()
    )
    grouped["host_set_relevance"] = (
        grouped["missing_posterior_planets"] * 3 + grouped["sum_candidate_score"] + grouped["candidate_planets"]
    )
    return grouped.sort_values(["host_set_relevance", "candidate_planets"], ascending=[False, False])


def write_markdown(
    path: Path,
    macro_audit: pd.DataFrame,
    category_audit: pd.DataFrame,
    implied_removed: pd.DataFrame,
    host_candidates: pd.DataFrame,
    candidate_path: Path,
) -> None:
    lines = [
        "# Sagear Total Mismatch Audit",
        "",
        "## Macro Total Audit",
        "",
        macro_audit.to_markdown(index=False),
        "",
        "## Category Audit",
        "",
        category_audit.to_markdown(index=False),
        "",
        "## Implied Manuscript Bucket",
        "",
        implied_removed.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "- Our current best sample has 2474 planets.",
        "- Sagear's `allthinplanets + allthickplanets` also equals 2474.",
        "- Sagear's `allplanets` and subgroup planet macros equal 2465.",
        "- The nine-planet mismatch is internal to the manuscript macros and is specifically in the thick-disk aggregate.",
        "- The manuscript source archive does not include an actual final planet list, so the exact nine IDs cannot be recovered from local source files alone.",
        "- Deleting nine planets from our current sample cannot make the subgroup counts match, because some categories are already below Sagear.",
        "",
        "## Candidate Hosts To Inspect",
        "",
        "These are not confirmed Sagear removals. They are current thick-disk hosts ranked by missing posterior / low SNR / shape-pathology proxies.",
        "",
        host_candidates.head(20).to_markdown(index=False, floatfmt=".4f"),
        "",
        f"Full candidate planet table: `{candidate_path}`",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
