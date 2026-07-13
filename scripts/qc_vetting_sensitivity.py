from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import output_dir
from hierarchical_rayleigh import (
    POPULATIONS,
    fit_from_mass_matrix,
    load_population_masses,
    recompute_grid_support_flags,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Test whether transparent fit-quality/vetting rules can recover Sagear-like Rayleigh values."
    )
    parser.add_argument(
        "--summary",
        default=None,
        help="Merged posterior summary. Defaults to paired-exact QC-primary output.",
    )
    parser.add_argument(
        "--sample",
        default=None,
        help="Canonical sample table used to add KOI reliability fields.",
    )
    parser.add_argument(
        "--leverage",
        default=None,
        help="Per-planet leverage table for top-fraction veto tests.",
    )
    parser.add_argument("--tag", default="QC_VETTING")
    parser.add_argument("--n-sigma", type=int, default=2000)
    args = parser.parse_args()

    out = output_dir()
    summary_path = Path(args.summary) if args.summary else out / "eccentricity_posterior_summary_merged_paired_exact_qcprimary.csv"
    sample_path = Path(args.sample) if args.sample else out / "canonical_sample_old_astropy_rawcc.csv"
    leverage_path = Path(args.leverage) if args.leverage else out / "rayleigh_per_planet_leverage_PAIRED_EXACT_qcprimary.csv"

    df = pd.read_csv(summary_path)
    df = recompute_grid_support_flags(df)
    df["population"] = df["disk"].astype(str) + "_" + df["system"].astype(str) + "s"
    if sample_path.exists():
        sample_cols = [
            "kepoi_name",
            "koi_score",
            "koi_disposition",
            "koi_pdisposition",
            "koi_prad",
            "koi_model_snr",
            "koi_num_transits",
            "koi_comment",
            "koi_fpflag_nt",
            "koi_fpflag_ss",
            "koi_fpflag_co",
            "koi_fpflag_ec",
        ]
        sample = pd.read_csv(sample_path, usecols=lambda c: c in sample_cols).drop_duplicates("kepoi_name")
        df = df.merge(sample, on="kepoi_name", how="left", suffixes=("", "_sample"))

    leverage = pd.read_csv(leverage_path) if leverage_path.exists() else pd.DataFrame()
    rules = build_rules(df, leverage)
    sigmas = np.linspace(1e-4, 1.0, args.n_sigma)

    count_rows = []
    fit_rows = []
    for rule_name, keep in rules.items():
        sub = df[keep].reset_index(drop=True)
        count_rows.extend(counts_for_rule(rule_name, df, sub))
        fit_rows.extend(fit_rule(rule_name, sub, sigmas))
        filtered_path = out / f"eccentricity_posterior_summary_filter_{rule_name}.csv"
        sub.drop(columns=["population"], errors="ignore").to_csv(filtered_path, index=False)

    count_table = pd.DataFrame(count_rows)
    fit_table = pd.DataFrame(fit_rows)
    count_path = out / f"qc_vetting_sensitivity_counts_{args.tag}.csv"
    fit_path = out / f"qc_vetting_sensitivity_rayleigh_{args.tag}.csv"
    pivot_path = out / f"qc_vetting_sensitivity_rayleigh_pivot_{args.tag}.csv"
    count_table.to_csv(count_path, index=False)
    fit_table.to_csv(fit_path, index=False)
    pivot = fit_table.pivot(index="rule", columns="population", values="expected_e").reset_index()
    pivot.to_csv(pivot_path, index=False)

    print(f"Wrote {count_path}")
    print(f"Wrote {fit_path}")
    print(f"Wrote {pivot_path}")
    print("\nCounts:")
    print(count_table.to_string(index=False))
    print("\nRayleigh expected <e>:")
    print(pivot.to_string(index=False))


def build_rules(df: pd.DataFrame, leverage: pd.DataFrame) -> dict[str, pd.Series]:
    base = pd.Series(True, index=df.index)
    rules: dict[str, pd.Series] = {"BASE_QC_PRIMARY": base}
    qcreasons = df["qc_reasons"].fillna("") if "qc_reasons" in df.columns else pd.Series("", index=df.index)

    rules["NO_ZETA_TAIL_REVIEW"] = base & ~qcreasons.str.contains("zeta_tail", regex=False)
    rules["NO_E84_GT_0P9_REVIEW"] = base & ~qcreasons.str.contains("e84_gt_0p9", regex=False)
    rules["NO_ZETA_OR_E84_REVIEW"] = (
        base
        & ~qcreasons.str.contains("zeta_tail", regex=False)
        & ~qcreasons.str.contains("e84_gt_0p9", regex=False)
    )

    impact_fit = pd.to_numeric(df.get("impact_fit_p50", np.nan), errors="coerce")
    impact_catalog = pd.to_numeric(df.get("catalog_impact", np.nan), errors="coerce")
    non_grazing = ((impact_fit.isna()) | (impact_fit < 0.9)) & ((impact_catalog.isna()) | (impact_catalog < 0.9))
    rules["NO_GRAZING_B_GE_0P9"] = base & non_grazing

    score = pd.to_numeric(df.get("koi_score", np.nan), errors="coerce")
    disposition = df.get("koi_disposition", pd.Series("", index=df.index)).fillna("").astype(str)
    reliable = (disposition == "CONFIRMED") | (score >= 0.5)
    rules["CONFIRMED_OR_SCORE_GE_0P5"] = base & reliable

    prad = pd.to_numeric(df.get("koi_prad", np.nan), errors="coerce")
    rules["PLANET_RADIUS_LT_8RE"] = base & ((prad.isna()) | (prad < 8.0))
    rules["PLANET_RADIUS_LT_4RE"] = base & ((prad.isna()) | (prad < 4.0))

    rules["TRANSPARENT_STRICT_FIT_QC"] = (
        base
        & non_grazing
        & ~qcreasons.str.contains("zeta_tail", regex=False)
        & ~qcreasons.str.contains("e84_gt_0p9", regex=False)
        & ~qcreasons.str.contains("period_mismatch", regex=False)
        & ~df.get("incomplete_system", pd.Series(False, index=df.index)).fillna(False).astype(bool)
    )
    rules["TRANSPARENT_STRICT_ASTRO_QC"] = (
        rules["TRANSPARENT_STRICT_FIT_QC"]
        & reliable
        & ((prad.isna()) | (prad < 8.0))
    )

    if not leverage.empty and "delta_loglike_map_minus_min" in leverage.columns:
        rules["DROP_TOP_2PCT_LEVERAGE_PER_POP"] = base & ~top_fraction_kois(df, leverage, 0.02)
        rules["DROP_TOP_5PCT_LEVERAGE_PER_POP"] = base & ~top_fraction_kois(df, leverage, 0.05)
        rules["DROP_TOP_10PCT_LEVERAGE_PER_POP"] = base & ~top_fraction_kois(df, leverage, 0.10)

    return rules


def top_fraction_kois(df: pd.DataFrame, leverage: pd.DataFrame, frac: float) -> pd.Series:
    lev = leverage.copy()
    if "population" not in lev:
        lev["population"] = lev["disk"].astype(str) + "_" + lev["system"].astype(str) + "s"
    drop_kois: set[str] = set()
    for population, g in lev.groupby("population"):
        n = max(1, int(np.ceil(frac * len(g))))
        top = g.sort_values("delta_loglike_map_minus_min", ascending=False).head(n)
        drop_kois.update(top["kepoi_name"].astype(str))
    return df["kepoi_name"].astype(str).isin(drop_kois)


def counts_for_rule(rule_name: str, original: pd.DataFrame, filtered: pd.DataFrame) -> list[dict]:
    rows = []
    for _, _, label in POPULATIONS:
        before = original[original["population"] == label]
        after = filtered[filtered["population"] == label]
        rows.append(
            {
                "rule": rule_name,
                "population": label,
                "n_before": int(len(before)),
                "n_after": int(len(after)),
                "n_dropped": int(len(before) - len(after)),
                "frac_dropped": float((len(before) - len(after)) / len(before)) if len(before) else np.nan,
            }
        )
    return rows


def fit_rule(rule_name: str, df: pd.DataFrame, sigmas: np.ndarray) -> list[dict]:
    rows = []
    for disk, system, label in POPULATIONS:
        sub = df[(df["disk"] == disk) & (df["system"] == system)].reset_index(drop=True)
        if len(sub) < 5:
            rows.append({"rule": rule_name, "population": label, "n": len(sub), "status": "skipped_n_lt_5"})
            continue
        masses, e_grid = load_population_masses(sub, apply_transit_selection=True)
        fit = fit_from_mass_matrix(masses, e_grid, sigmas, apply_transit_selection=True)
        rows.append({"rule": rule_name, "population": label, "n": len(sub), "status": "ok", **fit})
    return rows


if __name__ == "__main__":
    main()
