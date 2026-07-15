"""Audit ALDERAAN limb-darkening inputs against the bundled Kepler-Gaia catalog.

Sagear describes quadratic limb-darkening priors as system-level Normal priors
centered on values derived from Gaia/Berger stellar parameters and stellar
atmosphere models. The local ALDERAAN checkout includes a Kepler DR25/Gaia DR2
catalog with such-looking limb-darkening centers. This diagnostic compares the
catalogs used for our ALDERAAN runs against that bundled reference and checks
whether high-leverage eccentricity planets are concentrated in large LD-offset
systems.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import load_config, root_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", default="outputs/eccentricity_posterior_summary_merged_paired_exact_qcprimary.csv")
    parser.add_argument("--leverage", default="outputs/rayleigh_per_planet_leverage_PAIRED_EXACT_qcprimary.csv")
    parser.add_argument("--out-tag", default="QC_PRIMARY")
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parent


def resolve(path: str | Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else repo_root() / p


def catalog_paths(cfg: dict) -> list[tuple[str, Path]]:
    root = root_path(cfg, "research_root")
    if root is None:
        root = repo_root().parent
    return [
        ("missing_fixed", repo_root() / "cloud_missing_batch" / "sagear_missing_catalog_FIXED.csv"),
        ("missing_original", repo_root() / "cloud_missing_batch" / "sagear_missing_catalog.csv"),
        ("cloud_initial", repo_root() / "cloud_batch" / "sagear_cloud_catalog.csv"),
        ("needed_local", repo_root() / "alderaan_project" / "Catalogs" / "sagear_needed_catalog.csv"),
    ]


def load_reference_catalog(cfg: dict) -> pd.DataFrame:
    alderaan_repo = root_path(cfg, "alderaan_repo")
    if alderaan_repo is None:
        alderaan_repo = repo_root().parent / "external" / "alderaan"
    path = alderaan_repo / "Catalogs" / "kepler_dr25_gaia_dr2_crossmatch.csv"
    if not path.exists():
        raise FileNotFoundError(f"ALDERAAN bundled catalog not found: {path}")
    ref = pd.read_csv(path)
    ref = ref.copy()
    ref["kic_id"] = pd.to_numeric(ref["kic_id"], errors="coerce").astype("Int64")
    ref["period"] = pd.to_numeric(ref["period"], errors="coerce")
    ref["ref_limbdark_1"] = pd.to_numeric(ref["limbdark_1"], errors="coerce")
    ref["ref_limbdark_2"] = pd.to_numeric(ref["limbdark_2"], errors="coerce")
    return ref[["kic_id", "period", "ref_limbdark_1", "ref_limbdark_2", "planet_name", "koi_id"]]


def match_catalog_to_reference(name: str, path: Path, ref: pd.DataFrame) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    cat = pd.read_csv(path)
    cat = cat.copy()
    cat["catalog_name"] = name
    cat["catalog_path"] = str(path)
    cat["kic_id"] = pd.to_numeric(cat["kic_id"], errors="coerce").astype("Int64")
    cat["period"] = pd.to_numeric(cat["period"], errors="coerce")
    cat["limbdark_1"] = pd.to_numeric(cat["limbdark_1"], errors="coerce")
    cat["limbdark_2"] = pd.to_numeric(cat["limbdark_2"], errors="coerce")

    rows = []
    ref_groups = {int(k): g for k, g in ref.dropna(subset=["kic_id"]).groupby("kic_id")}
    for _, row in cat.iterrows():
        if pd.isna(row["kic_id"]) or not np.isfinite(row["period"]):
            continue
        candidates = ref_groups.get(int(row["kic_id"]))
        if candidates is None or candidates.empty:
            continue
        idx = (candidates["period"] - row["period"]).abs().idxmin()
        matched = candidates.loc[idx]
        period_delta = abs(float(matched["period"]) - float(row["period"]))
        if period_delta > 1e-3:
            continue
        rows.append(
            {
                "catalog_name": name,
                "koi_target": row.get("koi_id"),
                "kic_id": int(row["kic_id"]),
                "period": float(row["period"]),
                "limbdark_1": float(row["limbdark_1"]),
                "limbdark_2": float(row["limbdark_2"]),
                "ref_limbdark_1": float(matched["ref_limbdark_1"]),
                "ref_limbdark_2": float(matched["ref_limbdark_2"]),
                "ld_du1": float(row["limbdark_1"] - matched["ref_limbdark_1"]),
                "ld_du2": float(row["limbdark_2"] - matched["ref_limbdark_2"]),
                "ld_abs_max_delta": float(
                    max(
                        abs(row["limbdark_1"] - matched["ref_limbdark_1"]),
                        abs(row["limbdark_2"] - matched["ref_limbdark_2"]),
                    )
                ),
                "period_match_delta_days": period_delta,
                "ref_planet_name": matched["planet_name"],
            }
        )
    return pd.DataFrame(rows)


def summarize_offsets(matched: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, grp in matched.groupby("catalog_name"):
        rows.append(
            {
                "catalog_name": name,
                "matched_rows": len(grp),
                "median_du1": grp["ld_du1"].median(),
                "median_du2": grp["ld_du2"].median(),
                "q16_du1": grp["ld_du1"].quantile(0.16),
                "q84_du1": grp["ld_du1"].quantile(0.84),
                "q16_du2": grp["ld_du2"].quantile(0.16),
                "q84_du2": grp["ld_du2"].quantile(0.84),
                "frac_abs_delta_gt_0p05": float((grp["ld_abs_max_delta"] > 0.05).mean()),
                "frac_abs_delta_gt_0p10": float((grp["ld_abs_max_delta"] > 0.10).mean()),
                "max_abs_delta": grp["ld_abs_max_delta"].max(),
            }
        )
    return pd.DataFrame(rows).sort_values("catalog_name")


def attach_population_context(matched: pd.DataFrame, summary_path: Path, leverage_path: Path) -> pd.DataFrame:
    if matched.empty or not summary_path.exists():
        return pd.DataFrame()
    summary = pd.read_csv(summary_path)
    summary = summary.copy()
    summary["period"] = pd.to_numeric(summary["koi_period"], errors="coerce")
    summary["kic_id"] = pd.to_numeric(summary["kepid"], errors="coerce").astype("Int64")

    if leverage_path.exists():
        leverage = pd.read_csv(leverage_path)
        keep = ["kepoi_name", "delta_loglike_map_minus_min"]
        summary = summary.merge(leverage[keep], on="kepoi_name", how="left")

    rows = []
    matched_by_kic = {int(k): g for k, g in matched.groupby("kic_id")}
    for _, row in summary.iterrows():
        if pd.isna(row["kic_id"]) or not np.isfinite(row["period"]):
            continue
        candidates = matched_by_kic.get(int(row["kic_id"]))
        if candidates is None or candidates.empty:
            continue
        idx = (candidates["period"] - row["period"]).abs().idxmin()
        hit = candidates.loc[idx]
        if abs(float(hit["period"]) - float(row["period"])) > 1e-3:
            continue
        out = row.to_dict()
        for col in [
            "catalog_name",
            "limbdark_1",
            "limbdark_2",
            "ref_limbdark_1",
            "ref_limbdark_2",
            "ld_du1",
            "ld_du2",
            "ld_abs_max_delta",
        ]:
            out[col] = hit[col]
        out["population"] = f"{row['disk']}_{row['system']}s"
        rows.append(out)
    return pd.DataFrame(rows)


def summarize_by_population(context: pd.DataFrame) -> pd.DataFrame:
    if context.empty:
        return pd.DataFrame()
    rows = []
    for (catalog_name, population, source), grp in context.groupby(["catalog_name", "population", "posterior_source"]):
        rows.append(
            {
                "catalog_name": catalog_name,
                "population": population,
                "posterior_source": source,
                "n": len(grp),
                "median_e50": grp["e50"].median(),
                "median_zeta": grp["zeta_median"].median(),
                "median_ld_abs_max_delta": grp["ld_abs_max_delta"].median(),
                "frac_ld_abs_delta_gt_0p10": float((grp["ld_abs_max_delta"] > 0.10).mean()),
                "median_leverage": grp.get("delta_loglike_map_minus_min", pd.Series(dtype=float)).median(),
            }
        )
    return pd.DataFrame(rows).sort_values(["catalog_name", "posterior_source", "population"])


def main() -> None:
    args = parse_args()
    cfg = load_config(repo_root() / "config.json")
    out_dir = repo_root() / "outputs"
    out_dir.mkdir(exist_ok=True)

    ref = load_reference_catalog(cfg)
    all_matches = []
    for name, path in catalog_paths(cfg):
        matched = match_catalog_to_reference(name, path, ref)
        if not matched.empty:
            all_matches.append(matched)
    matches = pd.concat(all_matches, ignore_index=True) if all_matches else pd.DataFrame()

    tag = args.out_tag
    matches_path = out_dir / f"limb_darkening_catalog_offsets_{tag}.csv"
    summary_path = out_dir / f"limb_darkening_catalog_offset_summary_{tag}.csv"
    context_path = out_dir / f"limb_darkening_population_context_{tag}.csv"
    pop_path = out_dir / f"limb_darkening_population_summary_{tag}.csv"
    top_path = out_dir / f"limb_darkening_top_leverage_{tag}.csv"

    matches.to_csv(matches_path, index=False)
    offset_summary = summarize_offsets(matches)
    offset_summary.to_csv(summary_path, index=False)

    context = attach_population_context(matches, resolve(args.summary), resolve(args.leverage))
    context.to_csv(context_path, index=False)
    pop_summary = summarize_by_population(context)
    pop_summary.to_csv(pop_path, index=False)

    if not context.empty and "delta_loglike_map_minus_min" in context.columns:
        top = context.sort_values("delta_loglike_map_minus_min", ascending=False).head(100)
        top[
            [
                "kepoi_name",
                "koi_target",
                "population",
                "posterior_source",
                "e50",
                "zeta_median",
                "delta_loglike_map_minus_min",
                "catalog_name",
                "limbdark_1",
                "limbdark_2",
                "ref_limbdark_1",
                "ref_limbdark_2",
                "ld_du1",
                "ld_du2",
                "ld_abs_max_delta",
            ]
        ].to_csv(top_path, index=False)

    print("Wrote:")
    for path in [matches_path, summary_path, context_path, pop_path, top_path]:
        print(f"  {path}")
    print()
    print("Offset summary:")
    print(offset_summary.to_string(index=False))
    if not pop_summary.empty:
        print()
        print("Population summary:")
        print(pop_summary.to_string(index=False))


if __name__ == "__main__":
    main()
