from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import load_config, output_dir, root_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit mixed posterior sources and identify non-Sagear-equivalent archive rows."
    )
    parser.add_argument("--config", default=None)
    parser.add_argument(
        "--posterior",
        default=None,
        help="Merged posterior summary. Defaults to paired-exact QC-primary output.",
    )
    parser.add_argument(
        "--leverage",
        default=None,
        help="Per-planet leverage table. Defaults to paired-exact QC-primary leverage output.",
    )
    parser.add_argument("--tag", default="QC_PRIMARY")
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = output_dir()
    posterior_path = Path(args.posterior) if args.posterior else out / "eccentricity_posterior_summary_merged_paired_exact_qcprimary.csv"
    leverage_path = Path(args.leverage) if args.leverage else out / "rayleigh_per_planet_leverage_PAIRED_EXACT_qcprimary.csv"
    project = root_path(cfg, "alderaan_project")
    if project is None:
        project = Path(cfg["_root"]) / "sagear_reproduction" / "alderaan_project"

    post = pd.read_csv(posterior_path)
    post["population"] = post["disk"].astype(str) + "_" + post["system"].astype(str) + "s"
    post["rho_required_for_circular_over_current"] = np.where(
        np.isfinite(pd.to_numeric(post["zeta_median"], errors="coerce"))
        & (pd.to_numeric(post["zeta_median"], errors="coerce") > 0),
        np.power(pd.to_numeric(post["zeta_median"], errors="coerce"), -3),
        np.nan,
    )
    post["has_raw_alderaan_fits_in_project"] = post["koi_target"].map(
        lambda target: result_file_exists(project, str(target))
    )
    post["needs_raw_alderaan_to_be_sagear_equivalent"] = (
        (post["posterior_source"] == "existing_archive") | (post["impact_mode"] != "alderaan")
    ) & ~post["has_raw_alderaan_fits_in_project"]

    source_summary = summarize_sources(post)
    source_summary_path = out / f"posterior_source_population_summary_{args.tag}.csv"
    source_summary.to_csv(source_summary_path, index=False)

    source_counts = (
        post.groupby(["population", "posterior_source", "impact_mode"], dropna=False)
        .size()
        .reset_index(name="n_planets")
        .sort_values(["population", "posterior_source", "impact_mode"])
    )
    source_counts_path = out / f"posterior_source_population_counts_{args.tag}.csv"
    source_counts.to_csv(source_counts_path, index=False)

    target_manifest = archive_target_manifest(post)
    manifest_path = out / f"alderaan_needed_to_replace_archive_targets_{args.tag}.csv"
    target_manifest.to_csv(manifest_path, index=False)

    high_leverage = high_leverage_mix(post, leverage_path)
    leverage_path_out = out / f"posterior_source_high_leverage_mix_{args.tag}.csv"
    high_leverage.to_csv(leverage_path_out, index=False)

    print(f"Wrote {source_summary_path}")
    print(f"Wrote {source_counts_path}")
    print(f"Wrote {manifest_path}")
    print(f"Wrote {leverage_path_out}")
    print("\nSource summary:")
    print(source_summary.to_string(index=False))
    print("\nCounts:")
    print(source_counts.to_string(index=False))
    print("\nArchive targets needing raw ALDERAAN, top 20:")
    print(target_manifest.head(20).to_string(index=False))
    print("\nHigh leverage source mix, top 40:")
    print(high_leverage.head(40).to_string(index=False))


def result_file_exists(project: Path, target: str) -> bool:
    if not target or target == "nan":
        return False
    results_root = project / "Results"
    if not results_root.exists():
        return False
    return any(results_root.glob(f"*/{target}/{target}-results.fits"))


def summarize_sources(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (population, source, impact), g in df.groupby(["population", "posterior_source", "impact_mode"], dropna=False):
        rows.append(
            {
                "population": population,
                "posterior_source": source,
                "impact_mode": impact,
                "n_planets": int(len(g)),
                "n_targets": int(g["koi_target"].nunique()),
                "median_e50": q(g["e50"], 0.50),
                "p84_e50": q(g["e50"], 0.84),
                "median_zeta": q(g["zeta_median"], 0.50),
                "frac_e50_gt_0p3": frac(pd.to_numeric(g["e50"], errors="coerce") > 0.3),
                "frac_zeta_tail_review": frac(g["qc_reasons"].fillna("").str.contains("zeta_tail", regex=False)),
                "frac_e84_gt_0p9_review": frac(g["qc_reasons"].fillna("").str.contains("e84_gt_0p9", regex=False)),
                "frac_has_raw_alderaan_fits_in_project": frac(g["has_raw_alderaan_fits_in_project"]),
                "frac_needs_raw_alderaan_to_be_sagear_equivalent": frac(
                    g["needs_raw_alderaan_to_be_sagear_equivalent"]
                ),
                "median_required_rho_over_current": q(g["rho_required_for_circular_over_current"], 0.50),
            }
        )
    return pd.DataFrame(rows).sort_values(["population", "posterior_source", "impact_mode"])


def archive_target_manifest(df: pd.DataFrame) -> pd.DataFrame:
    archive = df[df["needs_raw_alderaan_to_be_sagear_equivalent"]].copy()
    if archive.empty:
        return pd.DataFrame()
    rows = []
    for target, g in archive.groupby("koi_target", dropna=False):
        rows.append(
            {
                "koi_target": target,
                "kepid": first(g, "kepid"),
                "n_archive_planets": int(len(g)),
                "populations": ";".join(sorted(map(str, g["population"].dropna().unique()))),
                "systems": ";".join(sorted(map(str, g["system"].dropna().unique()))),
                "max_e50": q(g["e50"], 1.0),
                "median_e50": q(g["e50"], 0.50),
                "min_zeta": q(g["zeta_median"], 0.0),
                "max_zeta": q(g["zeta_median"], 1.0),
                "n_zeta_tail_review": int(g["qc_reasons"].fillna("").str.contains("zeta_tail", regex=False).sum()),
                "n_e84_gt_0p9_review": int(g["qc_reasons"].fillna("").str.contains("e84_gt_0p9", regex=False).sum()),
                "has_raw_alderaan_fits_in_project": bool(g["has_raw_alderaan_fits_in_project"].any()),
                "priority_score": priority_score(g),
                "kepoi_names": ";".join(g["kepoi_name"].astype(str)),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(["priority_score", "max_e50", "n_archive_planets"], ascending=False)


def high_leverage_mix(post: pd.DataFrame, leverage_path: Path) -> pd.DataFrame:
    if not leverage_path.exists():
        return pd.DataFrame()
    lev = pd.read_csv(leverage_path)
    merged = lev.merge(
        post[
            [
                "kepoi_name",
                "posterior_source",
                "impact_mode",
                "qc_reasons",
                "has_raw_alderaan_fits_in_project",
                "needs_raw_alderaan_to_be_sagear_equivalent",
                "rho_required_for_circular_over_current",
            ]
        ],
        on="kepoi_name",
        how="left",
    )
    keep = [
        "kepoi_name",
        "koi_target",
        "kepid",
        "population",
        "e50",
        "zeta_median",
        "posterior_source",
        "impact_mode",
        "delta_loglike_map_minus_min",
        "qc_reasons",
        "has_raw_alderaan_fits_in_project",
        "needs_raw_alderaan_to_be_sagear_equivalent",
        "rho_required_for_circular_over_current",
    ]
    keep = [c for c in keep if c in merged.columns]
    return merged.sort_values("delta_loglike_map_minus_min", ascending=False)[keep].head(250)


def priority_score(g: pd.DataFrame) -> float:
    score = 0.0
    score += float(np.nanmax(pd.to_numeric(g["e50"], errors="coerce"))) * 10.0
    score += float(g["qc_reasons"].fillna("").str.contains("zeta_tail", regex=False).sum()) * 2.0
    if (g["system"] == "single").any():
        score += 2.0
    if (g["disk"] == "thin").any():
        score += 1.0
    return score


def q(values: pd.Series, quantile: float) -> float:
    x = pd.to_numeric(values, errors="coerce")
    x = x[np.isfinite(x)]
    return float(np.quantile(x, quantile)) if len(x) else np.nan


def frac(mask: pd.Series | np.ndarray) -> float:
    if len(mask) == 0:
        return np.nan
    m = pd.Series(mask).dropna()
    return float(m.mean()) if len(m) else np.nan


def first(df: pd.DataFrame, col: str):
    vals = df[col].dropna()
    if vals.empty:
        return np.nan
    return vals.iloc[0]


if __name__ == "__main__":
    main()
