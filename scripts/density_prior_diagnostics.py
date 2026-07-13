from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    load_config,
    output_dir,
    read_berger2018_stellar_table,
    read_berger_table2,
    root_path,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnose whether Berger+2018 radius-based density alternatives can explain "
            "the high-e photoeccentric tail."
        )
    )
    parser.add_argument("--config", default=None)
    parser.add_argument(
        "--posterior",
        default=None,
        help="Merged posterior summary. Defaults to paired-exact QC-primary output.",
    )
    parser.add_argument(
        "--sample",
        default=None,
        help="Canonical sample table used to add KOI reliability and planet-property fields.",
    )
    parser.add_argument(
        "--leverage",
        default=None,
        help="Per-planet leverage table. Defaults to paired-exact QC-primary leverage output.",
    )
    parser.add_argument(
        "--direct-audit",
        default=None,
        help="Optional post_model_formalism_audit.csv table.",
    )
    parser.add_argument("--tag", default="QC_PRIMARY")
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = output_dir()
    posterior_path = Path(args.posterior) if args.posterior else out / "eccentricity_posterior_summary_merged_paired_exact_qcprimary.csv"
    sample_path = Path(args.sample) if args.sample else out / "canonical_sample_old_astropy_rawcc.csv"
    leverage_path = Path(args.leverage) if args.leverage else out / "rayleigh_per_planet_leverage_PAIRED_EXACT_qcprimary.csv"
    direct_path = Path(args.direct_audit) if args.direct_audit else out / "post_model_formalism_audit.csv"

    if not posterior_path.exists():
        raise FileNotFoundError(f"Posterior summary not found: {posterior_path}")

    post = pd.read_csv(posterior_path)
    post["population"] = post["disk"].astype(str) + "_" + post["system"].astype(str) + "s"
    if sample_path.exists():
        sample_cols = [
            "kepoi_name",
            "koi_score",
            "koi_disposition",
            "koi_pdisposition",
            "koi_comment",
            "koi_fpflag_nt",
            "koi_fpflag_ss",
            "koi_fpflag_co",
            "koi_fpflag_ec",
            "koi_prad",
            "koi_model_snr",
            "koi_num_transits",
            "koi_duration",
            "koi_ror",
            "koi_impact",
        ]
        sample = pd.read_csv(sample_path, usecols=lambda c: c in sample_cols)
        sample = sample.drop_duplicates("kepoi_name")
        post = post.merge(sample, on="kepoi_name", how="left", suffixes=("", "_sample"))

    b20 = read_berger_table2(cfg)
    b18_path = root_path(cfg, "berger2018_stellar")
    if b18_path is None or not b18_path.exists():
        raise FileNotFoundError(f"Berger+2018 stellar table not found: {b18_path}")
    b18 = read_berger2018_stellar_table(b18_path)

    density = post.merge(
        b20[
            [
                "kepid",
                "berger_mass",
                "berger_rad",
                "rho_log",
                "rho_log_upper",
                "rho_log_lower",
                "berger_teff",
                "berger_logg",
                "berger_age",
            ]
        ],
        on="kepid",
        how="left",
        suffixes=("", "_b20"),
    ).merge(
        b18[
            [
                "kepid",
                "berger2018_teff",
                "berger2018_rad",
                "berger2018_rad_err_upper",
                "berger2018_rad_err_lower",
                "berger2018_evol",
                "berger2018_bin",
            ]
        ],
        on="kepid",
        how="left",
    )

    density["rho_current_solar"] = np.power(10.0, pd.to_numeric(density["rho_log"], errors="coerce"))
    density["rho_from_berger2020_mass_rad_solar"] = (
        pd.to_numeric(density["berger_mass"], errors="coerce")
        / np.power(pd.to_numeric(density["berger_rad"], errors="coerce"), 3)
    )
    density["rho_berger2018rad_berger2020mass_solar"] = (
        pd.to_numeric(density["berger_mass"], errors="coerce")
        / np.power(pd.to_numeric(density["berger2018_rad"], errors="coerce"), 3)
    )
    density["rho_berger2018rad_b20mass_upper_solar"] = (
        pd.to_numeric(density["berger_mass"], errors="coerce")
        / np.power(
            pd.to_numeric(density["berger2018_rad"], errors="coerce")
            - pd.to_numeric(density["berger2018_rad_err_lower"], errors="coerce"),
            3,
        )
    )
    density["rho_berger2018rad_b20mass_lower_solar"] = (
        pd.to_numeric(density["berger_mass"], errors="coerce")
        / np.power(
            pd.to_numeric(density["berger2018_rad"], errors="coerce")
            + pd.to_numeric(density["berger2018_rad_err_upper"], errors="coerce"),
            3,
        )
    )

    zeta = pd.to_numeric(density["zeta_median"], errors="coerce")
    density["rho_required_for_circular_over_current"] = np.where(
        np.isfinite(zeta) & (zeta > 0), np.power(zeta, -3), np.nan
    )
    density["rho_berger2020_massrad_over_table"] = (
        density["rho_from_berger2020_mass_rad_solar"] / density["rho_current_solar"]
    )
    density["rho_berger2018rad_b20mass_over_current"] = (
        density["rho_berger2018rad_berger2020mass_solar"] / density["rho_current_solar"]
    )
    density["rho_berger2018rad_b20mass_over_required"] = (
        density["rho_berger2018rad_b20mass_over_current"]
        / density["rho_required_for_circular_over_current"]
    )
    density["log_required_ratio"] = safe_log10(density["rho_required_for_circular_over_current"])
    density["log_b18like_ratio"] = safe_log10(density["rho_berger2018rad_b20mass_over_current"])
    density["density_shift_direction_matches_circular_need"] = direction_matches(
        density["rho_required_for_circular_over_current"],
        density["rho_berger2018rad_b20mass_over_current"],
    )
    density["b18like_density_within_factor_1p5_of_circular_need"] = within_factor(
        density["rho_berger2018rad_b20mass_over_current"],
        density["rho_required_for_circular_over_current"],
        1.5,
    )
    density["b18like_density_within_factor_2_of_circular_need"] = within_factor(
        density["rho_berger2018rad_b20mass_over_current"],
        density["rho_required_for_circular_over_current"],
        2.0,
    )

    if direct_path.exists():
        direct = pd.read_csv(direct_path)
        direct_cols = [
            "kepoi_name",
            "direct_status",
            "direct_e50",
            "direct_rho_true_solar",
            "rho_circ_e0_p50_solar",
            "rho_circ_e0_over_true",
            "audit_flags",
        ]
        direct_cols = [c for c in direct_cols if c in direct.columns]
        density = density.merge(direct[direct_cols], on="kepoi_name", how="left", suffixes=("", "_direct"))

    diagnostics_path = out / f"density_prior_reconstruction_diagnostics_{args.tag}.csv"
    density.to_csv(diagnostics_path, index=False)

    summary = summarize_density(density)
    summary_path = out / f"density_prior_reconstruction_summary_{args.tag}.csv"
    summary.to_csv(summary_path, index=False)

    top = top_table(density, leverage_path)
    top_path = out / f"density_prior_top_leverage_{args.tag}.csv"
    top.to_csv(top_path, index=False)

    plot_path = out / f"density_prior_required_vs_b18like_{args.tag}.png"
    make_plot(density, plot_path)

    print(f"Wrote {diagnostics_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {top_path}")
    print(f"Wrote {plot_path}")
    print("\nSummary:")
    print(summary.to_string(index=False))
    print("\nTop leverage density comparison:")
    print(top.head(20).to_string(index=False))


def safe_log10(values: pd.Series) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce")
    return pd.Series(np.where(np.isfinite(x) & (x > 0), np.log10(x), np.nan), index=values.index)


def direction_matches(required: pd.Series, candidate: pd.Series) -> pd.Series:
    req = pd.to_numeric(required, errors="coerce")
    cand = pd.to_numeric(candidate, errors="coerce")
    req_log = np.log(req.to_numpy(dtype=float))
    cand_log = np.log(cand.to_numpy(dtype=float))
    valid = np.isfinite(req_log) & np.isfinite(cand_log)
    out = np.zeros(len(req), dtype=bool)
    out[valid] = np.sign(req_log[valid]) == np.sign(cand_log[valid])
    near_circular = valid & (np.abs(req_log) < np.log(1.1))
    out[near_circular] = np.abs(cand_log[near_circular]) < np.log(1.1)
    return pd.Series(out, index=required.index)


def within_factor(candidate: pd.Series, target: pd.Series, factor: float) -> pd.Series:
    cand = pd.to_numeric(candidate, errors="coerce")
    targ = pd.to_numeric(target, errors="coerce")
    ratio = cand / targ
    valid = np.isfinite(ratio) & (ratio > 0)
    out = pd.Series(False, index=candidate.index)
    out[valid] = np.abs(np.log(ratio[valid])) <= math.log(factor)
    return out


def summarize_density(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for population, g in df.groupby("population", dropna=False):
        valid_b18 = np.isfinite(pd.to_numeric(g["rho_berger2018rad_b20mass_over_current"], errors="coerce"))
        valid_req = np.isfinite(pd.to_numeric(g["rho_required_for_circular_over_current"], errors="coerce"))
        valid_both = valid_b18 & valid_req
        rows.append(
            {
                "population": population,
                "n_planets": int(len(g)),
                "n_with_berger2018_radius": int(valid_b18.sum()),
                "median_e50": q(g["e50"], 0.50),
                "median_zeta": q(g["zeta_median"], 0.50),
                "median_required_rho_over_current": q(g["rho_required_for_circular_over_current"], 0.50),
                "p16_required_rho_over_current": q(g["rho_required_for_circular_over_current"], 0.16),
                "p84_required_rho_over_current": q(g["rho_required_for_circular_over_current"], 0.84),
                "median_b18like_rho_over_current": q(g["rho_berger2018rad_b20mass_over_current"], 0.50),
                "p16_b18like_rho_over_current": q(g["rho_berger2018rad_b20mass_over_current"], 0.16),
                "p84_b18like_rho_over_current": q(g["rho_berger2018rad_b20mass_over_current"], 0.84),
                "frac_need_density_factor_gt_1p5": frac(g["rho_required_for_circular_over_current"] > 1.5),
                "frac_need_density_factor_gt_2": frac(g["rho_required_for_circular_over_current"] > 2.0),
                "frac_b18like_direction_matches_need": frac(g.loc[valid_both, "density_shift_direction_matches_circular_need"]),
                "frac_b18like_within_factor_1p5": frac(g.loc[valid_both, "b18like_density_within_factor_1p5_of_circular_need"]),
                "frac_b18like_within_factor_2": frac(g.loc[valid_both, "b18like_density_within_factor_2_of_circular_need"]),
            }
        )
    rows.append(overall_row(df))
    return pd.DataFrame(rows)


def overall_row(df: pd.DataFrame) -> dict:
    valid_b18 = np.isfinite(pd.to_numeric(df["rho_berger2018rad_b20mass_over_current"], errors="coerce"))
    valid_req = np.isfinite(pd.to_numeric(df["rho_required_for_circular_over_current"], errors="coerce"))
    valid_both = valid_b18 & valid_req
    return {
        "population": "ALL",
        "n_planets": int(len(df)),
        "n_with_berger2018_radius": int(valid_b18.sum()),
        "median_e50": q(df["e50"], 0.50),
        "median_zeta": q(df["zeta_median"], 0.50),
        "median_required_rho_over_current": q(df["rho_required_for_circular_over_current"], 0.50),
        "p16_required_rho_over_current": q(df["rho_required_for_circular_over_current"], 0.16),
        "p84_required_rho_over_current": q(df["rho_required_for_circular_over_current"], 0.84),
        "median_b18like_rho_over_current": q(df["rho_berger2018rad_b20mass_over_current"], 0.50),
        "p16_b18like_rho_over_current": q(df["rho_berger2018rad_b20mass_over_current"], 0.16),
        "p84_b18like_rho_over_current": q(df["rho_berger2018rad_b20mass_over_current"], 0.84),
        "frac_need_density_factor_gt_1p5": frac(df["rho_required_for_circular_over_current"] > 1.5),
        "frac_need_density_factor_gt_2": frac(df["rho_required_for_circular_over_current"] > 2.0),
        "frac_b18like_direction_matches_need": frac(df.loc[valid_both, "density_shift_direction_matches_circular_need"]),
        "frac_b18like_within_factor_1p5": frac(df.loc[valid_both, "b18like_density_within_factor_1p5_of_circular_need"]),
        "frac_b18like_within_factor_2": frac(df.loc[valid_both, "b18like_density_within_factor_2_of_circular_need"]),
    }


def q(values: pd.Series, quantile: float) -> float:
    x = pd.to_numeric(values, errors="coerce")
    x = x[np.isfinite(x)]
    return float(np.quantile(x, quantile)) if len(x) else np.nan


def frac(mask: pd.Series | np.ndarray) -> float:
    if len(mask) == 0:
        return np.nan
    m = pd.Series(mask).dropna()
    return float(m.mean()) if len(m) else np.nan


def top_table(density: pd.DataFrame, leverage_path: Path) -> pd.DataFrame:
    if leverage_path.exists():
        lev = pd.read_csv(leverage_path)
        cols = [
            "kepoi_name",
            "population",
            "loglike_at_sigma_min",
            "loglike_at_sigma_map",
            "delta_loglike_map_minus_min",
        ]
        cols = [c for c in cols if c in lev.columns]
        merged = lev[cols].merge(density, on=["kepoi_name", "population"], how="left", suffixes=("", "_density"))
    else:
        merged = density.copy()
        merged["delta_loglike_map_minus_min"] = np.nan
    keep = [
        "kepoi_name",
        "kepid",
        "population",
        "posterior_source",
        "impact_mode",
        "e50",
        "zeta_median",
        "rho_required_for_circular_over_current",
        "rho_berger2018rad_b20mass_over_current",
        "rho_berger2018rad_b20mass_over_required",
        "berger_rad",
        "berger2018_rad",
        "berger2018_evol",
        "berger2018_bin",
        "koi_score",
        "koi_disposition",
        "koi_prad",
        "catalog_impact",
        "impact_fit_p50",
        "qc_reasons",
        "direct_e50",
        "rho_circ_e0_over_true",
        "audit_flags",
        "delta_loglike_map_minus_min",
    ]
    keep = [c for c in keep if c in merged.columns]
    merged = merged.sort_values("delta_loglike_map_minus_min", ascending=False, na_position="last")
    return merged[keep].head(120)


def make_plot(df: pd.DataFrame, path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"Skipping plot because matplotlib import failed: {exc}")
        return

    sub = df[
        np.isfinite(df["log_required_ratio"])
        & np.isfinite(df["log_b18like_ratio"])
        & df["population"].notna()
    ].copy()
    if sub.empty:
        return
    colors = {
        "thin_singles": "#008b8b",
        "thick_singles": "#990000",
        "thin_multis": "#66b2b2",
        "thick_multis": "#cc7777",
    }
    fig, ax = plt.subplots(figsize=(7, 6))
    for population, g in sub.groupby("population"):
        ax.scatter(
            g["log_required_ratio"],
            g["log_b18like_ratio"],
            s=14,
            alpha=0.45,
            label=population,
            color=colors.get(population, "0.4"),
            edgecolor="none",
        )
    lim = np.nanmax(np.abs(np.r_[sub["log_required_ratio"], sub["log_b18like_ratio"], [-1, 1]]))
    lim = min(max(lim, 0.5), 2.0)
    ax.plot([-lim, lim], [-lim, lim], color="black", lw=1.0, label="would exactly circularize")
    ax.axhline(0, color="0.5", lw=0.8, ls="--", label="no density change")
    ax.axvline(0, color="0.5", lw=0.8, ls=":")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log10(required rho / current rho) for circular zeta")
    ax.set_ylabel("log10(Berger2018-radius density / current rho)")
    ax.legend(fontsize=8, loc="best")
    ax.set_title("Can a Berger+2018-radius density shift explain high-e zeta?")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
