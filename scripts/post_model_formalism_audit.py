from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

from common import load_config, normalize_dynesty_weights, output_dir, root_path
from extract_eccentricity_posteriors import absolute_density_error, stable_seed


G_SI = 6.67430e-11
RHO_SUN_KG_M3 = 1408.0
DAY_S = 86400.0


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Audit whether the zeta-KDE eccentricity extractor agrees with a direct "
            "MacDougall/Sagear-style post-model importance sampler for high-leverage planets."
        )
    )
    parser.add_argument("--config", default=None)
    parser.add_argument("--summary", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--leverage", default=None)
    parser.add_argument("--run-id", default="sagear_missing")
    parser.add_argument("--top-per-pop", type=int, default=15)
    parser.add_argument("--n-importance", type=int, default=150_000)
    parser.add_argument("--e-max", type=float, default=0.95)
    parser.add_argument("--out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--plot-out", default=None)
    args = parser.parse_args()

    out_dir = output_dir()
    summary_path = Path(args.summary) if args.summary else out_dir / "eccentricity_posterior_summary_merged_paired_exact.csv"
    sample_path = Path(args.sample) if args.sample else out_dir / "canonical_sample_old_astropy_rawcc.csv"
    leverage_path = (
        Path(args.leverage) if args.leverage else out_dir / "rayleigh_per_planet_leverage_PAIRED_EXACT_qcprimary.csv"
    )
    out_path = Path(args.out) if args.out else out_dir / "post_model_formalism_audit.csv"
    summary_out = Path(args.summary_out) if args.summary_out else out_dir / "post_model_formalism_audit_summary.csv"
    plot_out = Path(args.plot_out) if args.plot_out else out_dir / "post_model_formalism_audit_current_vs_direct.png"

    cfg = load_config(args.config)
    project = root_path(cfg, "alderaan_project")
    if project is None:
        raise ValueError("alderaan_project must be configured for FITS-level formalism audit")

    summary = pd.read_csv(summary_path)
    sample = pd.read_csv(sample_path)
    leverage = pd.read_csv(leverage_path)

    audited_ids = top_leverage_ids(leverage, args.top_per_pop)
    rows = []
    for kepoi_name in audited_ids:
        row = build_base_row(kepoi_name, leverage, summary, sample)
        if row.get("posterior_source") != "new_alderaan":
            row["direct_status"] = "skipped_existing_archive_no_new_fits"
            rows.append(row)
            continue
        target = str(row["koi_target"])
        fits_path = project / "Results" / args.run_id / target / f"{target}-results.fits"
        if not fits_path.exists():
            row["direct_status"] = "missing_new_fits"
            row["fits_path"] = str(fits_path)
            rows.append(row)
            continue
        try:
            direct = direct_importance_for_planet(row, fits_path, args.n_importance, args.e_max)
            row.update(direct)
        except Exception as exc:  # keep audit robust: one bad FITS should not hide the pattern
            row["direct_status"] = f"error:{type(exc).__name__}:{exc}"
        rows.append(row)

    audit = pd.DataFrame(rows)
    audit = add_interpretation_flags(audit)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(out_path, index=False)

    group_summary = summarize_audit(audit)
    group_summary.to_csv(summary_out, index=False)
    make_plot(audit, plot_out)

    print(f"Wrote: {out_path}")
    print(f"Wrote: {summary_out}")
    print(f"Wrote: {plot_out}")
    print("\nFormalism audit summary:")
    print(group_summary.to_string(index=False))


def top_leverage_ids(leverage: pd.DataFrame, top_per_pop: int) -> list[str]:
    ids: list[str] = []
    for _, sub in leverage.sort_values("delta_loglike_map_minus_min", ascending=False).groupby("population", sort=False):
        for kepoi_name in sub.head(top_per_pop)["kepoi_name"].astype(str):
            if kepoi_name not in ids:
                ids.append(kepoi_name)
    return ids


def build_base_row(kepoi_name: str, leverage: pd.DataFrame, summary: pd.DataFrame, sample: pd.DataFrame) -> dict:
    lev = leverage[leverage["kepoi_name"].astype(str).eq(kepoi_name)].iloc[0].to_dict()
    summ = summary[summary["kepoi_name"].astype(str).eq(kepoi_name)]
    samp = sample[sample["kepoi_name"].astype(str).eq(kepoi_name)]
    row = dict(lev)
    if not summ.empty:
        for col in [
            "posterior_source",
            "impact_mode",
            "posterior_file",
            "koi_period",
            "alderaan_planet_index",
            "alderaan_period_days",
            "period_relative_difference",
            "e16",
            "e84",
            "zeta_p16",
            "zeta_p84",
            "duration_ratio",
            "catalog_duration_hr",
            "alderaan_duration_median_hr",
            "ror_ratio",
            "catalog_ror",
            "alderaan_ror_median",
            "catalog_impact",
            "impact_fit_p50",
            "impact_used_p50",
            "qc_reasons",
        ]:
            if col in summ.columns:
                row[col] = summ.iloc[0].get(col)
    if not samp.empty:
        for col in [
            "rho_log",
            "rho_log_upper",
            "rho_log_lower",
            "koi_duration",
            "koi_ror",
            "koi_impact",
            "koi_model_snr",
            "koi_num_transits",
            "koi_prad",
            "berger_rad",
            "berger_mass",
            "berger_age",
            "berger_feh",
        ]:
            if col in samp.columns:
                row[f"sample_{col}"] = samp.iloc[0].get(col)
    row["e_min_from_zeta_median"] = e_min_from_zeta(row.get("zeta_median", np.nan))
    return row


def e_min_from_zeta(zeta: object) -> float:
    try:
        g = float(zeta)
    except Exception:
        return np.nan
    if not np.isfinite(g) or g <= 0:
        return np.nan
    return float(abs(1.0 - g * g) / (1.0 + g * g))


def direct_importance_for_planet(row: dict, fits_path: Path, n_importance: int, e_max: float) -> dict:
    planet_index = int(row["alderaan_planet_index"])
    suffix = f"_{planet_index}"
    with fits.open(fits_path, memmap=False) as hdul:
        samples = pd.DataFrame(np.array(hdul["SAMPLES"].data).byteswap().newbyteorder())
        npl = int(hdul[0].header.get("NPL", -1))
        ttimes_name = f"TTIMES_{planet_index:02d}"
        n_ttimes = len(hdul[ttimes_name].data) if ttimes_name in hdul else np.nan
        n_out_flag = int(np.sum(hdul[ttimes_name].data["OUT_FLAG"] > 0)) if ttimes_name in hdul else np.nan

    needed = [f"ROR{suffix}", f"IMPACT{suffix}", f"DUR14{suffix}", "LN_WT"]
    missing = [c for c in needed if c not in samples.columns]
    if missing:
        return {"direct_status": f"missing_columns:{','.join(missing)}", "fits_path": str(fits_path)}

    weights = normalize_dynesty_weights(samples["LN_WT"].to_numpy(float))
    rng = np.random.default_rng(stable_seed("direct-formalism", row["kepoi_name"], n_importance))
    sample_idx = rng.choice(np.arange(len(samples)), size=n_importance, replace=True, p=weights)
    s = samples.iloc[sample_idx]

    period_days = float(row.get("koi_period", row.get("alderaan_period_days")))
    period_s = period_days * DAY_S
    dur_s = s[f"DUR14{suffix}"].to_numpy(float) * DAY_S
    ror = s[f"ROR{suffix}"].to_numpy(float)
    impact = s[f"IMPACT{suffix}"].to_numpy(float)
    e = rng.uniform(0.0, e_max, size=n_importance)
    omega = rng.uniform(-0.5 * np.pi, 1.5 * np.pi, size=n_importance)

    rho_true = 10.0 ** float(row["sample_rho_log"])
    err_hi = absolute_density_error(row.get("sample_rho_log_upper", np.nan))
    err_lo = absolute_density_error(row.get("sample_rho_log_lower", np.nan))
    if not np.isfinite(err_hi) or err_hi <= 0:
        err_hi = 0.13 * rho_true
    if not np.isfinite(err_lo) or err_lo <= 0:
        err_lo = 0.13 * rho_true

    rho_model = rho_circ_sample(period_s, dur_s, ror, impact, e, omega)
    rho_circ_e0 = rho_circ_sample(period_s, dur_s, ror, impact, np.zeros_like(e), np.zeros_like(e))
    sigma = np.where(rho_model >= rho_true, err_hi, err_lo)
    resid = (rho_model - rho_true) / sigma
    loglike = -0.5 * resid**2
    valid = (
        np.isfinite(loglike)
        & np.isfinite(e)
        & np.isfinite(rho_model)
        & np.isfinite(sigma)
        & (sigma > 0)
        & (rho_model > 0)
    )
    if int(valid.sum()) < 100:
        return {
            "direct_status": f"too_few_valid_importance_samples:{int(valid.sum())}",
            "fits_path": str(fits_path),
            "direct_valid_fraction": float(valid.mean()),
        }

    e_valid = e[valid]
    loglike = loglike[valid]
    loglike -= np.nanmax(loglike)
    w = np.exp(loglike)
    w_sum = np.sum(w)
    if not np.isfinite(w_sum) or w_sum <= 0:
        return {"direct_status": "zero_importance_weight", "fits_path": str(fits_path)}
    w /= w_sum
    q16, q50, q84 = weighted_quantile(e_valid, w, [0.16, 0.5, 0.84])
    ess = 1.0 / np.sum(w**2)
    rho_circ_finite = rho_circ_e0[np.isfinite(rho_circ_e0) & (rho_circ_e0 > 0)]
    rho_circ_q16, rho_circ_q50, rho_circ_q84 = (
        np.nanpercentile(rho_circ_finite, [16, 50, 84]) if len(rho_circ_finite) else (np.nan, np.nan, np.nan)
    )

    return {
        "direct_status": "ok",
        "fits_path": str(fits_path),
        "fits_n_samples": int(len(samples)),
        "fits_npl": npl,
        "fits_n_ttimes": n_ttimes,
        "fits_n_ttimes_out_flag": n_out_flag,
        "direct_valid_fraction": float(valid.mean()),
        "direct_effective_sample_size": float(ess),
        "direct_e16": float(q16),
        "direct_e50": float(q50),
        "direct_e84": float(q84),
        "direct_rho_true_solar": float(rho_true),
        "direct_rho_err_hi_solar": float(err_hi),
        "direct_rho_err_lo_solar": float(err_lo),
        "rho_circ_e0_p16_solar": float(rho_circ_q16),
        "rho_circ_e0_p50_solar": float(rho_circ_q50),
        "rho_circ_e0_p84_solar": float(rho_circ_q84),
        "rho_circ_e0_over_true": float(rho_circ_q50 / rho_true) if np.isfinite(rho_circ_q50) else np.nan,
    }


def rho_circ_sample(
    period_s: float,
    dur_s: np.ndarray,
    ror: np.ndarray,
    impact: np.ndarray,
    e: np.ndarray,
    omega: np.ndarray,
) -> np.ndarray:
    numerator = (1.0 + ror) ** 2 - impact**2
    g = (1.0 + e * np.sin(omega)) / np.sqrt(1.0 - e**2)
    arg = (dur_s * np.pi / period_s) * g
    sin2 = np.sin(arg) ** 2
    base = numerator / sin2 + impact**2
    rho = (3.0 * np.pi / (G_SI * period_s**2)) * np.power(base, 1.5)
    invalid = (numerator <= 0) | (sin2 <= 0) | (base <= 0)
    rho = rho / RHO_SUN_KG_M3
    rho[invalid] = np.nan
    return rho


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: list[float]) -> list[float]:
    order = np.argsort(values)
    v = np.asarray(values)[order]
    w = np.asarray(weights)[order]
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    return [float(np.interp(q, cdf, v)) for q in quantiles]


def add_interpretation_flags(audit: pd.DataFrame) -> pd.DataFrame:
    out = audit.copy()
    out["current_minus_direct_e50"] = out["e50"] - out.get("direct_e50", np.nan)
    flags = []
    for _, row in out.iterrows():
        f: list[str] = []
        if row.get("posterior_source") == "existing_archive":
            f.append("archive_row_not_reaudited_from_fits")
        if bool(pd.notna(row.get("direct_e50"))) and row.get("e50", np.nan) - row.get("direct_e50", np.nan) > 0.10:
            f.append("zeta_kde_higher_than_direct_importance")
        if bool(pd.notna(row.get("direct_e50"))) and row.get("direct_e50", np.nan) > 0.20:
            f.append("direct_importance_still_high_e")
        if abs(float(row.get("duration_ratio", np.nan)) - 1.0) > 0.25 if pd.notna(row.get("duration_ratio", np.nan)) else False:
            f.append("large_duration_shift")
        if abs(float(row.get("ror_ratio", np.nan)) - 1.0) > 0.35 if pd.notna(row.get("ror_ratio", np.nan)) else False:
            f.append("large_ror_shift")
        if float(row.get("impact_used_p50", np.nan)) > 0.85 if pd.notna(row.get("impact_used_p50", np.nan)) else False:
            f.append("grazing_or_supergrazing_fit")
        rho_ratio = float(row.get("rho_circ_e0_over_true", np.nan)) if pd.notna(row.get("rho_circ_e0_over_true", np.nan)) else np.nan
        if np.isfinite(rho_ratio) and (rho_ratio < 0.7 or rho_ratio > 1.3):
            f.append("circular_transit_density_vs_berger_mismatch")
        zeta = float(row.get("zeta_median", np.nan)) if pd.notna(row.get("zeta_median", np.nan)) else np.nan
        if np.isfinite(zeta) and (zeta < 0.7 or zeta > 1.3):
            f.append("zeta_tail")
        if float(row.get("period_relative_difference", np.nan)) > 1e-3 if pd.notna(row.get("period_relative_difference", np.nan)) else False:
            f.append("period_mismatch")
        flags.append(";".join(f))
    out["audit_flags"] = flags
    return out


def summarize_audit(audit: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pop, sub in audit.groupby("population", dropna=False):
        ok = sub[sub["direct_status"].eq("ok")]
        rows.append(
            {
                "population": pop,
                "audited_rows": int(len(sub)),
                "direct_ok_rows": int(len(ok)),
                "new_alderaan_rows": int(sub["posterior_source"].eq("new_alderaan").sum()),
                "existing_archive_rows": int(sub["posterior_source"].eq("existing_archive").sum()),
                "median_current_e50": float(np.nanmedian(sub["e50"])) if len(sub) else np.nan,
                "median_direct_e50": float(np.nanmedian(ok["direct_e50"])) if len(ok) else np.nan,
                "median_current_minus_direct_e50": float(np.nanmedian(ok["current_minus_direct_e50"])) if len(ok) else np.nan,
                "direct_high_e_rows_gt_0p2": int((ok["direct_e50"] > 0.2).sum()) if len(ok) else 0,
                "zeta_kde_higher_than_direct_gt_0p1": int((ok["current_minus_direct_e50"] > 0.1).sum()) if len(ok) else 0,
                "median_direct_ess": float(np.nanmedian(ok["direct_effective_sample_size"])) if len(ok) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def make_plot(audit: pd.DataFrame, plot_out: Path) -> None:
    ok = audit[audit["direct_status"].eq("ok")].copy()
    if ok.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 6))
    colors = {
        "thin_singles": "#0b8b8c",
        "thick_singles": "#a30d0d",
        "thin_multis": "#5bbfc0",
        "thick_multis": "#d27a7a",
    }
    for pop, sub in ok.groupby("population"):
        ax.scatter(
            sub["e50"],
            sub["direct_e50"],
            s=42,
            alpha=0.85,
            label=pop,
            color=colors.get(pop, "0.4"),
            edgecolor="white",
            linewidth=0.5,
        )
    ax.plot([0, 0.95], [0, 0.95], color="0.4", lw=1, ls="--")
    ax.set_xlabel("current zeta-KDE e50")
    ax.set_ylabel("direct post-model importance e50")
    ax.set_title("High-Leverage Formalism Audit")
    ax.set_xlim(0, 0.95)
    ax.set_ylim(0, 0.95)
    ax.legend(frameon=True, fontsize=9)
    fig.tight_layout()
    plot_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_out, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
