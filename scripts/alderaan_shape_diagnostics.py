from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

from common import load_config, normalize_dynesty_weights, output_dir, root_path
from extract_eccentricity_posteriors import period_from_ttimes


def main() -> None:
    parser = argparse.ArgumentParser(description="Diagnose ALDERAAN transit-shape inputs to photoeccentric inference.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--summary", default=None)
    parser.add_argument(
        "--results-dir",
        default=r"C:\Users\shres\Downloads\ALDERAAN_posteriors\ALDERAAN_posteriors",
        help="Flat directory containing ALDERAAN *-results.fits files.",
    )
    parser.add_argument("--max-planets", type=int, default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    sample_path = Path(args.sample) if args.sample else output_dir() / "canonical_sample_diagnostic.csv"
    summary_path = Path(args.summary) if args.summary else output_dir() / "eccentricity_posterior_summary.csv"
    sample = pd.read_csv(sample_path)
    summary = pd.read_csv(summary_path)
    if args.max_planets is not None:
        summary = summary.head(args.max_planets).copy()

    koi_cols = [
        "kepoi_name",
        "koi_duration",
        "koi_duration_err1",
        "koi_duration_err2",
        "koi_ror",
        "koi_impact",
        "koi_model_snr",
        "koi_prad",
        "rho_log",
        "rho_log_upper",
        "rho_log_lower",
        "P_thick",
    ]
    keep = sample[[c for c in koi_cols if c in sample.columns]].drop_duplicates("kepoi_name")
    df = summary.merge(keep, on="kepoi_name", how="left")

    rows = []
    results_dir = Path(args.results_dir)
    by_target = {k: v for k, v in df.groupby("koi_target", sort=False)}
    for i, (target, target_rows) in enumerate(by_target.items(), start=1):
        path = results_dir / f"{target}-results.fits"
        if not path.exists():
            continue
        rows.extend(process_target(path, target_rows))
        if i % 200 == 0:
            print(f"Processed {i}/{len(by_target)} targets")

    out = pd.DataFrame(rows)
    out_path = output_dir() / "alderaan_shape_diagnostics.csv"
    out.to_csv(out_path, index=False)

    flagged = flag_rows(out)
    flagged_path = output_dir() / "alderaan_shape_diagnostics_flagged.csv"
    flagged.to_csv(flagged_path, index=False)

    make_plots(out)
    write_markdown_summary(out, flagged)

    print(f"Wrote: {out_path}")
    print(f"Wrote: {flagged_path}")
    print("\nGroup summary:")
    print(group_summary(out).to_string(index=False))


def process_target(path: Path, target_rows: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    with fits.open(path, memmap=False) as hdul:
        samples = hdul["SAMPLES"].data
        weights = normalize_dynesty_weights(np.asarray(samples["LN_WT"], dtype=float))
        period_by_idx = {}
        for hdu in hdul:
            if hdu.name.startswith("TTIMES_"):
                idx = int(hdu.name.split("_")[1])
                period_by_idx[idx] = period_from_ttimes(hdu.data["TTIME"])

        for _, row in target_rows.iterrows():
            idx = int(row["alderaan_planet_index"])
            suffix = f"_{idx}"
            dur = np.asarray(samples[f"DUR14{suffix}"], dtype=float)
            ror = np.asarray(samples[f"ROR{suffix}"], dtype=float)
            impact = np.asarray(samples[f"IMPACT{suffix}"], dtype=float)
            finite = np.isfinite(dur) & np.isfinite(ror) & np.isfinite(impact)
            if not finite.any():
                continue
            w = weights[finite]
            w = w / w.sum()
            dur = dur[finite]
            ror = ror[finite]
            impact = impact[finite]

            dur_q = weighted_quantile(dur * 24.0, w, [0.16, 0.5, 0.84])
            ror_q = weighted_quantile(ror, w, [0.16, 0.5, 0.84])
            impact_q = weighted_quantile(impact, w, [0.16, 0.5, 0.84])
            koi_duration_hr = float(row.get("koi_duration", np.nan))
            koi_ror = float(row.get("koi_ror", np.nan))
            koi_impact = float(row.get("koi_impact", np.nan))
            rho = 10.0 ** float(row.get("rho_log", np.nan))
            rho_err_hi = absolute_density_error(row.get("rho_log_upper", np.nan))
            rho_err_lo = absolute_density_error(row.get("rho_log_lower", np.nan))
            rho_frac = np.nan
            if np.isfinite(rho) and rho > 0:
                rho_frac = np.nanmean([rho_err_hi, rho_err_lo]) / rho

            rows.append(
                {
                    "koi_target": row["koi_target"],
                    "kepoi_name": row["kepoi_name"],
                    "disk": row["disk"],
                    "system": row["system"],
                    "P_thick": row.get("P_thick", np.nan),
                    "koi_period": row["koi_period"],
                    "alderaan_period_days": period_by_idx.get(idx, np.nan),
                    "period_relative_difference": row["period_relative_difference"],
                    "koi_duration_hr": koi_duration_hr,
                    "alderaan_dur14_hr_p16": dur_q[0],
                    "alderaan_dur14_hr_med": dur_q[1],
                    "alderaan_dur14_hr_p84": dur_q[2],
                    "alderaan_to_koi_duration_ratio": dur_q[1] / koi_duration_hr if koi_duration_hr > 0 else np.nan,
                    "koi_ror": koi_ror,
                    "alderaan_ror_p16": ror_q[0],
                    "alderaan_ror_med": ror_q[1],
                    "alderaan_ror_p84": ror_q[2],
                    "alderaan_to_koi_ror_ratio": ror_q[1] / koi_ror if koi_ror > 0 else np.nan,
                    "koi_impact": koi_impact,
                    "alderaan_impact_p16": impact_q[0],
                    "alderaan_impact_med": impact_q[1],
                    "alderaan_impact_p84": impact_q[2],
                    "rho_solar": rho,
                    "rho_frac_err": rho_frac,
                    "e50": row["e50"],
                    "e16": row["e16"],
                    "e84": row["e84"],
                    "zeta_p16": row["zeta_p16"],
                    "zeta_median": row["zeta_median"],
                    "zeta_p84": row["zeta_p84"],
                    "n_zeta": row["n_zeta"],
                    "koi_model_snr": row.get("koi_model_snr", np.nan),
                    "koi_prad": row.get("koi_prad", np.nan),
                }
            )
    return rows


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: list[float]) -> np.ndarray:
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cdf = np.cumsum(w)
    cdf = cdf / cdf[-1]
    return np.interp(quantiles, cdf, v)


def absolute_density_error(value: object) -> float:
    value = float(value)
    if not np.isfinite(value):
        return np.nan
    return 10.0**value


def flag_rows(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["flag_short_zeta"] = out["zeta_median"] < 0.7
    out["flag_long_zeta"] = out["zeta_median"] > 1.3
    out["flag_duration_shift"] = np.abs(np.log(out["alderaan_to_koi_duration_ratio"])) > np.log(1.25)
    out["flag_ror_shift"] = np.abs(np.log(out["alderaan_to_koi_ror_ratio"])) > np.log(1.25)
    out["flag_high_impact"] = out["alderaan_impact_med"] > 0.85
    out["flag_large_rho_uncertainty"] = out["rho_frac_err"] > 0.35
    out["flag_high_e50"] = out["e50"] > 0.5
    flag_cols = [c for c in out.columns if c.startswith("flag_")]
    out["n_flags"] = out[flag_cols].sum(axis=1)
    return out[out["n_flags"] > 0].sort_values(["n_flags", "e50"], ascending=False)


def group_summary(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["disk", "system"])
        .agg(
            n=("kepoi_name", "count"),
            median_zeta=("zeta_median", "median"),
            q16_zeta=("zeta_median", lambda x: x.quantile(0.16)),
            q84_zeta=("zeta_median", lambda x: x.quantile(0.84)),
            median_e50=("e50", "median"),
            median_duration_ratio=("alderaan_to_koi_duration_ratio", "median"),
            median_rho_frac_err=("rho_frac_err", "median"),
            high_e50_frac=("e50", lambda x: float((x > 0.5).mean())),
            short_zeta_frac=("zeta_median", lambda x: float((x < 0.7).mean())),
            long_zeta_frac=("zeta_median", lambda x: float((x > 1.3).mean())),
        )
        .reset_index()
    )


def make_plots(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    colors = {("thin", "single"): "#008b8b", ("thick", "single"): "#8b0000", ("thin", "multi"): "#66c2c2", ("thick", "multi"): "#c46a6a"}
    for (disk, system), sub in df.groupby(["disk", "system"]):
        label = f"{disk} {system}"
        color = colors.get((disk, system), None)
        axes[0].hist(sub["zeta_median"], bins=np.linspace(0, 2.5, 80), histtype="step", density=True, color=color, label=label)
        axes[1].hist(sub["e50"], bins=np.linspace(0, 0.95, 80), histtype="step", density=True, color=color, label=label)
        axes[2].scatter(sub["zeta_median"], sub["e50"], s=8, alpha=0.35, color=color, label=label)
    axes[0].axvline(1.0, color="0.3", lw=1, ls="--")
    axes[0].set_xlabel(r"median $\zeta = T_{obs}/T_{circ}$")
    axes[0].set_ylabel("density")
    axes[1].set_xlabel("posterior median e")
    axes[1].set_ylabel("density")
    axes[2].set_xlabel(r"median $\zeta$")
    axes[2].set_ylabel("posterior median e")
    axes[2].set_xlim(0, 2.5)
    axes[2].set_ylim(0, 0.95)
    axes[0].legend(fontsize=8)
    fig.savefig(output_dir() / "alderaan_shape_diagnostics.png", dpi=200)
    plt.close(fig)


def write_markdown_summary(df: pd.DataFrame, flagged: pd.DataFrame) -> None:
    path = output_dir() / "alderaan_shape_diagnostics.md"
    lines = [
        "# ALDERAAN Shape Diagnostics",
        "",
        "This diagnostic compares KOI catalog transit-shape quantities to the existing ALDERAAN FITS posteriors used in the current photoeccentric extraction.",
        "",
        "## Group Summary",
        "",
        group_summary(df).to_markdown(index=False, floatfmt=".4f"),
        "",
        "## Most Flagged Planets",
        "",
        flagged[
            [
                "kepoi_name",
                "disk",
                "system",
                "n_flags",
                "e50",
                "zeta_median",
                "alderaan_to_koi_duration_ratio",
                "alderaan_impact_med",
                "rho_frac_err",
            ]
        ]
        .head(30)
        .to_markdown(index=False, floatfmt=".4f"),
        "",
        "## Interpretation Prompt",
        "",
        "If the high-e thin-single signal is real in these posteriors, it should appear as a coherent short-duration zeta tail rather than only catalog/ALDERAAN duration mismatches, extreme impact parameters, or very broad stellar-density priors.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
