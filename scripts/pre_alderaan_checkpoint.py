from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import output_dir


SAGEAR_E = {
    "thin_singles": 0.022,
    "thick_singles": 0.066,
    "thin_multis": 0.030,
    "thick_multis": 0.033,
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the current best pre-ALDERAAN checkpoint products.")
    parser.add_argument("--sample", default=None)
    parser.add_argument("--summary", default=None)
    parser.add_argument("--rayleigh", default=None)
    parser.add_argument("--shape", default=None)
    parser.add_argument("--prefix", default="pre_alderaan_best")
    args = parser.parse_args()

    out = output_dir()
    sample_path = Path(args.sample) if args.sample else out / "canonical_sample_old_astropy_rawcc.csv"
    summary_path = Path(args.summary) if args.summary else out / "eccentricity_posterior_summary_old_astropy_rawcc.csv"
    rayleigh_path = Path(args.rayleigh) if args.rayleigh else out / "rayleigh_population_fit_transit_selection_old_astropy_rawcc.csv"
    shape_path = Path(args.shape) if args.shape else out / "alderaan_shape_diagnostics.csv"

    sample = pd.read_csv(sample_path)
    summary = pd.read_csv(summary_path)
    rayleigh = pd.read_csv(rayleigh_path)
    shape = pd.read_csv(shape_path) if shape_path.exists() else pd.DataFrame()

    pop = population_table(sample, summary, rayleigh)
    pop_path = out / f"{args.prefix}_population_table.csv"
    pop.to_csv(pop_path, index=False)

    missing = missing_posteriors(sample, summary)
    missing_path = out / f"{args.prefix}_missing_posteriors.csv"
    missing.to_csv(missing_path, index=False)

    outliers = top_outliers(summary, shape)
    outlier_path = out / f"{args.prefix}_top_outliers.csv"
    outliers.to_csv(outlier_path, index=False)

    fig_path = out / f"{args.prefix}_eccentricity_diagnostics.png"
    make_eccentricity_plot(summary, fig_path)

    md_path = out / f"{args.prefix}_summary.md"
    write_markdown(md_path, pop, missing, outliers, sample_path, summary_path)

    print(pop.to_string(index=False))
    print(f"\nWrote: {pop_path}")
    print(f"Wrote: {missing_path}")
    print(f"Wrote: {outlier_path}")
    print(f"Wrote: {fig_path}")
    print(f"Wrote: {md_path}")


def population_table(sample: pd.DataFrame, summary: pd.DataFrame, rayleigh: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for disk, system, label in [
        ("thin", "single", "thin_singles"),
        ("thick", "single", "thick_singles"),
        ("thin", "multi", "thin_multis"),
        ("thick", "multi", "thick_multis"),
    ]:
        sub = sample[(sample["disk"] == disk) & (sample["system"] == system)]
        post = summary[(summary["disk"] == disk) & (summary["system"] == system)]
        fit = rayleigh[rayleigh["population"] == label].iloc[0]
        rows.append(
            {
                "population": label,
                "sample_planets": int(len(sub)),
                "sample_hosts": int(sub["kepid"].nunique()),
                "posterior_planets": int(len(post)),
                "missing_posteriors": int(len(sub) - len(post)),
                "coverage_fraction": len(post) / len(sub) if len(sub) else np.nan,
                "median_e50": float(np.nanmedian(post["e50"])) if len(post) else np.nan,
                "frac_e50_gt_0p5": float((post["e50"] > 0.5).mean()) if len(post) else np.nan,
                "rayleigh_expected_e": float(fit["expected_e"]),
                "rayleigh_e_lo": float(fit["expected_e_lo"]),
                "rayleigh_e_hi": float(fit["expected_e_hi"]),
                "sagear_rayleigh_e": SAGEAR_E[label],
                "delta_rayleigh_e": float(fit["expected_e"]) - SAGEAR_E[label],
            }
        )
    return pd.DataFrame(rows)


def missing_posteriors(sample: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    matched = set(summary["kepoi_name"].astype(str))
    missing = sample[~sample["kepoi_name"].astype(str).isin(matched)].copy()
    keep = [
        "kepid",
        "kepoi_name",
        "disk",
        "system",
        "koi_period",
        "koi_model_snr",
        "koi_prad",
        "koi_count",
        "P_thick",
        "system_after_all_cuts",
    ]
    return missing[[c for c in keep if c in missing.columns]].sort_values(["disk", "system", "koi_model_snr"], ascending=[True, True, False])


def top_outliers(summary: pd.DataFrame, shape: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "kepoi_name",
        "alderaan_to_koi_duration_ratio",
        "alderaan_to_koi_ror_ratio",
        "alderaan_impact_med",
        "rho_frac_err",
        "koi_model_snr",
    ]
    if not shape.empty:
        merged = summary.merge(shape[[c for c in cols if c in shape.columns]], on="kepoi_name", how="left")
    else:
        merged = summary.copy()
    merged["duration_shift_abs_log"] = np.abs(np.log(merged.get("alderaan_to_koi_duration_ratio", np.nan)))
    merged["ror_shift_abs_log"] = np.abs(np.log(merged.get("alderaan_to_koi_ror_ratio", np.nan)))
    merged["quality_flag_count"] = (
        (merged["zeta_median"] < 0.7).astype(int)
        + (merged["zeta_median"] > 1.3).astype(int)
        + (merged["duration_shift_abs_log"] > np.log(1.25)).astype(int)
        + (merged["ror_shift_abs_log"] > np.log(1.25)).astype(int)
        + (merged.get("alderaan_impact_med", 0) > 0.85).astype(int)
    )
    keep = [
        "kepid",
        "kepoi_name",
        "disk",
        "system",
        "koi_period",
        "e16",
        "e50",
        "e84",
        "zeta_median",
        "quality_flag_count",
        "alderaan_to_koi_duration_ratio",
        "alderaan_to_koi_ror_ratio",
        "alderaan_impact_med",
        "koi_model_snr",
        "posterior_file",
    ]
    return merged[[c for c in keep if c in merged.columns]].sort_values(["system", "e50"], ascending=[False, False]).head(80)


def make_eccentricity_plot(summary: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    colors = {
        ("thin", "single"): "#008b8b",
        ("thick", "single"): "#9b1b1b",
        ("thin", "multi"): "#59bfc1",
        ("thick", "multi"): "#c75d5d",
    }
    labels = {
        ("thin", "single"): "thin singles",
        ("thick", "single"): "thick singles",
        ("thin", "multi"): "thin multis",
        ("thick", "multi"): "thick multis",
    }
    bins = np.linspace(0, 0.95, 70)
    for key, sub in summary.groupby(["disk", "system"]):
        color = colors.get(key)
        label = labels.get(key, f"{key[0]} {key[1]}")
        e = np.sort(sub["e50"].dropna().to_numpy())
        if len(e) == 0:
            continue
        cdf = np.arange(1, len(e) + 1) / len(e)
        axes[0, 0].hist(e, bins=bins, histtype="step", density=True, color=color, label=label, linewidth=1.8)
        axes[0, 1].plot(e, cdf, color=color, label=label, linewidth=1.8)
        axes[1, 0].hist(sub["zeta_median"], bins=np.linspace(0, 2.4, 70), histtype="step", density=True, color=color, linewidth=1.8)
        axes[1, 1].scatter(sub["zeta_median"], sub["e50"], s=10, alpha=0.35, color=color, label=label)
    axes[0, 0].set_xlabel("posterior median e")
    axes[0, 0].set_ylabel("density")
    axes[0, 1].set_xlabel("posterior median e")
    axes[0, 1].set_ylabel("CDF")
    axes[1, 0].set_xlabel(r"median $\zeta$")
    axes[1, 0].set_ylabel("density")
    axes[1, 1].set_xlabel(r"median $\zeta$")
    axes[1, 1].set_ylabel("posterior median e")
    axes[1, 1].set_xlim(0, 2.4)
    axes[1, 1].set_ylim(0, 0.95)
    axes[0, 0].legend(fontsize=8)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_markdown(path: Path, pop: pd.DataFrame, missing: pd.DataFrame, outliers: pd.DataFrame, sample_path: Path, summary_path: Path) -> None:
    lines = [
        "# Pre-ALDERAAN Checkpoint",
        "",
        f"Sample: `{sample_path}`",
        f"Posterior summary: `{summary_path}`",
        "",
        "## Population Table",
        "",
        pop.to_markdown(index=False, floatfmt=".4f"),
        "",
        "## Missing Posterior Counts",
        "",
        missing.groupby(["disk", "system"]).size().reset_index(name="missing").to_markdown(index=False),
        "",
        "## Highest-e Existing Posterior Rows",
        "",
        outliers.head(30).to_markdown(index=False, floatfmt=".4f"),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
