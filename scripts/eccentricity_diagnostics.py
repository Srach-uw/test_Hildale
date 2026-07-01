from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import load_config, output_dir, root_path, write_json


POPULATIONS = [
    ("thin", "single", "thin_singles"),
    ("thick", "single", "thick_singles"),
    ("thin", "multi", "thin_multis"),
    ("thick", "multi", "thick_multis"),
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Diagnose eccentricity distributions and outlier sensitivity.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--eccentricity-table", default=None)
    parser.add_argument("--out-prefix", default="eccentricity_diagnostics")
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = output_dir()
    sample_path = Path(args.sample) if args.sample else out_dir / "canonical_sample_diagnostic.csv"
    sample = pd.read_csv(sample_path)
    ecc, ecc_meta = load_eccentricities(cfg, args.eccentricity_table)
    merged = merge_eccentricities(sample, ecc)

    summary = population_summary(merged, cfg)
    sensitivity = outlier_sensitivity(merged, cfg)
    quality = quality_cut_sensitivity(merged)
    top = top_outliers(merged, disk="thin", system="single", n=40)
    overlap = label_overlap(sample, ecc)

    summary_path = out_dir / f"{args.out_prefix}_population_summary.csv"
    sensitivity_path = out_dir / f"{args.out_prefix}_outlier_sensitivity.csv"
    quality_path = out_dir / f"{args.out_prefix}_quality_sensitivity.csv"
    top_path = out_dir / f"{args.out_prefix}_thin_single_top_outliers.csv"
    overlap_path = out_dir / f"{args.out_prefix}_label_overlap.csv"
    plot_path = out_dir / f"{args.out_prefix}_thin_single_distribution.png"

    summary.to_csv(summary_path, index=False)
    sensitivity.to_csv(sensitivity_path, index=False)
    quality.to_csv(quality_path, index=False)
    top.to_csv(top_path, index=False)
    overlap.to_csv(overlap_path, index=False)
    make_distribution_plot(merged, cfg, plot_path)

    metadata = {
        "sample": str(sample_path),
        "eccentricity_source": ecc_meta,
        "merged_planets_with_e": int(merged["e_value"].notna().sum()),
        "warning": (
            "If eccentricity_source.kind is photoeccentric_point_estimate, these diagnostics are for sample/outlier "
            "triage only and are not a replacement for ALDERAAN posterior hierarchical inference."
        ),
        "outputs": {
            "population_summary": str(summary_path),
            "outlier_sensitivity": str(sensitivity_path),
            "quality_sensitivity": str(quality_path),
            "thin_single_top_outliers": str(top_path),
            "label_overlap": str(overlap_path),
            "plot": str(plot_path),
        },
    }
    write_json(out_dir / f"{args.out_prefix}_metadata.json", metadata)

    print("=== Eccentricity Source ===")
    print(ecc_meta)
    print("\n=== Population Summary ===")
    print(summary.to_string(index=False))
    print(f"\nWrote: {summary_path}")
    print(f"Wrote: {sensitivity_path}")
    print(f"Wrote: {quality_path}")
    print(f"Wrote: {top_path}")
    print(f"Wrote: {plot_path}")


def load_eccentricities(cfg: dict, explicit_path: str | None) -> tuple[pd.DataFrame, dict[str, str]]:
    root = Path(cfg["_root"])
    candidates = []
    if explicit_path:
        candidates.append(Path(explicit_path))
    candidates.extend(
        [
            output_dir() / "eccentricity_posterior_summary.csv",
            root / "kepler_disk_eccentricity_analysis.csv",
            root / "clean_planets_qc.csv",
            root / "merged_planets_full.csv",
        ]
    )

    for path in candidates:
        if not path.exists():
            continue
        df = pd.read_csv(path)
        cols = set(df.columns)
        if {"kepoi_name", "e50"}.issubset(cols):
            out = df[["kepoi_name", "e50", "e16", "e84"]].copy()
            out = out.rename(columns={"e50": "e_value"})
            out["e_source_kind"] = "alderaan_posterior_median"
            return out, {"kind": "alderaan_posterior_median", "path": str(path)}
        if {"kepoi_name", "e_photo"}.issubset(cols):
            keep = [c for c in ["kepoi_name", "kepid", "e_photo", "disk", "system"] if c in df.columns]
            out = df[keep].copy().rename(columns={"e_photo": "e_value", "disk": "old_disk", "system": "old_system"})
            out["e_source_kind"] = "photoeccentric_point_estimate"
            return out, {"kind": "photoeccentric_point_estimate", "path": str(path), "column": "e_photo"}
        if {"kepoi_name", "e_photoeccentric"}.issubset(cols):
            keep = [
                c
                for c in ["kepoi_name", "kepid", "e_photoeccentric", "disk", "system_type", "duration_ratio"]
                if c in df.columns
            ]
            out = df[keep].copy().rename(
                columns={
                    "e_photoeccentric": "e_value",
                    "disk": "old_disk",
                    "system_type": "old_system",
                }
            )
            out["e_source_kind"] = "photoeccentric_point_estimate"
            return out, {"kind": "photoeccentric_point_estimate", "path": str(path), "column": "e_photoeccentric"}
    raise FileNotFoundError("No supported eccentricity table found.")


def merge_eccentricities(sample: pd.DataFrame, ecc: pd.DataFrame) -> pd.DataFrame:
    cols_to_drop = [c for c in ["old_disk", "old_system"] if c in sample.columns]
    sample = sample.drop(columns=cols_to_drop)
    merged = sample.merge(ecc, on="kepoi_name", how="left", suffixes=("", "_ecc"))
    merged["e_value"] = pd.to_numeric(merged["e_value"], errors="coerce")
    merged.loc[(merged["e_value"] < 0) | (merged["e_value"] >= 1), "e_value"] = np.nan
    return merged


def population_summary(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    rows = []
    for disk, system, label in POPULATIONS:
        sub = finite_population(df, disk, system)
        e = sub["e_value"].to_numpy(float)
        sagear_e = cfg["sagear_targets"]["rayleigh_expected_e"].get(label)
        rows.append(
            {
                "population": label,
                "planets_with_e": int(len(e)),
                "hosts_with_e": int(sub["kepid"].nunique()),
                "mean_e": safe_stat(e, np.mean),
                "median_e": safe_stat(e, np.median),
                "q16_e": safe_percentile(e, 16),
                "q84_e": safe_percentile(e, 84),
                "q95_e": safe_percentile(e, 95),
                "max_e": safe_stat(e, np.max),
                "frac_e_gt_0p05": frac_gt(e, 0.05),
                "frac_e_gt_0p10": frac_gt(e, 0.10),
                "frac_e_gt_0p20": frac_gt(e, 0.20),
                "sagear_rayleigh_expected_e": sagear_e,
                "mean_minus_sagear": safe_stat(e, np.mean) - sagear_e if sagear_e is not None and len(e) else np.nan,
                "median_minus_sagear": safe_stat(e, np.median) - sagear_e if sagear_e is not None and len(e) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def outlier_sensitivity(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    rows = []
    for disk, system, label in POPULATIONS:
        sub = finite_population(df, disk, system).sort_values("e_value", ascending=False).reset_index(drop=True)
        e_all = sub["e_value"].to_numpy(float)
        removals = [0, 1, 3, 5, 10, 20, 50]
        if len(e_all):
            removals.extend(sorted({int(np.ceil(len(e_all) * f)) for f in [0.01, 0.03, 0.05, 0.10]}))
        for remove_top_n in sorted(set(r for r in removals if r < len(e_all))):
            e = e_all[remove_top_n:]
            rows.append(
                {
                    "population": label,
                    "rule": f"remove_top_{remove_top_n}_by_e",
                    "n_remaining": int(len(e)),
                    "mean_e": safe_stat(e, np.mean),
                    "median_e": safe_stat(e, np.median),
                    "q95_e": safe_percentile(e, 95),
                    "sagear_rayleigh_expected_e": cfg["sagear_targets"]["rayleigh_expected_e"].get(label),
                }
            )
        for threshold in [0.05, 0.10, 0.15, 0.20, 0.30]:
            e = e_all[e_all <= threshold]
            rows.append(
                {
                    "population": label,
                    "rule": f"keep_e_le_{threshold:.2f}",
                    "n_remaining": int(len(e)),
                    "mean_e": safe_stat(e, np.mean),
                    "median_e": safe_stat(e, np.median),
                    "q95_e": safe_percentile(e, 95),
                    "sagear_rayleigh_expected_e": cfg["sagear_targets"]["rayleigh_expected_e"].get(label),
                }
            )
    return pd.DataFrame(rows)


def quality_cut_sensitivity(df: pd.DataFrame) -> pd.DataFrame:
    base = finite_population(df, "thin", "single")
    rules: list[tuple[str, pd.Series]] = [("none", pd.Series(True, index=base.index))]
    numeric_rules = [
        ("impact_lt_0p9", "koi_impact", lambda s: s < 0.9),
        ("impact_lt_0p8", "koi_impact", lambda s: s < 0.8),
        ("snr_ge_10", "koi_model_snr", lambda s: s >= 10),
        ("snr_ge_20", "koi_model_snr", lambda s: s >= 20),
        ("snr_ge_30", "koi_model_snr", lambda s: s >= 30),
        ("num_transits_ge_3", "koi_num_transits", lambda s: s >= 3),
        ("num_transits_ge_5", "koi_num_transits", lambda s: s >= 5),
        ("period_le_50d", "koi_period", lambda s: s <= 50),
        ("period_le_30d", "koi_period", lambda s: s <= 30),
        ("rp_lt_3p5", "koi_prad", lambda s: s < 3.5),
        ("has_rv_only", "has_rv", lambda s: s.astype(bool)),
        ("furlan_contam_le_1pct", "furlan_contam", lambda s: s.isna() | (s <= 0.01)),
    ]
    for name, col, fn in numeric_rules:
        if col in base.columns:
            rules.append((name, fn(pd.to_numeric(base[col], errors="coerce"))))
    rows = []
    for name, mask in rules:
        sub = base[mask.fillna(False)]
        e = sub["e_value"].to_numpy(float)
        rows.append(
            {
                "population": "thin_singles",
                "quality_rule": name,
                "n": int(len(e)),
                "hosts": int(sub["kepid"].nunique()),
                "mean_e": safe_stat(e, np.mean),
                "median_e": safe_stat(e, np.median),
                "q95_e": safe_percentile(e, 95),
                "frac_e_gt_0p10": frac_gt(e, 0.10),
            }
        )
    return pd.DataFrame(rows)


def top_outliers(df: pd.DataFrame, disk: str, system: str, n: int) -> pd.DataFrame:
    sub = finite_population(df, disk, system).sort_values("e_value", ascending=False).head(n).copy()
    cols = [
        "kepid",
        "kepoi_name",
        "kepler_name",
        "e_value",
        "koi_period",
        "koi_duration",
        "koi_impact",
        "koi_model_snr",
        "koi_num_transits",
        "koi_prad",
        "P_thick",
        "furlan_contam",
        "ruwe",
        "has_rv",
        "berger_feh",
        "berger_age",
    ]
    return sub[[c for c in cols if c in sub.columns]]


def label_overlap(sample: pd.DataFrame, ecc: pd.DataFrame) -> pd.DataFrame:
    if not {"old_disk", "old_system"}.intersection(ecc.columns):
        return pd.DataFrame([{"note": "No old disk/system labels in eccentricity source."}])
    merged = sample[["kepoi_name", "disk", "system"]].merge(
        ecc[[c for c in ["kepoi_name", "old_disk", "old_system"] if c in ecc.columns]],
        on="kepoi_name",
        how="inner",
    )
    if "old_disk" not in merged:
        return pd.DataFrame([{"note": "No old disk labels in eccentricity source."}])
    rows = []
    for (old_disk, new_disk), sub in merged.groupby(["old_disk", "disk"]):
        rows.append({"old_disk": old_disk, "new_disk": new_disk, "planets": len(sub)})
    return pd.DataFrame(rows).sort_values(["old_disk", "new_disk"])


def make_distribution_plot(df: pd.DataFrame, cfg: dict, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    colors = {
        "thin_singles": "#008b8b",
        "thick_singles": "#9b0000",
        "thin_multis": "#76b7b2",
        "thick_multis": "#c76f6f",
    }
    for disk, system, label in POPULATIONS:
        e = finite_population(df, disk, system)["e_value"].to_numpy(float)
        if len(e) == 0:
            continue
        axes[0].hist(e, bins=np.linspace(0, min(0.5, max(0.06, np.nanmax(e))), 35), histtype="step", density=True, color=colors[label], label=f"{label} (n={len(e)})")
        xs = np.sort(e)
        ys = np.arange(1, len(xs) + 1) / len(xs)
        axes[1].plot(xs, ys, color=colors[label], label=f"{label} (n={len(e)})")

    sagear_thin = cfg["sagear_targets"]["rayleigh_expected_e"]["thin_singles"]
    for ax in axes:
        ax.axvline(sagear_thin, color="#008b8b", linestyle="--", alpha=0.75, label="Sagear thin-single <e>")
        ax.set_xlabel("eccentricity estimate")
        ax.grid(alpha=0.2)
    axes[0].set_title("Eccentricity Density")
    axes[0].set_ylabel("density")
    axes[1].set_title("Eccentricity CDF")
    axes[1].set_ylabel("cumulative fraction")
    axes[1].set_xlim(0, 0.5)
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8, loc="lower right")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def finite_population(df: pd.DataFrame, disk: str, system: str) -> pd.DataFrame:
    sub = df[(df["disk"] == disk) & (df["system"] == system)].copy()
    sub = sub[np.isfinite(pd.to_numeric(sub["e_value"], errors="coerce"))]
    return sub


def safe_stat(values: np.ndarray, fn) -> float:
    return float(fn(values)) if len(values) else np.nan


def safe_percentile(values: np.ndarray, q: float) -> float:
    return float(np.nanpercentile(values, q)) if len(values) else np.nan


def frac_gt(values: np.ndarray, threshold: float) -> float:
    return float(np.mean(values > threshold)) if len(values) else np.nan


if __name__ == "__main__":
    main()
