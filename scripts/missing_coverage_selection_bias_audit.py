from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


METRICS = [
    "koi_model_snr",
    "koi_period",
    "koi_prad",
    "koi_impact",
    "P_thick",
    "berger_logg",
    "berger_rad",
    "rho_log",
]


def build_target_level_frame(sample: pd.DataFrame, missing: pd.DataFrame) -> pd.DataFrame:
    required = {"koi_target", "disk", "system"} | set(METRICS)
    missing_columns = sorted(required - set(sample.columns))
    if missing_columns:
        raise ValueError(f"sample is missing required columns: {missing_columns}")
    missing_targets = set(missing["koi_target"].astype(str))
    frame = sample.copy()
    frame["koi_target"] = frame["koi_target"].astype(str)
    for metric in METRICS:
        frame[metric] = pd.to_numeric(frame[metric], errors="coerce")
    # A system is the independent ALDERAAN target. Summaries use finite system
    # values and preserve the population label from the exact inventory.
    aggregations = {
        metric: (metric, "median") for metric in METRICS
    }
    aggregations["max_snr"] = ("koi_model_snr", "max")
    aggregations["min_impact"] = ("koi_impact", "min")
    target = (
        frame.groupby("koi_target", as_index=False)
        .agg(
            kepid=("kepid", "first"),
            disk=("disk", "first"),
            system=("system", "first"),
            **aggregations,
        )
    )
    target["coverage_state"] = np.where(
        target["koi_target"].isin(missing_targets), "missing", "available"
    )
    target["missing_target_reason"] = np.where(
        target["coverage_state"].eq("missing"),
        "no_usable_result",
        "available",
    )
    return target


def compare_target_groups(targets: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (disk, system), group in targets.groupby(["disk", "system"]):
        missing = group[group["coverage_state"].eq("missing")]
        available = group[group["coverage_state"].eq("available")]
        for metric in METRICS + ["max_snr", "min_impact"]:
            m = pd.to_numeric(missing[metric], errors="coerce").dropna()
            a = pd.to_numeric(available[metric], errors="coerce").dropna()
            row: dict[str, object] = {
                "disk": disk,
                "system": system,
                "metric": metric,
                "missing_targets": len(missing),
                "available_targets": len(available),
                "missing_median": m.median() if len(m) else np.nan,
                "available_median": a.median() if len(a) else np.nan,
                "missing_mean": m.mean() if len(m) else np.nan,
                "available_mean": a.mean() if len(a) else np.nan,
            }
            if len(m) and len(a):
                row["median_difference_missing_minus_available"] = (
                    float(m.median() - a.median())
                )
                try:
                    from scipy.stats import ks_2samp

                    result = ks_2samp(m, a, alternative="two-sided", mode="auto")
                    row["ks_statistic"] = float(result.statistic)
                    row["ks_pvalue"] = float(result.pvalue)
                except ImportError:
                    row["ks_statistic"] = np.nan
                    row["ks_pvalue"] = np.nan
            else:
                row["median_difference_missing_minus_available"] = np.nan
                row["ks_statistic"] = np.nan
                row["ks_pvalue"] = np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def write_markdown(targets: pd.DataFrame, comparison: pd.DataFrame, path: Path) -> None:
    lines = [
        "# Missing-posterior selection-bias audit",
        "",
        "This audit compares published-inventory systems with and without usable",
        "ALDERAAN results using only pre-fit catalog metadata. It does not infer",
        "eccentricity for missing systems.",
        "",
        "| disk | system | missing targets | available targets | median S/N missing | median S/N available | median impact missing | median impact available |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for (disk, system), group in targets.groupby(["disk", "system"]):
        missing = group[group["coverage_state"].eq("missing")]
        available = group[group["coverage_state"].eq("available")]
        lines.append(
            f"| {disk} | {system} | {len(missing)} | {len(available)} | "
            f"{missing['max_snr'].median():.2f} | {available['max_snr'].median():.2f} | "
            f"{missing['min_impact'].median():.3f} | {available['min_impact'].median():.3f} |"
        )
    lines.extend(
        [
            "",
            "A strong difference means the current posterior sample is not a",
            "random subset of the published inventory. The completion run is then",
            "needed before interpreting population-level differences.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit missing-target selection bias.")
    parser.add_argument("--sample", default="outputs/sagear2026_planet_inventory_pre_visual_qc.csv")
    parser.add_argument("--missing", default="outputs/sagear2026_planets_missing_alderaan_target_results.csv")
    parser.add_argument("--out", default="outputs/missing_coverage_selection_bias.csv")
    args = parser.parse_args()

    sample = pd.read_csv(args.sample)
    missing = pd.read_csv(args.missing)
    targets = build_target_level_frame(sample, missing)
    comparison = compare_target_groups(targets)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    targets.to_csv(out.with_name(f"{out.stem}_target_level.csv"), index=False)
    comparison.to_csv(out, index=False)
    write_markdown(targets, comparison, out.with_suffix(".md"))
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
