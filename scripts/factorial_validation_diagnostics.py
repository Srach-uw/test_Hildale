"""Secondary diagnostics for paired ALDERAAN factorial-validation outputs.

This module operates on the compact CSVs produced by
``compare_factorial_validation.py``. It does not read light curves, rerun
ALDERAAN, or regenerate eccentricity posteriors.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


REQUIRED = {
    "comparison_id",
    "effect",
    "koi_target",
    "kepoi_name",
    "disk_baseline",
    "system_baseline",
    "e16_baseline",
    "e50_baseline",
    "e84_baseline",
    "e16_comparison",
    "e50_comparison",
    "e84_comparison",
    "impact_p50_baseline",
    "impact_p50_comparison",
    "delta_e",
    "delta_impact",
    "delta_t14_hr",
    "fractional_delta_rp_over_rs",
    "qc_primary_exclude_baseline",
    "qc_primary_exclude_comparison",
}


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def clean_pairs(frame: pd.DataFrame) -> pd.DataFrame:
    require_columns(frame, REQUIRED, "paired-planets table")
    result = frame.copy()
    result["qc_exclude"] = as_bool(result["qc_primary_exclude_baseline"]) | as_bool(
        result["qc_primary_exclude_comparison"]
    )
    numeric = [
        "e16_baseline",
        "e50_baseline",
        "e84_baseline",
        "e16_comparison",
        "e50_comparison",
        "e84_comparison",
        "impact_p50_baseline",
        "impact_p50_comparison",
        "delta_e",
        "delta_impact",
        "delta_t14_hr",
        "fractional_delta_rp_over_rs",
    ]
    for column in numeric:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    result["e_width_baseline"] = result["e84_baseline"] - result["e16_baseline"]
    result["e_width_comparison"] = result["e84_comparison"] - result["e16_comparison"]
    result["abs_delta_e"] = result["delta_e"].abs()
    result["population"] = (
        result["disk_baseline"].astype(str).str.lower()
        + "_"
        + result["system_baseline"].astype(str).str.lower()
    )
    finite = np.isfinite(result["delta_e"])
    return result[finite & ~result["qc_exclude"]].reset_index(drop=True)


def cluster_bootstrap(
    frame: pd.DataFrame,
    column: str,
    *,
    n_bootstrap: int,
    seed: int,
) -> tuple[float, float]:
    systems = frame["koi_target"].drop_duplicates().to_numpy()
    if len(systems) == 0:
        return np.nan, np.nan
    grouped_values = [
        frame.loc[frame["koi_target"].eq(system), column].to_numpy(float)
        for system in systems
    ]
    rng = np.random.default_rng(seed)
    values = np.empty(n_bootstrap, dtype=float)
    for index in range(n_bootstrap):
        sampled = rng.integers(0, len(systems), size=len(systems))
        values[index] = np.nanmedian(np.concatenate([grouped_values[item] for item in sampled]))
    return tuple(np.nanquantile(values, [0.025, 0.975]))


def population_summaries(
    pairs: pd.DataFrame,
    *,
    n_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_columns = ["comparison_id", "effect", "population"]
    for group_index, (keys, group) in enumerate(pairs.groupby(group_columns, sort=True)):
        comparison_id, effect, population = keys
        low, high = cluster_bootstrap(
            group,
            "delta_e",
            n_bootstrap=n_bootstrap,
            seed=seed + group_index,
        )
        rows.append(
            {
                "comparison_id": comparison_id,
                "effect": effect,
                "population": population,
                "n_planets": len(group),
                "n_systems": group["koi_target"].nunique(),
                "median_delta_e": group["delta_e"].median(),
                "median_delta_e_ci_low": low,
                "median_delta_e_ci_high": high,
                "mean_delta_e": group["delta_e"].mean(),
                "median_abs_delta_e": group["abs_delta_e"].median(),
                "p95_abs_delta_e": group["abs_delta_e"].quantile(0.95),
                "fraction_positive": group["delta_e"].gt(0).mean(),
                "median_baseline_e_width": group["e_width_baseline"].median(),
                "median_abs_delta_over_e_width": np.nanmedian(
                    group["abs_delta_e"] / group["e_width_baseline"]
                ),
                "interpretation": (
                    "underpowered_fewer_than_5_systems"
                    if group["koi_target"].nunique() < 5
                    else "descriptive_targeted_sample"
                ),
            }
        )
    return pd.DataFrame(rows)


def leverage_table(pairs: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for comparison_id, group in pairs.groupby("comparison_id", sort=True):
        full_median = group["delta_e"].median()
        system_scores = (
            group.groupby("koi_target", as_index=False)
            .agg(
                n_planets=("kepoi_name", "size"),
                median_delta_e=("delta_e", "median"),
                max_abs_delta_e=("abs_delta_e", "max"),
            )
            .sort_values("max_abs_delta_e", ascending=False)
        )
        for rank, row in enumerate(system_scores.itertuples(index=False), start=1):
            retained = group[~group["koi_target"].eq(row.koi_target)]
            leave_out = retained["delta_e"].median() if not retained.empty else np.nan
            rows.append(
                {
                    "comparison_id": comparison_id,
                    "koi_target": row.koi_target,
                    "rank_by_max_abs_delta_e": rank,
                    "n_planets": row.n_planets,
                    "system_median_delta_e": row.median_delta_e,
                    "system_max_abs_delta_e": row.max_abs_delta_e,
                    "full_median_delta_e": full_median,
                    "leave_system_out_median_delta_e": leave_out,
                    "leave_system_out_shift": leave_out - full_median,
                }
            )
    return pd.DataFrame(rows)


def correlation_table(pairs: pd.DataFrame, selection: pd.DataFrame | None) -> pd.DataFrame:
    joined = pairs.copy()
    if selection is not None and "kepoi_name" in selection:
        optional = [
            column
            for column in ("kepoi_name", "ld_abs_max_delta", "catalog_impact", "selection_reason")
            if column in selection
        ]
        joined = joined.merge(
            selection[optional].drop_duplicates("kepoi_name"),
            on="kepoi_name",
            how="left",
            validate="many_to_one",
        )
    features = [
        "e50_baseline",
        "e_width_baseline",
        "impact_p50_baseline",
        "delta_impact",
        "delta_t14_hr",
        "fractional_delta_rp_over_rs",
        "ld_abs_max_delta",
        "catalog_impact",
    ]
    rows: list[dict[str, object]] = []
    for comparison_id, group in joined.groupby("comparison_id", sort=True):
        for response in ("delta_e", "abs_delta_e"):
            for feature in features:
                if feature not in group:
                    continue
                x = pd.to_numeric(group[feature], errors="coerce")
                y = pd.to_numeric(group[response], errors="coerce")
                valid = np.isfinite(x) & np.isfinite(y)
                n = int(valid.sum())
                if n < 5 or x[valid].nunique() < 2 or y[valid].nunique() < 2:
                    rho = p_value = np.nan
                else:
                    rho, p_value = spearmanr(x[valid], y[valid])
                rows.append(
                    {
                        "comparison_id": comparison_id,
                        "response": response,
                        "feature": feature,
                        "n_planets": n,
                        "spearman_rho": rho,
                        "nominal_p_value": p_value,
                        "interpretation": "exploratory_multiple_tests_not_corrected",
                    }
                )
    return pd.DataFrame(rows)


def robustness_table(pairs: pd.DataFrame, *, seed: int, trials: int = 1000) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    rng = np.random.default_rng(seed)
    for comparison_id, group in pairs.groupby("comparison_id", sort=True):
        full = group["delta_e"].median()
        systems = group["koi_target"].drop_duplicates().to_numpy()
        ranked = (
            group.groupby("koi_target")["abs_delta_e"].max().sort_values(ascending=False).index
        )
        variants = {
            "all_qc_pass": group,
            "exclude_highest_leverage_system": group[~group["koi_target"].eq(ranked[0])],
            "impact_below_0p8_both_arms": group[
                group["impact_p50_baseline"].lt(0.8)
                & group["impact_p50_comparison"].lt(0.8)
            ],
            "e84_below_0p9_both_arms": group[
                group["e84_baseline"].lt(0.9) & group["e84_comparison"].lt(0.9)
            ],
        }
        for label, subset in variants.items():
            rows.append(
                {
                    "comparison_id": comparison_id,
                    "analysis": label,
                    "n_planets": len(subset),
                    "n_systems": subset["koi_target"].nunique(),
                    "median_delta_e": subset["delta_e"].median() if not subset.empty else np.nan,
                    "shift_from_full_median": (
                        subset["delta_e"].median() - full if not subset.empty else np.nan
                    ),
                }
            )
        if len(systems) >= 5:
            n_drop = max(1, int(np.ceil(0.1 * len(systems))))
            shifts = np.empty(trials)
            for index in range(trials):
                dropped = set(rng.choice(systems, size=n_drop, replace=False))
                retained = group[~group["koi_target"].isin(dropped)]
                shifts[index] = retained["delta_e"].median() - full
            rows.append(
                {
                    "comparison_id": comparison_id,
                    "analysis": "leave_10pct_systems_out_trials",
                    "n_planets": len(group),
                    "n_systems": len(systems),
                    "median_delta_e": full,
                    "shift_from_full_median": np.nanmedian(np.abs(shifts)),
                    "shift_p95": np.nanquantile(np.abs(shifts), 0.95),
                    "shift_max": np.nanmax(np.abs(shifts)),
                    "n_trials": trials,
                }
            )
    return pd.DataFrame(rows)


def make_plot(pairs: pd.DataFrame, output: Path) -> None:
    comparisons = list(pairs["comparison_id"].drop_duplicates())
    comparison = comparisons[0]
    group = pairs[pairs["comparison_id"].eq(comparison)].copy()
    colors = {
        "thin_single": "#008b8b",
        "thin_multi": "#4c78a8",
        "thick_single": "#a51616",
        "thick_multi": "#e07b39",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    axes[0, 0].hist(group["delta_e"], bins=min(12, max(5, len(group) // 2)), color="#666666")
    axes[0, 0].axvline(0, color="black", linewidth=1)
    axes[0, 0].set(xlabel="paired delta e", ylabel="planets", title=comparison)
    for population, subset in group.groupby("population"):
        color = colors.get(population, "#666666")
        axes[0, 1].scatter(
            subset["e50_baseline"], subset["e50_comparison"], label=population, color=color
        )
        axes[1, 0].scatter(
            subset["impact_p50_baseline"], subset["delta_e"], label=population, color=color
        )
    limit = max(group["e50_baseline"].max(), group["e50_comparison"].max()) * 1.05
    axes[0, 1].plot([0, limit], [0, limit], color="black", linewidth=1)
    axes[0, 1].set(xlabel="baseline e50", ylabel="comparison e50", xlim=(0, limit), ylim=(0, limit))
    axes[1, 0].axhline(0, color="black", linewidth=1)
    axes[1, 0].set(xlabel="baseline impact median", ylabel="paired delta e")
    order = [name for name in colors if name in set(group["population"])]
    axes[1, 1].boxplot(
        [group.loc[group["population"].eq(name), "delta_e"] for name in order],
        tick_labels=[name.replace("_", "\n") for name in order],
    )
    axes[1, 1].axhline(0, color="black", linewidth=1)
    axes[1, 1].set(ylabel="paired delta e", title="Targeted population subsets")
    axes[0, 1].legend(fontsize=8)
    fig.suptitle("ALDERAAN factorial validation diagnostics")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_report(
    path: Path,
    pairs: pd.DataFrame,
    populations: pd.DataFrame,
    leverage: pd.DataFrame,
    correlations: pd.DataFrame,
    robustness: pd.DataFrame,
    discovery: pd.DataFrame | None,
) -> None:
    lines = [
        "# Factorial validation secondary diagnostics",
        "",
        "These are paired, targeted validation diagnostics. They are not population estimates or",
        "formal discovery tests. Nominal correlation p-values are uncorrected and are reported only",
        "to identify engineering relationships worth checking in the completed matrix.",
        "",
        "## Coverage",
        "",
        f"- QC-passing paired planets: {len(pairs)}",
        f"- Paired target systems: {pairs['koi_target'].nunique()}",
        f"- Available comparisons: {pairs['comparison_id'].nunique()}",
    ]
    if discovery is not None:
        present = int(discovery["status"].eq("present").sum())
        lines += [f"- Discovered FITS: {present}/{len(discovery)}"]
    lines += [
        "",
        "## Population summaries",
        "",
        populations.to_markdown(index=False),
        "",
        "Any row with fewer than five independent systems is explicitly underpowered.",
        "",
        "## Highest-leverage systems",
        "",
        leverage.sort_values(["comparison_id", "rank_by_max_abs_delta_e"]).groupby(
            "comparison_id", group_keys=False
        ).head(5).to_markdown(index=False),
        "",
        "## Robustness",
        "",
        robustness.to_markdown(index=False),
        "",
        "## Largest exploratory correlations",
        "",
        correlations.assign(abs_rho=correlations["spearman_rho"].abs())
        .sort_values("abs_rho", ascending=False)
        .head(12)
        .drop(columns="abs_rho")
        .to_markdown(index=False),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def run(
    paired_path: Path,
    output_dir: Path,
    *,
    selection_path: Path | None = None,
    discovery_path: Path | None = None,
    n_bootstrap: int = 10_000,
    seed: int = 20260713,
) -> dict[str, Path]:
    pairs = clean_pairs(pd.read_csv(paired_path))
    if pairs.empty:
        raise ValueError("No QC-passing paired planets are available")
    selection = pd.read_csv(selection_path) if selection_path and selection_path.is_file() else None
    discovery = pd.read_csv(discovery_path) if discovery_path and discovery_path.is_file() else None
    populations = population_summaries(pairs, n_bootstrap=n_bootstrap, seed=seed)
    leverage = leverage_table(pairs)
    correlations = correlation_table(pairs, selection)
    robustness = robustness_table(pairs, seed=seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "populations": output_dir / "factorial_validation_population_diagnostics.csv",
        "leverage": output_dir / "factorial_validation_system_leverage.csv",
        "correlations": output_dir / "factorial_validation_correlations.csv",
        "robustness": output_dir / "factorial_validation_robustness.csv",
        "plot": output_dir / "factorial_validation_secondary_diagnostics.png",
        "report": output_dir / "factorial_validation_secondary_diagnostics.md",
    }
    populations.to_csv(paths["populations"], index=False)
    leverage.to_csv(paths["leverage"], index=False)
    correlations.to_csv(paths["correlations"], index=False)
    robustness.to_csv(paths["robustness"], index=False)
    make_plot(pairs, paths["plot"])
    write_report(paths["report"], pairs, populations, leverage, correlations, robustness, discovery)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--paired",
        required=True,
        help="Detailed factorial_validation_paired_planets.csv table.",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--selection", default=None)
    parser.add_argument("--discovery", default=None)
    parser.add_argument("--bootstrap-replicates", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260713)
    args = parser.parse_args()
    if args.bootstrap_replicates <= 0:
        parser.error("--bootstrap-replicates must be positive")
    paths = run(
        Path(args.paired),
        Path(args.output_dir),
        selection_path=Path(args.selection) if args.selection else None,
        discovery_path=Path(args.discovery) if args.discovery else None,
        n_bootstrap=args.bootstrap_replicates,
        seed=args.seed,
    )
    for label, path in paths.items():
        print(f"{label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
