from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from common import load_config, read_berger_table2
from extract_eccentricity_posteriors_direct import _add_density, process_target


@dataclass(frozen=True)
class Variant:
    name: str
    nested_weight_mode: str = "dynesty"
    density_error_mode: str = "symmetric-average"
    density_error_scale: float = 1.0
    density_offset_dex: float = 0.0
    e_max: float = 0.95


VARIANTS = (
    Variant("baseline_50k"),
    Variant("equal_nested_weights", nested_weight_mode="equal"),
    Variant("density_error_half", density_error_scale=0.5),
    Variant("density_error_double", density_error_scale=2.0),
    Variant("density_offset_minus0p1dex", density_offset_dex=-0.1),
    Variant("density_offset_plus0p1dex", density_offset_dex=0.1),
    Variant("emax_0p92", e_max=0.92),
    Variant("split_density_errors", density_error_mode="split"),
)


def summarize_deltas(paired: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (variant, population), group in paired.groupby(["variant", "population"], sort=False):
        delta = pd.to_numeric(group["delta_e50"], errors="coerce").dropna().to_numpy(float)
        if len(delta) == 0:
            continue
        rows.append(
            {
                "variant": variant,
                "population": population,
                "n_planets": len(delta),
                "median_delta_e50": float(np.median(delta)),
                "median_abs_delta_e50": float(np.median(np.abs(delta))),
                "p95_abs_delta_e50": float(np.quantile(np.abs(delta), 0.95)),
                "max_abs_delta_e50": float(np.max(np.abs(delta))),
            }
        )
    return pd.DataFrame(rows)


def result_path(
    target: str,
    source: str,
    archive_results: Path,
    cloud_results: Path,
) -> Path:
    if source == "original_alderaan_archive":
        return archive_results / f"{target}-results.fits"
    if source == "cloud_missing_run":
        return cloud_results / target / f"{target}-results.fits"
    raise ValueError(f"Unknown transit_fit_source for {target}: {source}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Re-extract the highest-leverage systems under explicit posterior sensitivity conventions."
    )
    parser.add_argument("--config", default="config.json")
    parser.add_argument("--sample", default="outputs/canonical_sample_old_astropy_rawcc.csv")
    parser.add_argument("--leverage", default="outputs/uniform_paired_top100_leverage_catalog_audit.csv")
    parser.add_argument("--canonical-summary", default="outputs/eccentricity_posterior_summary_uniform_paired_full.csv")
    parser.add_argument(
        "--archive-results",
        default="inputs/alderaan_archive_raw_20260612/ALDERAAN_posteriors",
    )
    parser.add_argument(
        "--cloud-results",
        default="alderaan_project/Results/sagear_missing",
    )
    parser.add_argument("--output-dir", default="outputs/targeted_posterior_sensitivity")
    parser.add_argument("--top-targets", type=int, default=100)
    parser.add_argument("--n-proposals", type=int, default=50_000)
    parser.add_argument("--min-importance-ess", type=float, default=100.0)
    args = parser.parse_args()

    if args.top_targets <= 0 or args.n_proposals <= 0:
        parser.error("--top-targets and --n-proposals must be positive")

    cfg = load_config(args.config)
    leverage = pd.read_csv(args.leverage).head(args.top_targets).copy()
    required = {"koi_target", "transit_fit_source", "population"}
    missing = sorted(required - set(leverage))
    if missing:
        raise ValueError(f"Leverage table missing required columns: {missing}")
    source_counts = leverage.groupby("koi_target")["transit_fit_source"].nunique()
    if (source_counts > 1).any():
        bad = source_counts[source_counts > 1].index.tolist()
        raise ValueError(f"Selected targets map to multiple FITS sources: {bad}")

    target_source = leverage.drop_duplicates("koi_target").set_index("koi_target")["transit_fit_source"]
    selected_targets = target_source.index.astype(str).tolist()
    sample = _add_density(pd.read_csv(args.sample), read_berger_table2(cfg))
    sample = sample[sample["koi_target"].astype(str).isin(selected_targets)].copy()
    missing_targets = sorted(set(selected_targets) - set(sample["koi_target"].astype(str)))
    if missing_targets:
        raise ValueError(f"Selected leverage targets missing from canonical sample: {missing_targets}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    archive_results = Path(args.archive_results)
    cloud_results = Path(args.cloud_results)
    grid_size = int(cfg["alderaan"]["eccentricity_grid_size"])
    omega_grid = np.linspace(
        0.0,
        2.0 * np.pi,
        int(cfg["alderaan"]["omega_grid_size"]),
        endpoint=False,
    )

    all_summaries: list[pd.DataFrame] = []
    all_exclusions: list[pd.DataFrame] = []
    for variant in VARIANTS:
        print(f"Running {variant.name} on {len(selected_targets)} systems")
        e_grid = np.linspace(0.0, variant.e_max, grid_size)
        posterior_dir = output_dir / variant.name / "posteriors"
        posterior_dir.mkdir(parents=True, exist_ok=True)
        summaries: list[dict[str, object]] = []
        exclusions: list[dict[str, object]] = []
        for target in selected_targets:
            planets = (
                sample[sample["koi_target"].astype(str).eq(target)]
                .sort_values("koi_period")
                .reset_index(drop=True)
            )
            path = result_path(
                target,
                str(target_source.loc[target]),
                archive_results,
                cloud_results,
            )
            if not path.exists():
                raise FileNotFoundError(f"Selected result FITS is missing: {path}")
            target_summary, target_exclusions = process_target(
                path,
                planets,
                posterior_dir,
                e_grid,
                omega_grid,
                n_proposals=args.n_proposals,
                e_max=variant.e_max,
                density_error_mode=variant.density_error_mode,
                period_tol=0.01,
                min_importance_ess=args.min_importance_ess,
                allow_density_error_fallback=False,
                nested_weight_mode=variant.nested_weight_mode,
                density_error_scale=variant.density_error_scale,
                density_offset_dex=variant.density_offset_dex,
            )
            summaries.extend(target_summary)
            exclusions.extend(target_exclusions)

        variant_summary = pd.DataFrame(summaries)
        variant_summary["variant"] = variant.name
        variant_summary["transit_fit_source"] = variant_summary["koi_target"].map(target_source)
        variant_exclusions = pd.DataFrame(exclusions)
        if not variant_exclusions.empty:
            variant_exclusions["variant"] = variant.name
        variant_summary.to_csv(output_dir / variant.name / "summary.csv", index=False)
        variant_exclusions.to_csv(output_dir / variant.name / "exclusions.csv", index=False)
        all_summaries.append(variant_summary)
        all_exclusions.append(variant_exclusions)

    combined = pd.concat(all_summaries, ignore_index=True)
    baseline = combined[combined["variant"].eq("baseline_50k")][
        ["kepoi_name", "e50", "disk", "system"]
    ].rename(columns={"e50": "e50_baseline_50k"})
    paired = combined.merge(baseline, on=["kepoi_name", "disk", "system"], how="inner")
    paired["population"] = paired["disk"].astype(str) + "_" + paired["system"].astype(str)
    paired["delta_e50"] = paired["e50"] - paired["e50_baseline_50k"]
    paired.to_csv(output_dir / "paired_variant_deltas.csv", index=False)
    summarize_deltas(paired).to_csv(output_dir / "variant_delta_summary.csv", index=False)

    canonical = pd.read_csv(args.canonical_summary)[["kepoi_name", "e50"]].rename(
        columns={"e50": "e50_canonical_150k"}
    )
    calibration = baseline.merge(canonical, on="kepoi_name", how="left")
    calibration["delta_50k_minus_150k"] = (
        calibration["e50_baseline_50k"] - calibration["e50_canonical_150k"]
    )
    calibration.to_csv(output_dir / "monte_carlo_calibration.csv", index=False)

    if all_exclusions:
        pd.concat(all_exclusions, ignore_index=True).to_csv(
            output_dir / "all_exclusions.csv",
            index=False,
        )
    print(f"Wrote targeted sensitivity products under {output_dir.resolve()}")


if __name__ == "__main__":
    main()
