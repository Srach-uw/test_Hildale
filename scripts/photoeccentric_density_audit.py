from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits


REQUIRED_COLUMNS = {
    "kepoi_name",
    "disk",
    "system",
    "rho_circular_p16_solar",
    "rho_circular_median_solar",
    "rho_circular_p84_solar",
    "rho_catalog_solar",
    "e50",
}


def normalized_log_weights(log_weights: np.ndarray) -> np.ndarray:
    values = np.asarray(log_weights, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        raise ValueError("no finite nested-sampling weights")
    weights = np.zeros_like(values)
    shifted = values[finite] - np.max(values[finite])
    weights[finite] = np.exp(shifted)
    total = weights.sum()
    if not np.isfinite(total) or total <= 0:
        raise ValueError("nested-sampling weights do not normalize")
    return weights / total


def weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    quantiles: list[float],
) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights >= 0)
    values = values[valid]
    weights = weights[valid]
    if values.size == 0 or weights.sum() <= 0:
        return np.full(len(quantiles), np.nan)
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights) - 0.5 * weights
    cumulative /= weights.sum()
    return np.interp(quantiles, cumulative, values)


def prepare_planets(summary: pd.DataFrame) -> pd.DataFrame:
    missing = sorted(REQUIRED_COLUMNS - set(summary.columns))
    if missing:
        raise ValueError(f"missing required columns: {', '.join(missing)}")

    planets = summary.copy()
    numeric = [
        "rho_circular_p16_solar",
        "rho_circular_median_solar",
        "rho_circular_p84_solar",
        "rho_catalog_solar",
        "e50",
    ]
    for column in numeric:
        planets[column] = pd.to_numeric(planets[column], errors="coerce")

    valid = (
        np.isfinite(planets[numeric]).all(axis=1)
        & (planets["rho_circular_p16_solar"] > 0)
        & (planets["rho_circular_median_solar"] > 0)
        & (planets["rho_circular_p84_solar"] > 0)
        & (planets["rho_catalog_solar"] > 0)
        & (planets["rho_circular_p16_solar"] <= planets["rho_circular_median_solar"])
        & (planets["rho_circular_median_solar"] <= planets["rho_circular_p84_solar"])
    )
    planets = planets.loc[valid].copy()
    planets["population"] = (
        planets["disk"].astype(str).str.lower()
        + "_"
        + planets["system"].astype(str).str.lower()
    )
    planets["delta_log10_rho"] = np.log10(
        planets["rho_circular_median_solar"] / planets["rho_catalog_solar"]
    )
    planets["abs_delta_log10_rho"] = planets["delta_log10_rho"].abs()
    planets["rho_circular_width_dex"] = np.log10(
        planets["rho_circular_p84_solar"] / planets["rho_circular_p16_solar"]
    )
    return planets


def merge_qc_audit(
    summary: pd.DataFrame,
    qc_audit: pd.DataFrame,
) -> pd.DataFrame:
    if summary["kepoi_name"].duplicated().any():
        raise ValueError("summary has duplicate kepoi_name rows")
    if qc_audit["kepoi_name"].duplicated().any():
        raise ValueError("QC audit has duplicate kepoi_name rows")
    overlap = sorted(
        (set(summary.columns) & set(qc_audit.columns)) - {"kepoi_name"}
    )
    if overlap:
        raise ValueError(
            "QC audit duplicates summary columns: " + ", ".join(overlap)
        )
    merged = summary.merge(
        qc_audit,
        on="kepoi_name",
        how="left",
        validate="one_to_one",
    )
    if len(merged) != len(qc_audit) or merged["raw_result_file"].isna().any():
        raise ValueError("summary and QC audit do not contain the same planet IDs")
    return merged


def add_impact_widths(planets: pd.DataFrame) -> pd.DataFrame:
    required = {"raw_result_file", "alderaan_planet_index"}
    missing = sorted(required - set(planets.columns))
    if missing:
        raise ValueError(f"missing FITS provenance columns: {', '.join(missing)}")

    result = planets.copy()
    result["impact_p16"] = np.nan
    result["impact_p50"] = np.nan
    result["impact_p84"] = np.nan

    for raw_path, group in result.groupby("raw_result_file", sort=True):
        path = Path(str(raw_path))
        if not path.is_file():
            raise FileNotFoundError(f"raw ALDERAAN result not found: {path}")
        with fits.open(path, memmap=True) as hdul:
            samples = hdul["SAMPLES"].data
            names = set(samples.names)
            if "LN_WT" not in names:
                raise ValueError(f"{path} has no LN_WT column")
            weights = normalized_log_weights(np.asarray(samples["LN_WT"], dtype=float))
            for row_index, row in group.iterrows():
                column = f"IMPACT_{int(row['alderaan_planet_index'])}"
                if column not in names:
                    raise ValueError(f"{path} has no {column} column")
                q16, q50, q84 = weighted_quantile(
                    np.asarray(samples[column], dtype=float),
                    weights,
                    [0.16, 0.5, 0.84],
                )
                result.loc[row_index, ["impact_p16", "impact_p50", "impact_p84"]] = [
                    q16,
                    q50,
                    q84,
                ]
    result["impact_width"] = result["impact_p84"] - result["impact_p16"]
    return result


def summarize(planets: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for population, group in planets.groupby("population", sort=True):
        rows.append(
            {
                "population": population,
                "n": len(group),
                "median_delta_log10_rho": group["delta_log10_rho"].median(),
                "median_abs_delta_log10_rho": group["abs_delta_log10_rho"].median(),
                "median_rho_circular_width_dex": group[
                    "rho_circular_width_dex"
                ].median(),
                "median_planet_e50": group["e50"].median(),
                "mean_planet_e50": group["e50"].mean(),
            }
        )
    return pd.DataFrame(rows)


def summarize_impact_width(
    planets: pd.DataFrame,
    threshold: float = 0.4,
) -> pd.DataFrame:
    required = {"population", "impact_width", "abs_delta_log10_rho"}
    missing = sorted(required - set(planets.columns))
    if missing:
        raise ValueError(f"missing impact-audit columns: {', '.join(missing)}")
    rows = []
    for population, group in planets.groupby("population", sort=True):
        for label, subset in [
            ("constrained", group[group["impact_width"] <= threshold]),
            ("broad", group[group["impact_width"] > threshold]),
        ]:
            rows.append(
                {
                    "population": population,
                    "impact_width_class": label,
                    "impact_width_threshold": threshold,
                    "n": len(subset),
                    "median_impact_width": subset["impact_width"].median(),
                    "median_abs_delta_log10_rho": subset[
                        "abs_delta_log10_rho"
                    ].median(),
                    "median_rho_circular_width_dex": subset[
                        "rho_circular_width_dex"
                    ].median(),
                    "median_planet_e50": subset["e50"].median(),
                }
            )
    return pd.DataFrame(rows)


def fit_impact_constrained_populations(
    planets: pd.DataFrame,
    threshold: float = 0.4,
) -> pd.DataFrame:
    from hierarchical_rayleigh import (
        POPULATIONS,
        fit_from_mass_matrix,
        load_population_masses,
    )

    if "qc_primary_exclude" not in planets:
        raise ValueError(
            "impact-constrained hierarchy requires qc_primary_exclude"
        )
    primary_exclude = planets["qc_primary_exclude"].fillna(True).astype(bool)
    constrained = planets[
        (planets["impact_width"] <= threshold) & ~primary_exclude
    ].copy()
    sigmas = np.linspace(1e-4, 1.0, 2000)
    rows = []
    for disk, system, population in POPULATIONS:
        subset = constrained[
            (constrained["disk"] == disk) & (constrained["system"] == system)
        ].reset_index(drop=True)
        if len(subset) < 5:
            continue
        for selection_mode in [
            "legacy_forward_norm",
            "manuscript_reciprocal",
        ]:
            masses, e_grid = load_population_masses(
                subset,
                apply_transit_selection=True,
                selection_mode=selection_mode,
            )
            fit = fit_from_mass_matrix(
                masses,
                e_grid,
                sigmas,
                apply_transit_selection=True,
                selection_mode=selection_mode,
            )
            rows.append(
                {
                    "population": population,
                    "selection_mode": selection_mode,
                    "impact_width_threshold": threshold,
                    "n": len(subset),
                    "n_primary_qc_excluded_before_fit": int(
                        (
                            (planets["disk"] == disk)
                            & (planets["system"] == system)
                            & (planets["impact_width"] <= threshold)
                            & primary_exclude
                        ).sum()
                    ),
                    "expected_e": fit["expected_e"],
                    "expected_e_lo": fit["expected_e_lo"],
                    "expected_e_hi": fit["expected_e_hi"],
                    "boundary_flag": fit["boundary_flag"],
                }
            )
    return pd.DataFrame(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize the circular-density mismatch from a paired-impact "
            "ALDERAAN posterior audit."
        )
    )
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument(
        "--qc-audit",
        type=Path,
        help=(
            "Optional one-row-per-planet raw-FITS audit containing circular "
            "density and raw_result_file fields. It is joined by kepoi_name."
        ),
    )
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument(
        "--impact-output-csv",
        type=Path,
        help=(
            "Optional impact-width stratification. Requires raw_result_file and "
            "alderaan_planet_index in the summary and access to the raw FITS."
        ),
    )
    parser.add_argument(
        "--impact-fit-output-csv",
        type=Path,
        help=(
            "Optional Rayleigh fits for the impact-constrained subsets. "
            "Requires --impact-output-csv and accessible posterior_file paths."
        ),
    )
    parser.add_argument("--impact-width-threshold", type=float, default=0.4)
    return parser


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = build_parser().parse_args()
    source = pd.read_csv(args.summary)
    audit_source = None
    if args.qc_audit:
        audit_source = pd.read_csv(args.qc_audit)
        source = merge_qc_audit(source, audit_source)
    planets = prepare_planets(source)
    result = summarize(planets)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_csv, index=False)
    if args.impact_output_csv:
        if not np.isfinite(args.impact_width_threshold) or args.impact_width_threshold <= 0:
            raise ValueError("--impact-width-threshold must be positive and finite")
        planets = add_impact_widths(planets)
        impact = summarize_impact_width(planets, args.impact_width_threshold)
        args.impact_output_csv.parent.mkdir(parents=True, exist_ok=True)
        impact.to_csv(args.impact_output_csv, index=False)
        if args.impact_fit_output_csv:
            impact_fit = fit_impact_constrained_populations(
                planets,
                args.impact_width_threshold,
            )
            args.impact_fit_output_csv.parent.mkdir(parents=True, exist_ok=True)
            impact_fit.to_csv(args.impact_fit_output_csv, index=False)
    elif args.impact_fit_output_csv:
        raise ValueError("--impact-fit-output-csv requires --impact-output-csv")
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "input_file": args.summary.name,
            "input_sha256": sha256_file(args.summary),
            "input_rows": int(len(source)),
            "valid_rows": int(len(planets)),
            "populations": result.to_dict(orient="records"),
        }
        if args.qc_audit:
            payload["qc_audit_file"] = args.qc_audit.name
            payload["qc_audit_sha256"] = sha256_file(args.qc_audit)
        if args.impact_output_csv:
            payload["unique_raw_result_files"] = int(
                planets["raw_result_file"].nunique()
            )
            payload["impact_width_threshold"] = args.impact_width_threshold
            payload["impact_width_populations"] = impact.to_dict(orient="records")
            if args.impact_fit_output_csv:
                payload["impact_constrained_population_fits"] = impact_fit.to_dict(
                    orient="records"
                )
        args.output_json.write_text(
            json.dumps(payload, indent=2) + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
