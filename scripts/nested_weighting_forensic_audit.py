"""Compare proper dynesty weighting with the equal-row forensic branch."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
PUBLISHED = {
    "rayleigh": {
        "thick_singles": 0.066,
        "thin_singles": 0.022,
        "thick_multis": 0.033,
        "thin_multis": 0.030,
    },
    "beta": {
        "thick_singles": 0.058,
        "thin_singles": 0.023,
        "thick_multis": 0.042,
        "thin_multis": 0.030,
    },
    "monotonic_beta": {
        "thick_singles": 0.041,
        "thin_singles": 0.022,
        "thick_multis": 0.025,
        "thin_multis": 0.026,
    },
    "half_gaussian": {
        "thick_singles": 0.066,
        "thin_singles": 0.025,
        "thick_multis": 0.037,
        "thin_multis": 0.033,
    },
}


def comparison_rows(
    branch: str, rayleigh: pd.DataFrame, shapes: pd.DataFrame
) -> list[dict]:
    rows: list[dict] = []
    for row in rayleigh.itertuples(index=False):
        published = PUBLISHED["rayleigh"][row.population]
        rows.append(
            {
                "branch": branch,
                "population": row.population,
                "model": "rayleigh",
                "n_planets": row.n,
                "local_value": row.expected_e,
                "local_lower": row.expected_e_lo,
                "local_upper": row.expected_e_hi,
                "published_value": published,
                "absolute_difference": abs(row.expected_e - published),
                "interval_contains_published": row.expected_e_lo <= published <= row.expected_e_hi,
            }
        )
    for row in shapes.itertuples(index=False):
        value = row.sigma if row.model == "half_gaussian" else row.mean_e
        published = PUBLISHED[row.model][row.population]
        rows.append(
            {
                "branch": branch,
                "population": row.population,
                "model": row.model,
                "n_planets": row.n,
                "local_value": value,
                "local_lower": np.nan,
                "local_upper": np.nan,
                "published_value": published,
                "absolute_difference": abs(value - published),
                "interval_contains_published": np.nan,
            }
        )
    return rows


def bic_rows(branch: str, rayleigh: pd.DataFrame, shapes: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    for row in rayleigh.itertuples(index=False):
        rows.append(
            {
                "branch": branch,
                "population": row.population,
                "model": "rayleigh",
                "bic": np.log(row.n) - 2.0 * row.ll_max,
            }
        )
    for row in shapes.itertuples(index=False):
        parameter_count = 1 if row.model == "half_gaussian" else 2
        rows.append(
            {
                "branch": branch,
                "population": row.population,
                "model": row.model,
                "bic": parameter_count * np.log(row.n) + 2.0 * row.nll,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "nested_weighting_forensic")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    inputs = {
        "dynesty_weighted": (
            ROOT / "outputs" / "rayleigh_population_fit_transit_selection_manuscript_reciprocal_EXACT_SIR1000_20260808.csv",
            ROOT / "outputs" / "table3_model_order_diagnostic_EXACT_SIR1000_20260808.csv",
        ),
        "equal_nested_rows_diagnostic": (
            ROOT / "outputs" / "rayleigh_population_fit_transit_selection_manuscript_reciprocal_EXACT_EQUAL_SIR1000_20260808.csv",
            ROOT / "outputs" / "table3_model_order_diagnostic_EXACT_EQUAL_SIR1000_20260808.csv",
        ),
    }
    comparison: list[dict] = []
    bics: list[dict] = []
    for branch, (rayleigh_path, shape_path) in inputs.items():
        if not rayleigh_path.exists() or not shape_path.exists():
            raise FileNotFoundError(f"Missing branch input for {branch}")
        rayleigh = pd.read_csv(rayleigh_path)
        shapes = pd.read_csv(shape_path)
        comparison.extend(comparison_rows(branch, rayleigh, shapes))
        bics.extend(bic_rows(branch, rayleigh, shapes))
    comparison_frame = pd.DataFrame(comparison)
    bic_frame = pd.DataFrame(bics)
    bic_frame["delta_bic"] = bic_frame.bic - bic_frame.groupby(
        ["branch", "population"]
    ).bic.transform("min")
    comparison_frame.to_csv(args.output_dir / "table3_comparison.csv", index=False)
    bic_frame.to_csv(args.output_dir / "local_bic_comparison.csv", index=False)
    rayleigh = comparison_frame[comparison_frame.model.eq("rayleigh")]
    summary = (
        rayleigh.groupby("branch")
        .agg(
            mean_absolute_error=("absolute_difference", "mean"),
            interval_overlaps=("interval_contains_published", "sum"),
        )
        .reset_index()
    )
    summary.to_csv(args.output_dir / "rayleigh_branch_summary.csv", index=False)
    lines = [
        "# Nested-weighting forensic audit",
        "",
        "ALDERAAN stores raw nested points and `LN_WT`; equal-row treatment is therefore "
        "a forensic implementation test, not a statistically valid posterior branch.",
        "",
        "## Rayleigh branch summary",
        "",
        summary.to_markdown(index=False),
        "",
        "## Table 3 comparison",
        "",
        comparison_frame.to_markdown(index=False),
        "",
        "## Local BIC values",
        "",
        bic_frame.sort_values(["branch", "population", "bic"]).to_markdown(index=False),
        "",
        "The local BIC calculation is a diagnostic MAP comparison. It does not reproduce "
        "Sagear's private NumPyro posterior or unpublished planet-level sample.",
    ]
    (args.output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print(summary.to_string(index=False))
    print(bic_frame.sort_values(["branch", "population", "bic"]).to_string(index=False))


if __name__ == "__main__":
    main()
