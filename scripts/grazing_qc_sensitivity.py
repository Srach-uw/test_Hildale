"""Quantify the effect of a post-fit grazing-QC sensitivity on the hierarchy.

This is deliberately not the canonical replication fit.  It asks whether the
large eccentricity result is driven by systems whose ALDERAAN posterior places
more than five percent of its mass beyond the grazing boundary.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from hierarchical_rayleigh import fit_from_mass_matrix, load_population_masses


POPULATIONS = (
    ("thick", "single", "thick_singles"),
    ("thin", "single", "thin_singles"),
    ("thick", "multi", "thick_multis"),
    ("thin", "multi", "thin_multis"),
)


def fit_subset(summary: pd.DataFrame, sigma_count: int) -> list[dict[str, object]]:
    sigmas = np.linspace(1e-4, 1.0, sigma_count)
    rows: list[dict[str, object]] = []
    for disk, system, population in POPULATIONS:
        subset = summary[(summary["disk"] == disk) & (summary["system"] == system)]
        row: dict[str, object] = {"population": population, "n": len(subset)}
        if len(subset) < 5:
            row["status"] = "too_few_rows"
            rows.append(row)
            continue
        masses, e_grid = load_population_masses(
            subset,
            apply_transit_selection=True,
            selection_mode="manuscript_reciprocal",
        )
        fit = fit_from_mass_matrix(
            masses,
            e_grid,
            sigmas,
            apply_transit_selection=True,
            selection_mode="manuscript_reciprocal",
        )
        row.update({"status": "ok", **fit})
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary",
        default="outputs/eccentricity_posterior_summary_uniform_paired_full_qcprimary.csv",
    )
    parser.add_argument("--audit", default="outputs/postfit_qc_convention_audit.csv")
    parser.add_argument("--output", default="outputs/grazing_qc_sensitivity.csv")
    parser.add_argument("--sigma-count", type=int, default=300)
    args = parser.parse_args()

    summary = pd.read_csv(args.summary)
    audit = pd.read_csv(args.audit)
    required = {"kepoi_name", "weighted_grazing_exclude"}
    missing = sorted(required - set(audit.columns))
    if missing:
        raise ValueError(f"Grazing audit is missing required columns: {missing}")
    merged = summary.merge(
        audit[["kepoi_name", "weighted_grazing_exclude"]],
        on="kepoi_name",
        how="left",
        validate="one_to_one",
    )
    if merged["weighted_grazing_exclude"].isna().any():
        missing_count = int(merged["weighted_grazing_exclude"].isna().sum())
        raise ValueError(f"Missing grazing QC decision for {missing_count} posterior rows")
    merged["weighted_grazing_exclude"] = merged["weighted_grazing_exclude"].astype(bool)

    rows: list[dict[str, object]] = []
    for label, subset in (
        ("canonical", merged),
        ("weighted_grazing_filtered", merged[~merged["weighted_grazing_exclude"]]),
    ):
        for row in fit_subset(subset, args.sigma_count):
            row["sample"] = label
            rows.append(row)
    result = pd.DataFrame(rows)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    print(result[["sample", "population", "n", "status", "expected_e", "expected_e_lo", "expected_e_hi"]].to_string(index=False))
    print(f"weighted_grazing_excluded,{int(merged['weighted_grazing_exclude'].sum())}")
    print(f"wrote,{output.resolve()}")


if __name__ == "__main__":
    main()
