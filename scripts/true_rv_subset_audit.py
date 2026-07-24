"""Build the Sagear true-radial-velocity posterior subset.

The published paper explicitly repeats the eccentricity analysis using only
hosts with measured (rather than inferred) radial velocities.  This script
joins that flag to the extracted posterior summary and writes a strict,
reproducible hierarchy input for the existing Rayleigh fitter.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def build(summary_path: Path, inventory_path: Path, output_path: Path) -> pd.DataFrame:
    summary = pd.read_csv(summary_path)
    inventory = pd.read_csv(inventory_path)
    required = {"kepid", "has_measured_velocity"}
    missing = required - set(inventory.columns)
    if missing:
        raise ValueError(f"Inventory is missing required columns: {sorted(missing)}")
    flag = inventory[["kepid", "has_measured_velocity"]].drop_duplicates("kepid")
    merged = summary.merge(flag, on="kepid", how="left", validate="many_to_one")
    if merged["has_measured_velocity"].isna().any():
        raise ValueError("Some posterior rows did not receive a measured-velocity flag")
    result = merged[merged["has_measured_velocity"].astype(bool)].copy()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", default="outputs/eccentricity_posterior_summary_dynesty_published_inventory_pre_visual_qc.csv")
    parser.add_argument("--inventory", default="outputs/sagear2026_planet_inventory_pre_visual_qc.csv")
    parser.add_argument("--output", default="outputs/eccentricity_posterior_summary_true_rv_dynesty.csv")
    args = parser.parse_args()
    result = build(Path(args.summary), Path(args.inventory), Path(args.output))
    print(f"true_rv_posterior_rows,{len(result)}")
    print(result.groupby(["disk", "system"]).size().to_string())


if __name__ == "__main__":
    main()
