from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize the circular-density mismatch from a paired-impact "
            "ALDERAAN posterior audit."
        )
    )
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    planets = prepare_planets(pd.read_csv(args.summary))
    result = summarize(planets)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_csv, index=False)
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "input_rows": int(len(pd.read_csv(args.summary, usecols=["kepoi_name"]))),
            "valid_rows": int(len(planets)),
            "populations": result.to_dict(orient="records"),
        }
        args.output_json.write_text(
            json.dumps(payload, indent=2) + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
