from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import output_dir
from hierarchical_rayleigh import fit_rayleigh


def main() -> None:
    parser = argparse.ArgumentParser(description="Rayleigh population-fit sensitivity to ALDERAAN shape outliers.")
    parser.add_argument("--summary", default=None)
    parser.add_argument("--shape", default=None)
    parser.add_argument("--sigma-min", type=float, default=0.005)
    parser.add_argument("--sigma-max", type=float, default=0.25)
    parser.add_argument("--n-sigma", type=int, default=400)
    args = parser.parse_args()

    summary_path = Path(args.summary) if args.summary else output_dir() / "eccentricity_posterior_summary.csv"
    shape_path = Path(args.shape) if args.shape else output_dir() / "alderaan_shape_diagnostics.csv"
    summary = pd.read_csv(summary_path)
    shape = pd.read_csv(shape_path)
    df = summary.merge(
        shape[
            [
                "kepoi_name",
                "zeta_median",
                "alderaan_to_koi_duration_ratio",
                "alderaan_impact_med",
                "rho_frac_err",
                "koi_model_snr",
            ]
        ],
        on="kepoi_name",
        how="left",
        suffixes=("", "_shape"),
    )
    sigma_grid = np.linspace(args.sigma_min, args.sigma_max, args.n_sigma)

    filters = {
        "all": lambda x: np.ones(len(x), dtype=bool),
        "zeta_ge_0p7": lambda x: x["zeta_median_shape"].fillna(x["zeta_median"]) >= 0.7,
        "zeta_0p7_to_1p3": lambda x: (x["zeta_median_shape"].fillna(x["zeta_median"]) >= 0.7)
        & (x["zeta_median_shape"].fillna(x["zeta_median"]) <= 1.3),
        "no_duration_shift_25pct": lambda x: np.abs(np.log(x["alderaan_to_koi_duration_ratio"])) <= np.log(1.25),
        "alderaan_b_lt_0p85": lambda x: x["alderaan_impact_med"] < 0.85,
        "rho_frac_err_lt_0p35": lambda x: x["rho_frac_err"] < 0.35,
        "snr_ge_20": lambda x: x["koi_model_snr"] >= 20,
        "snr_ge_50": lambda x: x["koi_model_snr"] >= 50,
        "strict_clean": lambda x: (x["zeta_median_shape"].fillna(x["zeta_median"]) >= 0.7)
        & (x["zeta_median_shape"].fillna(x["zeta_median"]) <= 1.3)
        & (np.abs(np.log(x["alderaan_to_koi_duration_ratio"])) <= np.log(1.25))
        & (x["alderaan_impact_med"] < 0.85)
        & (x["rho_frac_err"] < 0.35),
    }

    rows = []
    for filter_name, fn in filters.items():
        filtered = df[fn(df)].copy()
        for disk, system, label in [
            ("thick", "single", "thick_singles"),
            ("thin", "single", "thin_singles"),
            ("thick", "multi", "thick_multis"),
            ("thin", "multi", "thin_multis"),
        ]:
            sub = filtered[(filtered["disk"] == disk) & (filtered["system"] == system)]
            row = {"filter": filter_name, "population": label, "n": len(sub)}
            if len(sub) >= 5:
                row.update(fit_rayleigh(sub["posterior_file"].tolist(), sigma_grid, apply_transit_selection=True))
            else:
                row["status"] = "skipped_n_lt_5"
            rows.append(row)

    out = pd.DataFrame(rows)
    path = output_dir() / "rayleigh_sensitivity_diagnostics.csv"
    out.to_csv(path, index=False)
    print(out[["filter", "population", "n", "expected_e", "expected_e_lo", "expected_e_hi"]].to_string(index=False))
    print(f"\nWrote: {path}")


if __name__ == "__main__":
    main()
