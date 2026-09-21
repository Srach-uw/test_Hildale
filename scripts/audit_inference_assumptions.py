"""Measure density-model and grid-resolution limits without fitting the data."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import trapezoid
from extract_eccentricity_posteriors_direct import density_log_likelihood
from hierarchical_rayleigh import rayleigh_grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root / "scripts/config.json").read_text())
    count = config["alderaan"]["eccentricity_grid_size"]
    e = np.linspace(0.0, 0.95, count)
    scales = np.array([0.0001, 0.001, 0.022 / np.sqrt(np.pi / 2), 0.1])
    rays, _ = rayleigh_grid(e, scales, False, "none")
    actual = trapezoid(rays * e[:, None], e, axis=0)
    analytic = scales * np.sqrt(np.pi / 2)
    near_peak = np.array([1.0 - 1e-9, 1.0 + 1e-9])
    jumps = {}
    for mode in ["split", "split-continuous"]:
        logp = density_log_likelihood(near_peak, 1.0, 0.2, 0.1, mode)
        jumps[mode] = float(np.exp(logp[1] - logp[0]))
    densities = pd.read_csv(root / "metadata/public_reconstruction_20260727/photoeccentric_density_audit.csv")
    result = {
        "scope": "Deterministic numerical-contract audit, not a new population fit",
        "density_peak_right_over_left": jumps,
        "production_grid_size": count,
        "production_grid_step": float(e[1] - e[0]),
        "rayleigh_grid": [
            {"sigma": float(s), "analytic_mean": float(a), "discrete_mean": float(d),
             "relative_mean_error": float(d / a - 1.0)}
            for s, a, d in zip(scales, analytic, actual)
        ],
        "archived_density_summary": densities.to_dict("records"),
        "interpretation": [
            "Median offsets and median widths cannot establish per-planet significance.",
            "The legacy split model has equal side mass and a discontinuous peak.",
            "The continuous two-piece normal has side mass proportional to side width.",
            "Neither asymmetric model is uniquely determined by reported 16/50/84 quantiles.",
            "The default symmetric-average branch is unchanged.",
            "Sub-grid Rayleigh means must not be interpreted as resolved measurements.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
