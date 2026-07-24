"""Audit the posterior-export contract in archived Sagear and ALDERAAN files.

This is deliberately descriptive. It does not promote any weighting convention
to the canonical population result.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits


def legacy_record(path: Path, target: str) -> dict[str, object]:
    equal_path = path / "jomnest_post_equal_weights.dat"
    pkl_path = path / "posteriors.pkl"
    posterior = pickle.load(pkl_path.open("rb"))
    samples = posterior["posterior_samples"]["unnamed"]
    equal = np.loadtxt(equal_path)
    if equal.ndim == 1:
        equal = equal[None, :]
    parameter_block_matches = equal[:, : samples.shape[1]]
    same_parameter_rows = bool(
        parameter_block_matches.shape == samples.shape
        and np.allclose(parameter_block_matches, samples, rtol=1e-10, atol=1e-12)
    )
    return {
        "source": "Sagear legacy photoeccentric",
        "target": target,
        "equal_file": equal_path.name,
        "equal_rows": int(equal.shape[0]),
        "equal_columns": int(equal.shape[1]),
        "pkl_rows": int(samples.shape[0]),
        "pkl_columns": int(samples.shape[1]),
        "equal_parameter_block_matches_pkl": same_parameter_rows,
        "has_nested_weight_column": False,
        "weight_contract": "post_equal_weights export; rows are intended equal-weight draws",
    }


def current_record(path: Path) -> dict[str, object]:
    with fits.open(path, memmap=False) as hdul:
        table = hdul["SAMPLES"].data
        names = set(table.names or [])
        weight = np.asarray(table["LN_WT"], dtype=float)
        finite = np.isfinite(weight)
        stable = weight[finite] - np.nanmax(weight[finite])
        weights = np.exp(np.clip(stable, -745.0, 0.0))
        weights /= weights.sum()
        ess = float(1.0 / np.sum(weights**2))
        return {
            "source": "Current ALDERAAN FITS",
            "target": path.stem.replace("-results", ""),
            "equal_file": "",
            "equal_rows": "",
            "equal_columns": "",
            "pkl_rows": "",
            "pkl_columns": "",
            "equal_matches_pkl_shape": "",
            "has_nested_weight_column": "LN_WT" in names,
            "weight_contract": "raw nested rows with LN_WT; equal raw rows are not posterior draws",
            "raw_rows": int(len(table)),
            "finite_weight_rows": int(finite.sum()),
            "nested_weight_ess": ess,
            "weight_min": float(np.nanmin(weight[finite])),
            "weight_max": float(np.nanmax(weight[finite])),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="outputs/legacy_posterior_export_audit.csv")
    parser.add_argument(
        "--fits",
        default="alderaan_project/Results/sagear_missing",
        help="Directory containing current ALDERAAN result FITS files.",
    )
    args = parser.parse_args()

    root = Path("external/ssagear_photoeccentric/docs/source/tutorials")
    rows = [
        legacy_record(root / "tutorial01/818.01", "818.01"),
        legacy_record(root / "tutorial02/254.01", "254.01"),
    ]
    fits_paths = sorted(Path(args.fits).glob("*/**/*-results.fits"))
    for path in fits_paths[:10]:
        rows.append(current_record(path))

    frame = pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
