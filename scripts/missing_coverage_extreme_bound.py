from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from hierarchical_rayleigh import (
    POPULATIONS,
    fit_from_mass_matrix,
    load_population_masses,
    transit_probability_weight,
)


PAPER_COUNTS = {
    "thick_singles": 275,
    "thin_singles": 1121,
    "thick_multis": 207,
    "thin_multis": 862,
}
PAPER_EXPECTED_E = {
    "thick_singles": 0.066,
    "thin_singles": 0.022,
    "thick_multis": 0.033,
    "thin_multis": 0.030,
}


def synthetic_near_circular_mass(
    e_grid: np.ndarray,
    omega_grid: np.ndarray,
    scale: float,
    selection_mode: str,
) -> np.ndarray:
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be positive and finite")
    e_pdf = np.exp(-0.5 * np.square(e_grid / scale))
    e_pdf /= e_pdf.sum()
    posterior = np.repeat(
        (e_pdf / len(omega_grid))[:, None], len(omega_grid), axis=1
    )
    if selection_mode == "manuscript_reciprocal":
        return np.sum(
            posterior / transit_probability_weight(e_grid, omega_grid), axis=1
        )
    if selection_mode == "none":
        return posterior.sum(axis=1)
    raise ValueError(f"unsupported diagnostic selection mode: {selection_mode}")


def run_bounds(
    summary: pd.DataFrame,
    scales: list[float],
    selection_mode: str = "manuscript_reciprocal",
) -> pd.DataFrame:
    sigmas = np.linspace(1e-4, 1.0, 2000)
    rows: list[dict[str, object]] = []
    for disk, system, population in POPULATIONS:
        sub = summary.loc[
            (summary["disk"] == disk)
            & (summary["system"] == system)
            & ~summary["qc_primary_exclude"].fillna(True).astype(bool)
        ].reset_index(drop=True)
        masses, e_grid = load_population_masses(
            sub,
            apply_transit_selection=True,
            selection_mode=selection_mode,
        )
        with np.load(sub["posterior_file"].iloc[0], allow_pickle=False) as data:
            omega_grid = np.asarray(data["omega_grid"], dtype=float)
        n_missing = PAPER_COUNTS[population] - len(sub)
        if n_missing < 0:
            raise ValueError(
                f"{population} has {len(sub)} rows, above paper count "
                f"{PAPER_COUNTS[population]}"
            )
        baseline = fit_from_mass_matrix(
            masses,
            e_grid,
            sigmas,
            apply_transit_selection=True,
            selection_mode=selection_mode,
        )
        rows.append(
            {
                "population": population,
                "scenario": "observed_available_only",
                "near_circular_scale": np.nan,
                "n_observed": len(sub),
                "n_synthetic_missing": 0,
                "n_total": len(sub),
                "expected_e": baseline["expected_e"],
                "expected_e_lo": baseline["expected_e_lo"],
                "expected_e_hi": baseline["expected_e_hi"],
                "sagear_expected_e": PAPER_EXPECTED_E[population],
            }
        )
        for scale in scales:
            synthetic = synthetic_near_circular_mass(
                e_grid, omega_grid, scale, selection_mode
            )
            filled = np.vstack(
                [masses, np.repeat(synthetic[None, :], n_missing, axis=0)]
            )
            fit = fit_from_mass_matrix(
                filled,
                e_grid,
                sigmas,
                apply_transit_selection=True,
                selection_mode=selection_mode,
            )
            rows.append(
                {
                    "population": population,
                    "scenario": "all_missing_artificially_near_circular",
                    "near_circular_scale": scale,
                    "n_observed": len(sub),
                    "n_synthetic_missing": n_missing,
                    "n_total": len(filled),
                    "expected_e": fit["expected_e"],
                    "expected_e_lo": fit["expected_e_lo"],
                    "expected_e_hi": fit["expected_e_hi"],
                    "sagear_expected_e": PAPER_EXPECTED_E[population],
                }
            )
    out = pd.DataFrame(rows)
    out["delta_from_sagear"] = out["expected_e"] - out["sagear_expected_e"]
    out["reaches_sagear_central_or_lower"] = out["expected_e"] <= out[
        "sagear_expected_e"
    ]
    return out


def write_markdown(frame: pd.DataFrame, path: Path) -> None:
    lines = [
        "# Missing-posterior extreme-bound diagnostic",
        "",
        "This is a deliberately nonphysical sensitivity bound. Every unavailable",
        "planet needed to reach the paper's final category count is assigned the",
        "same extremely near-circular posterior. It asks whether missing coverage",
        "could possibly account for the residual without claiming that these are",
        "the planets' real posteriors.",
        "",
        "| population | observed / paper | synthetic scale | fitted mean e | "
        "Sagear mean e | reaches Sagear center |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for _, row in frame.iterrows():
        scale = (
            "none"
            if pd.isna(row["near_circular_scale"])
            else f"{row['near_circular_scale']:.4f}"
        )
        lines.append(
            f"| {row['population']} | {int(row['n_observed'])} / "
            f"{int(row['n_total'])} | {scale} | {row['expected_e']:.4f} | "
            f"{row['sagear_expected_e']:.3f} | "
            f"{'yes' if row['reaches_sagear_central_or_lower'] else 'no'} |"
        )
    lines.extend(
        [
            "",
            "A positive result would only show mathematical possibility. A negative",
            "result rules out incomplete posterior coverage under an assumption",
            "that is already more favorable to low eccentricity than realistic",
            "ALDERAAN posteriors.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bound the effect of missing posterior coverage."
    )
    parser.add_argument(
        "--summary",
        default=(
            "outputs/eccentricity_posterior_summary_equal_nested_"
            "published_inventory_pre_visual_qc.csv"
        ),
    )
    parser.add_argument(
        "--scales",
        default="0.0025,0.005,0.01,0.02",
        help="Comma-separated half-Gaussian eccentricity scales.",
    )
    parser.add_argument(
        "--out",
        default="outputs/missing_coverage_extreme_bound.csv",
    )
    args = parser.parse_args()

    scales = [float(value) for value in args.scales.split(",")]
    summary = pd.read_csv(args.summary)
    frame = run_bounds(summary, scales)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)
    write_markdown(frame, out.with_suffix(".md"))
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
