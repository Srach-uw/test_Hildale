"""Plot the versioned Table 3 comparison used in the public reconstruction."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


LABELS = {
    "thin_singles": "Thin singles",
    "thick_singles": "Thick singles",
    "thin_multis": "Thin multis",
    "thick_multis": "Thick multis",
}
ORDER = tuple(LABELS)


def load_comparison(path: Path) -> pd.DataFrame:
    """Read and validate the compact population comparison table."""
    table = pd.read_csv(path)
    required = {
        "population",
        "qc_n",
        "qc_expected_e",
        "qc_expected_e_lo",
        "qc_expected_e_hi",
        "paper_n",
        "paper_expected_e",
        "paper_expected_e_lo",
        "paper_expected_e_hi",
    }
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"comparison table missing columns: {missing}")
    table = table.set_index("population").reindex(ORDER).reset_index()
    if table["population"].isna().any():
        raise ValueError("comparison table is missing a required population")
    return table


def make_plot(table: pd.DataFrame, output: Path) -> None:
    """Draw paper and reconstruction Rayleigh means with intervals."""
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.8), constrained_layout=True)
    for axis, (_, row) in zip(axes.flat, table.iterrows()):
        paper = float(row["paper_expected_e"])
        recovered = float(row["qc_expected_e"])
        axis.errorbar(
            paper,
            0,
            xerr=[
                [paper - float(row["paper_expected_e_lo"])],
                [float(row["paper_expected_e_hi"]) - paper],
            ],
            fmt="o",
            color="#9b1b1b",
            capsize=3,
            label="Sagear et al. (2026)",
        )
        axis.errorbar(
            recovered,
            1,
            xerr=[
                [recovered - float(row["qc_expected_e_lo"])],
                [float(row["qc_expected_e_hi"]) - recovered],
            ],
            fmt="o",
            color="#007c83",
            capsize=3,
            label="Public reconstruction",
        )
        axis.set_title(f"{LABELS[row['population']]}\n(n={int(row['qc_n'])}; paper n={int(row['paper_n'])})")
        axis.set_xlim(0, 0.28)
        axis.set_yticks([0, 1], ["Paper", "Reconstruction"])
        axis.set_xlabel("Rayleigh mean eccentricity")
        axis.grid(axis="x", alpha=0.25)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=repo / "metadata/public_reconstruction_20260727/population_comparison.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=repo / "figures/population_comparison.png",
    )
    args = parser.parse_args()
    table = load_comparison(args.input)
    make_plot(table, args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
