from __future__ import annotations

from pathlib import Path

import pandas as pd


PUBLISHED = {
    "thick_singles": (275, 0.066, 0.045, 0.096),
    "thin_singles": (1121, 0.022, 0.017, 0.029),
    "thick_multis": (207, 0.033, 0.015, 0.065),
    "thin_multis": (862, 0.030, 0.023, 0.031),
}

ORDER = ["thick_singles", "thin_singles", "thick_multis", "thin_multis"]


BRANCHES = [
    (
        "dynesty_weighted_exact_inventory",
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "exact_dynesty_weighted_20260724.csv",
        "Statistically weighted ALDERAAN transit posterior; exact reconstructed inventory.",
    ),
    (
        "sagear_literal_equal_raw_exact_inventory",
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "exact_equal_raw_diagnostic_20260724.csv",
        "Literal source-comment diagnostic: raw ALDERAAN rows treated equally; same inventory.",
    ),
]


def main() -> None:
    root = Path(__file__).resolve().parent
    output_dir = root / "outputs"
    rows: list[pd.DataFrame] = []

    for method, filename, description in BRANCHES:
        path = output_dir / filename
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path)
        required = [
            "population",
            "n",
            "expected_e",
            "expected_e_lo",
            "expected_e_hi",
            "boundary_flag",
        ]
        missing = [column for column in required if column not in frame.columns]
        if missing:
            raise ValueError(f"{path} is missing columns: {missing}")
        frame = frame[required].copy()
        frame.insert(0, "method", method)
        frame.insert(1, "description", description)
        rows.append(frame)

    published = pd.DataFrame(
        [
            {
                "population": population,
                "sagear_n": values[0],
                "sagear_expected_e": values[1],
                "sagear_expected_e_lo": values[2],
                "sagear_expected_e_hi": values[3],
            }
            for population, values in PUBLISHED.items()
        ]
    )
    comparison = pd.concat(rows, ignore_index=True).merge(
        published, on="population", validate="many_to_one"
    )
    comparison["ratio_to_sagear"] = comparison["expected_e"] / comparison["sagear_expected_e"]
    comparison["absolute_difference"] = (
        comparison["expected_e"] - comparison["sagear_expected_e"]
    ).abs()
    comparison["one_sigma_intervals_overlap"] = (
        comparison["expected_e_hi"].ge(comparison["sagear_expected_e_lo"])
        & comparison["expected_e_lo"].le(comparison["sagear_expected_e_hi"])
    )
    comparison["population"] = pd.Categorical(
        comparison["population"], ORDER, ordered=True
    )
    comparison = comparison.sort_values(["method", "population"])

    csv_path = output_dir / "sagear_replication_contract_comparison_exact_inventory.csv"
    comparison.to_csv(csv_path, index=False)

    lines = [
        "# Sagear Replication Contract Comparison",
        "",
        "Both branches use the same exact reconstructed inventory family and the same",
        "paired-impact, exact-duration, Berger-2020, and reciprocal-selection settings.",
        "Branch-specific valid and QC exclusions can therefore produce small N differences.",
        "They differ only in how ALDERAAN nested rows are exported into the",
        "post-model calculation.",
        "",
        "| branch | population | N | expected e (16-84%) | Sagear | ratio | overlap |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for row in comparison.itertuples(index=False):
        lines.append(
            f"| {row.method} | {row.population} | {row.n} | "
            f"{row.expected_e:.3f} ({row.expected_e_lo:.3f}-{row.expected_e_hi:.3f}) | "
            f"{row.sagear_expected_e:.3f} "
            f"({row.sagear_expected_e_lo:.3f}-{row.sagear_expected_e_hi:.3f}) | "
            f"{row.ratio_to_sagear:.2f} | "
            f"{'yes' if row.one_sigma_intervals_overlap else 'no'} |"
        )
    (output_dir / "sagear_replication_contract_comparison_exact_inventory.md").write_text(
        "\n".join(lines) + "\n", encoding="ascii"
    )


if __name__ == "__main__":
    main()
