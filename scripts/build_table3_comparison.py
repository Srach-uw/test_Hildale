from __future__ import annotations

from pathlib import Path

import pandas as pd


PUBLISHED = {
    "thick_singles": (275, 0.066, 0.045, 0.096),
    "thin_singles": (1121, 0.022, 0.017, 0.029),
    "thick_multis": (207, 0.033, 0.015, 0.065),
    "thin_multis": (862, 0.030, 0.023, 0.031),
}


def load_fit(path: Path, method: str) -> pd.DataFrame:
    columns = ["population", "n", "expected_e", "expected_e_lo", "expected_e_hi", "boundary_flag"]
    frame = pd.read_csv(path)[columns].copy()
    frame.insert(0, "method", method)
    return frame


def main() -> None:
    root = Path(__file__).resolve().parent
    outputs = root / "outputs"
    fits = pd.concat(
        [
            load_fit(
                outputs
                / "rayleigh_population_fit_transit_selection_manuscript_reciprocal_UNIFORM_PAIRED_FULL_QC_PRIMARY.csv",
                "canonical_dynesty_weights",
            ),
            load_fit(
                outputs
                / "rayleigh_population_fit_transit_selection_manuscript_reciprocal_UNIFORM_EQUAL_NESTED_50K_STRICT.csv",
                "diagnostic_equal_nested_rows",
            ),
        ],
        ignore_index=True,
    )
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
    comparison = fits.merge(published, on="population", validate="many_to_one")
    comparison["ratio_to_sagear"] = comparison["expected_e"] / comparison["sagear_expected_e"]
    comparison["one_sigma_intervals_overlap"] = (
        comparison["expected_e_hi"].ge(comparison["sagear_expected_e_lo"])
        & comparison["expected_e_lo"].le(comparison["sagear_expected_e_hi"])
    )
    comparison["absolute_difference"] = (
        comparison["expected_e"] - comparison["sagear_expected_e"]
    ).abs()
    comparison.to_csv(outputs / "sagear_table3_rayleigh_comparison.csv", index=False)

    order = ["thick_singles", "thin_singles", "thick_multis", "thin_multis"]
    comparison["population"] = pd.Categorical(comparison["population"], order, ordered=True)
    comparison = comparison.sort_values(["method", "population"])
    lines = [
        "# Sagear Table 3 Rayleigh Comparison",
        "",
        "The equal-row result is a diagnostic reconstruction, not a valid dynesty",
        "posterior treatment. ALDERAAN stores raw nested samples and requires LN_WT.",
        "",
        "| method | population | N | mean e (16th-84th) | Sagear | ratio | 1-sigma overlap |",
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
    (outputs / "sagear_table3_rayleigh_comparison.md").write_text(
        "\n".join(lines) + "\n", encoding="ascii"
    )


if __name__ == "__main__":
    main()
