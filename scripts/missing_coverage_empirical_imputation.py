from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from hierarchical_rayleigh import (
    POPULATIONS,
    load_population_masses,
    posterior_weights_from_ll,
    rayleigh_grid,
    weighted_quantile,
)
from missing_coverage_extreme_bound import PAPER_COUNTS, PAPER_EXPECTED_E


def duplicate_analog_results(
    summary: pd.DataFrame,
    selection_mode: str = "manuscript_reciprocal",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sigmas = np.linspace(1e-4, 1.0, 2000)
    expected_grid = sigmas * np.sqrt(np.pi / 2.0)
    detail_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

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
        rays, _ = rayleigh_grid(
            e_grid,
            sigmas,
            apply_transit_selection=True,
            selection_mode=selection_mode,
            outlier_floor=0.0,
        )
        log_terms = np.log(np.clip(masses @ rays, 1e-300, None))
        full_ll = log_terms.sum(axis=0)
        n_missing = PAPER_COUNTS[population] - len(sub)
        paper_value = PAPER_EXPECTED_E[population]

        values = []
        for index in range(len(sub)):
            weights = posterior_weights_from_ll(
                sigmas, full_ll + n_missing * log_terms[index]
            )
            expected_e = weighted_quantile(
                expected_grid, weights, [0.5]
            )[0]
            values.append(expected_e)
            row = sub.iloc[index]
            detail_rows.append(
                {
                    "population": population,
                    "analog_kepoi_name": row["kepoi_name"],
                    "analog_koi_target": row.get("koi_target", ""),
                    "analog_e50": row.get("e50", np.nan),
                    "analog_zeta_median": row.get("zeta_median", np.nan),
                    "analog_posterior_source": row.get(
                        "posterior_source", ""
                    ),
                    "n_observed": len(sub),
                    "n_missing_duplicated": n_missing,
                    "expected_e_after_imputation": expected_e,
                    "sagear_expected_e": paper_value,
                    "absolute_error": abs(expected_e - paper_value),
                }
            )
        values_array = np.asarray(values)
        closest_index = int(np.argmin(np.abs(values_array - paper_value)))
        closest = sub.iloc[closest_index]
        summary_rows.append(
            {
                "population": population,
                "n_observed": len(sub),
                "paper_count": PAPER_COUNTS[population],
                "n_missing": n_missing,
                "imputed_expected_e_min": float(values_array.min()),
                "imputed_expected_e_max": float(values_array.max()),
                "sagear_expected_e": paper_value,
                "sagear_within_empirical_analog_range": bool(
                    values_array.min() <= paper_value <= values_array.max()
                ),
                "closest_analog_kepoi_name": closest["kepoi_name"],
                "closest_analog_e50": closest.get("e50", np.nan),
                "closest_imputed_expected_e": float(
                    values_array[closest_index]
                ),
                "closest_absolute_error": float(
                    abs(values_array[closest_index] - paper_value)
                ),
            }
        )
    return pd.DataFrame(summary_rows), pd.DataFrame(detail_rows)


def write_markdown(frame: pd.DataFrame, path: Path) -> None:
    lines = [
        "# Missing-posterior empirical analog diagnostic",
        "",
        "For each population, every missing posterior is replaced by copies of",
        "one real observed posterior from the same population. Repeating this",
        "over every available planet gives an empirical extreme range. This is",
        "not an imputation estimate; it tests whether observed posterior shapes",
        "are capable of moving the incomplete population to Sagear's value.",
        "",
        "| population | observed / paper | analog-imputed range | Sagear | "
        "Sagear inside range | closest real analog | analog e50 |",
        "|---|---:|---:|---:|---|---|---:|",
    ]
    for _, row in frame.iterrows():
        lines.append(
            f"| {row['population']} | {int(row['n_observed'])} / "
            f"{int(row['paper_count'])} | "
            f"{row['imputed_expected_e_min']:.4f}-"
            f"{row['imputed_expected_e_max']:.4f} | "
            f"{row['sagear_expected_e']:.3f} | "
            f"{'yes' if row['sagear_within_empirical_analog_range'] else 'no'} | "
            f"{row['closest_analog_kepoi_name']} | "
            f"{row['closest_analog_e50']:.3f} |"
        )
    lines.extend(
        [
            "",
            "If Sagear lies outside this range, missing coverage cannot explain",
            "the bin even under the extreme assumption that all missing systems",
            "resemble the most favorable observed analog. If it lies inside, the",
            "missing systems must be fitted before attributing the residual to a",
            "group-label implementation error.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Empirically bound missing-posterior effects."
    )
    parser.add_argument(
        "--summary",
        default=(
            "outputs/eccentricity_posterior_summary_equal_nested_"
            "published_inventory_pre_visual_qc.csv"
        ),
    )
    parser.add_argument(
        "--out",
        default="outputs/missing_coverage_empirical_imputation.csv",
    )
    args = parser.parse_args()

    summary, details = duplicate_analog_results(pd.read_csv(args.summary))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out, index=False)
    details.to_csv(
        out.with_name(f"{out.stem}_all_analogs.csv"), index=False
    )
    write_markdown(summary, out.with_suffix(".md"))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
