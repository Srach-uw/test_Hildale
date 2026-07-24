from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SAGEAR = {
    "thick_singles": (0.066, 0.045, 0.096),
    "thin_singles": (0.022, 0.017, 0.029),
    "thick_multis": (0.033, 0.015, 0.065),
    "thin_multis": (0.030, 0.023, 0.031),
}


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def comparison_rows(
    fit: pd.DataFrame, hypothesis: str
) -> pd.DataFrame:
    required = {
        "population",
        "n",
        "expected_e",
        "expected_e_lo",
        "expected_e_hi",
        "boundary_flag",
    }
    require_columns(fit, required, hypothesis)
    if set(fit["population"]) != set(SAGEAR):
        raise ValueError(
            f"{hypothesis} populations differ from Sagear: "
            f"{sorted(fit['population'].astype(str).unique())}"
        )
    out = fit.copy()
    out["hypothesis"] = hypothesis
    out["sagear_expected_e"] = out["population"].map(
        lambda population: SAGEAR[population][0]
    )
    out["sagear_lo"] = out["population"].map(
        lambda population: SAGEAR[population][1]
    )
    out["sagear_hi"] = out["population"].map(
        lambda population: SAGEAR[population][2]
    )
    out["delta"] = out["expected_e"] - out["sagear_expected_e"]
    out["ratio"] = out["expected_e"] / out["sagear_expected_e"]
    out["interval_overlap"] = (
        (out["expected_e_hi"] >= out["sagear_lo"])
        & (out["expected_e_lo"] <= out["sagear_hi"])
    )

    ours_sigma = 0.5 * (out["expected_e_hi"] - out["expected_e_lo"])
    paper_sigma = np.where(
        out["delta"] >= 0,
        out["sagear_hi"] - out["sagear_expected_e"],
        out["sagear_expected_e"] - out["sagear_lo"],
    )
    out["combined_sigma_approx"] = np.sqrt(ours_sigma**2 + paper_sigma**2)
    out["standardized_residual_approx"] = (
        out["delta"] / out["combined_sigma_approx"]
    )
    return out


def score_hypotheses(comparison: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for hypothesis, group in comparison.groupby("hypothesis", sort=False):
        rows.append(
            {
                "hypothesis": hypothesis,
                "populations": len(group),
                "interval_overlaps": int(group["interval_overlap"].sum()),
                "rmse_expected_e": float(
                    np.sqrt(np.mean(np.square(group["delta"])))
                ),
                "mean_absolute_error": float(np.mean(np.abs(group["delta"]))),
                "approx_chi_square": float(
                    np.sum(np.square(group["standardized_residual_approx"]))
                ),
                "max_absolute_standardized_residual": float(
                    np.max(np.abs(group["standardized_residual_approx"]))
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("rmse_expected_e")


def make_markdown(
    count_audit: pd.DataFrame,
    comparison: pd.DataFrame,
    scores: pd.DataFrame,
    posterior_rows: int,
    posterior_hosts: int,
) -> str:
    total = count_audit.loc[
        (count_audit["disk"] == "all") & (count_audit["system"] == "all")
    ].iloc[0]
    lines = [
        "# Exact-inventory Sagear replication audit",
        "",
        "## Inventory reconstruction",
        "",
        "Joining the published 1,888-host table to non-false-positive KOIs with",
        "periods from 1 to 100 days, while preserving the pre-cut `koi_count`,",
        "produces 2,491 pre-visual-QC planets. The paper contains 2,465.",
        "",
        "| disk | system | pre-visual | paper | implied final removals |",
        "|---|---|---:|---:|---:|",
    ]
    for _, row in count_audit.loc[count_audit["disk"] != "all"].iterrows():
        lines.append(
            f"| {row['disk']} | {row['system']} | "
            f"{int(row['pre_visual_planets'])} | {int(row['paper_planets'])} | "
            f"{int(row['implied_final_rejections'])} |"
        )
    lines.extend(
        [
            "",
            f"The available ALDERAAN FITS provide {posterior_rows} planet",
            f"posteriors on {posterior_hosts} published hosts. The two single",
            "categories already equal the paper counts before the unpublished",
            "final rejection stage. The 26 count-implied removals are 21 thin",
            "multis and 5 thick multis.",
            "",
            "## Population-label hypotheses",
            "",
            "| hypothesis | paper population | N | ours (16th-84th) | paper | "
            "ratio | interval overlap |",
            "|---|---|---:|---:|---:|---:|---|",
        ]
    )
    for _, row in comparison.sort_values(["hypothesis", "population"]).iterrows():
        lines.append(
            f"| {row['hypothesis']} | {row['population']} | {int(row['n'])} | "
            f"{row['expected_e']:.4f} "
            f"[{row['expected_e_lo']:.4f}, {row['expected_e_hi']:.4f}] | "
            f"{row['sagear_expected_e']:.3f} | {row['ratio']:.2f} | "
            f"{'yes' if row['interval_overlap'] else 'no'} |"
        )
    lines.extend(
        [
            "",
            "## Aggregate agreement",
            "",
            "| hypothesis | overlaps / 4 | RMSE | mean absolute error | "
            "approx chi-square | max absolute standardized residual |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in scores.iterrows():
        lines.append(
            f"| {row['hypothesis']} | {int(row['interval_overlaps'])} / 4 | "
            f"{row['rmse_expected_e']:.4f} | "
            f"{row['mean_absolute_error']:.4f} | "
            f"{row['approx_chi_square']:.2f} | "
            f"{row['max_absolute_standardized_residual']:.2f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The inverted-disk diagnostics are not scientifically permissible",
            "relabeling. Published Table 1 explicitly defines Pthick and its disk",
            "association, and the paper text and figures consistently state that",
            "thick-disk singles are hotter. The diagnostic instead tests whether",
            "the posterior arrays or group labels may have been inverted between",
            "Table 1 and the population-analysis code. In combination with equal",
            "treatment of raw dynesty rows, the inversion reproduces all four",
            "published Rayleigh intervals. The corresponding dynesty-weighted",
            "controls remain far from the paper whether labels are swapped or not.",
            "",
            "This is strong evidence for a possible implementation-level label",
            "orientation error, but it is not proof. Proof requires the authors'",
            "planet-level inclusion table, post-model posterior files, and",
            "population-analysis code.",
            "",
            f"Count contract: {int(total['pre_visual_planets'])} pre-visual planets, "
            f"{int(total['paper_planets'])} published planets, "
            f"{int(total['implied_final_rejections'])} implied removals.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare exact published-host replication hypotheses."
    )
    parser.add_argument("--outputs", default="outputs")
    args = parser.parse_args()
    out = Path(args.outputs).expanduser().resolve()

    count_path = out / "sagear2026_planet_inventory_count_audit.csv"
    summary_path = (
        out
        / "eccentricity_posterior_summary_equal_nested_published_inventory_pre_visual_qc.csv"
    )
    equal_normal_path = (
        out
        / "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "EQUAL_NESTED_PUBLISHED_INVENTORY_PREVISUAL.csv"
    )
    equal_swap_path = (
        out
        / "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "EQUAL_NESTED_PUBLISHED_INVENTORY_DISK_SWAP_DIAGNOSTIC.csv"
    )
    dynesty_normal_path = (
        out
        / "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "DYNESTY_PUBLISHED_INVENTORY_PREVISUAL.csv"
    )
    dynesty_swap_path = (
        out
        / "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "DYNESTY_PUBLISHED_INVENTORY_DISK_SWAP_DIAGNOSTIC.csv"
    )
    for path in (
        count_path,
        summary_path,
        equal_normal_path,
        equal_swap_path,
        dynesty_normal_path,
        dynesty_swap_path,
    ):
        if not path.exists():
            raise FileNotFoundError(path)

    count_audit = pd.read_csv(count_path)
    summary = pd.read_csv(summary_path)
    equal_normal = comparison_rows(
        pd.read_csv(equal_normal_path),
        "equal_rows_published_labels",
    )
    equal_swapped = comparison_rows(
        pd.read_csv(equal_swap_path),
        "equal_rows_inverted_disk_diagnostic",
    )
    dynesty_normal = comparison_rows(
        pd.read_csv(dynesty_normal_path),
        "dynesty_weights_published_labels",
    )
    dynesty_swapped = comparison_rows(
        pd.read_csv(dynesty_swap_path),
        "dynesty_weights_inverted_disk_diagnostic",
    )
    comparison = pd.concat(
        [equal_normal, equal_swapped, dynesty_normal, dynesty_swapped],
        ignore_index=True,
    )
    scores = score_hypotheses(comparison)

    comparison.to_csv(
        out / "sagear2026_replication_hypothesis_comparison.csv", index=False
    )
    scores.to_csv(
        out / "sagear2026_replication_hypothesis_scores.csv", index=False
    )
    report = make_markdown(
        count_audit,
        comparison,
        scores,
        posterior_rows=len(summary),
        posterior_hosts=summary["kepid"].nunique(),
    )
    (out / "sagear2026_exact_inventory_replication_audit.md").write_text(
        report, encoding="utf-8"
    )
    print(scores.to_string(index=False))


if __name__ == "__main__":
    main()
