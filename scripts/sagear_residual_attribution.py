"""Locate where the Sagear population mismatch enters the recovered sample.

This audit keeps sample membership, radius definitions, and post-fit quality
cuts separate.  It is intentionally diagnostic: the recovered 2,465-planet
set has the published total but cannot reproduce the unpublished visual-QC
identities or final four-way category counts.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from gilbert_real_alderaan_control import source_prefit_mask
from hierarchical_rayleigh import fit_rayleigh, load_population_masses
from hierarchical_table3_order_diagnostic import fit_model, starts_for


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUMMARY = (
    ROOT
    / "tmp"
    / "exact_inventory_berger2020_fixed_dynesty_pairedperiod_rms_periastron_sir1000_20260807"
    / "exact_inventory_summary_density_draw.csv"
)
DEFAULT_MANIFEST = (
    ROOT / "outputs" / "eccentricity_posterior_manifest_equal_nested_published_inventory_pre_visual_qc.csv"
)
DEFAULT_GILBERT = (
    ROOT / "tmp" / "gjgilbert_kepler_ecc_rp" / "Catalogs" / "kepler_dr25_gaia_dr2_crossmatch.csv"
)
DEFAULT_POSTFIT = ROOT / "outputs" / "gilbert_real_alderaan_control" / "planet_qc_ledger.csv"

PAPER_FINAL_COUNTS = {
    ("thick", "single"): 275,
    ("thin", "single"): 1121,
    ("thick", "multi"): 207,
    ("thin", "multi"): 862,
}
PAPER_SMALL_COUNTS = {
    ("thick", "single"): 238,
    ("thin", "single"): 961,
    ("thick", "multi"): 199,
    ("thin", "multi"): 786,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--gilbert-catalog", type=Path, default=DEFAULT_GILBERT)
    parser.add_argument("--postfit-ledger", type=Path, default=DEFAULT_POSTFIT)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "outputs" / "sagear_residual_attribution"
    )
    return parser.parse_args()


def require_unique(frame: pd.DataFrame, key: str, label: str) -> None:
    if key not in frame:
        raise ValueError(f"{label} is missing {key}")
    if frame[key].duplicated().any():
        raise ValueError(f"{label} contains duplicate {key} values")


def group_count_rows(frame: pd.DataFrame, source: str, radius_definition: str) -> list[dict]:
    rows: list[dict] = []
    counts = frame.groupby(["disk", "system"]).size().to_dict()
    for group, paper_count in PAPER_SMALL_COUNTS.items():
        rows.append(
            {
                "source": source,
                "radius_definition": radius_definition,
                "disk": group[0],
                "system": group[1],
                "count": int(counts.get(group, 0)),
                "paper_small_count": paper_count,
                "difference": int(counts.get(group, 0)) - paper_count,
            }
        )
    return rows


def fit_stage(stage: str, frame: pd.DataFrame) -> list[dict]:
    groups = {
        "all": frame,
        "thick_singles": frame[(frame.disk == "thick") & (frame.system == "single")],
        "thin_singles": frame[(frame.disk == "thin") & (frame.system == "single")],
        "thick_multis": frame[(frame.disk == "thick") & (frame.system == "multi")],
        "thin_multis": frame[(frame.disk == "thin") & (frame.system == "multi")],
    }
    rows: list[dict] = []
    sigma_grid = np.linspace(1e-4, 1.0, 2000)
    for population, group in groups.items():
        if len(group) < 2:
            continue
        rayleigh = fit_rayleigh(
            group.posterior_file.tolist(),
            sigma_grid,
            apply_transit_selection=True,
            selection_mode="manuscript_reciprocal",
        )
        rows.append(
            {
                "stage": stage,
                "population": population,
                "model": "rayleigh",
                "n_planets": len(group),
                "estimate": rayleigh["expected_e"],
                "lower": rayleigh["expected_e_lo"],
                "upper": rayleigh["expected_e_hi"],
            }
        )
        masses, e_grid = load_population_masses(
            group, apply_transit_selection=True, selection_mode="manuscript_reciprocal"
        )
        beta = fit_model(masses, e_grid, "beta", starts_for("beta"))
        rows.append(
            {
                "stage": stage,
                "population": population,
                "model": "beta_map",
                "n_planets": len(group),
                "estimate": beta["mean_e"],
                "lower": np.nan,
                "upper": np.nan,
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    for path in (args.summary, args.manifest, args.gilbert_catalog, args.postfit_ledger):
        if not path.exists():
            raise FileNotFoundError(path)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary = pd.read_csv(args.summary)
    manifest = pd.read_csv(args.manifest)
    gilbert = pd.read_csv(args.gilbert_catalog, index_col=0)
    postfit = pd.read_csv(args.postfit_ledger)
    require_unique(summary, "kepoi_name", "posterior summary")
    require_unique(manifest, "kepoi_name", "pre-visual manifest")
    require_unique(gilbert, "planet_name", "Gilbert catalog")
    require_unique(postfit, "kepoi_name", "Gilbert post-fit ledger")
    if len(summary) != 2465:
        raise ValueError(f"Expected 2,465 recovered posterior rows, found {len(summary)}")

    manifest_columns = [
        "kepoi_name", "koi_prad", "koi_ror", "berger_rad", "disk", "system"
    ]
    joined = summary.merge(
        manifest[manifest_columns],
        on="kepoi_name",
        how="left",
        validate="one_to_one",
        suffixes=("", "_manifest"),
    )
    if joined.koi_prad.isna().any():
        raise ValueError("Recovered summary contains planets absent from the pre-visual manifest")

    gilbert = gilbert.copy()
    gilbert["gilbert_small_prefit"] = source_prefit_mask(gilbert, small_planets=True)
    joined = joined.merge(
        gilbert[["planet_name", "gilbert_small_prefit", "rp"]],
        left_on="kepoi_name",
        right_on="planet_name",
        how="left",
        validate="one_to_one",
    )
    postfit_flags = [
        column
        for column in postfit.columns
        if column == "gilbert_control_exclude" or column.startswith("exclude_")
    ]
    joined = joined.merge(
        postfit[["kepoi_name", *postfit_flags]],
        on="kepoi_name",
        how="left",
        validate="one_to_one",
    )
    joined["rp_koi"] = pd.to_numeric(joined.koi_prad, errors="coerce")
    joined["rp_berger_koi_ror"] = (
        pd.to_numeric(joined.koi_ror, errors="coerce")
        * pd.to_numeric(joined.berger_rad, errors="coerce")
        * 109.076
    )

    recovered_counts = joined.groupby(["disk", "system"]).size().to_dict()
    category_rows = []
    for group, expected in PAPER_FINAL_COUNTS.items():
        category_rows.append(
            {
                "disk": group[0],
                "system": group[1],
                "recovered_count": int(recovered_counts.get(group, 0)),
                "paper_final_count": expected,
                "difference": int(recovered_counts.get(group, 0)) - expected,
            }
        )
    category_counts = pd.DataFrame(category_rows)

    radius_rows: list[dict] = []
    radius_rows += group_count_rows(
        manifest[manifest.koi_prad.lt(3.5)], "pre_visual_manifest", "DR25 koi_prad < 3.5"
    )
    derived = manifest.koi_ror * manifest.berger_rad * 109.076
    radius_rows += group_count_rows(
        manifest[derived.lt(3.5)], "pre_visual_manifest", "DR25 koi_ror x Berger Rstar < 3.5"
    )
    radius_rows += group_count_rows(
        joined[joined.rp_koi.lt(3.5)], "recovered_2465", "DR25 koi_prad < 3.5"
    )
    radius_rows += group_count_rows(
        joined[joined.rp_berger_koi_ror.lt(3.5)],
        "recovered_2465",
        "DR25 koi_ror x Berger Rstar < 3.5",
    )
    radius_counts = pd.DataFrame(radius_rows)

    prefit_mask = joined.gilbert_small_prefit.fillna(False).astype(bool)
    sequential = prefit_mask.copy()
    stages = {
        "recovered_all": joined,
        "recovered_koi_prad_lt_3p5": joined[joined.rp_koi.lt(3.5)],
        "recovered_berger_ror_lt_3p5": joined[joined.rp_berger_koi_ror.lt(3.5)],
        "gilbert_small_prefit_intersection": joined[prefit_mask],
        "gilbert_small_postfit_intersection": joined[
            joined.gilbert_control_exclude.eq(False)
        ],
    }
    for flag in (
        "exclude_qc_unavailable",
        "exclude_grazing",
        "exclude_planet_radius_precision",
        "exclude_density_asymmetry",
        "exclude_importance_viability",
    ):
        if flag not in joined:
            raise ValueError(f"Post-fit ledger is missing {flag}")
        sequential &= ~joined[flag].fillna(True).astype(bool)
        stages[f"gilbert_cumulative_through_{flag.removeprefix('exclude_')}"] = joined[
            sequential
        ]
    distribution_rows = []
    for stage, frame in stages.items():
        for (disk, system), group in frame.groupby(["disk", "system"]):
            distribution_rows.append(
                {
                    "stage": stage,
                    "disk": disk,
                    "system": system,
                    "n_planets": len(group),
                    "median_planet_e50": group.e50.median(),
                    "mean_planet_e50": group.e50.mean(),
                    "fraction_planet_e50_gt_0p2": group.e50.gt(0.2).mean(),
                }
            )
    distribution = pd.DataFrame(distribution_rows)

    hierarchy_rows: list[dict] = []
    for stage, frame in stages.items():
        hierarchy_rows.extend(fit_stage(stage, frame))
    hierarchy = pd.DataFrame(hierarchy_rows)

    category_counts.to_csv(args.output_dir / "recovered_category_counts.csv", index=False)
    radius_counts.to_csv(args.output_dir / "small_planet_radius_count_audit.csv", index=False)
    distribution.to_csv(args.output_dir / "posterior_distribution_by_stage.csv", index=False)
    hierarchy.to_csv(args.output_dir / "hierarchical_attribution.csv", index=False)
    joined.to_csv(args.output_dir / "planet_stage_ledger.csv", index=False)

    lines = [
        "# Sagear residual attribution",
        "",
        "The recovered set has the published total but not the unpublished final membership. "
        "Its category differences are reported rather than silently reconciled.",
        "",
        "## Recovered category counts",
        "",
        category_counts.to_markdown(index=False),
        "",
        "## Small-planet radius definitions",
        "",
        radius_counts.to_markdown(index=False),
        "",
        "## Hierarchical cut ladder",
        "",
        hierarchy.to_markdown(index=False),
        "",
        "The Gilbert post-fit intersection is an external positive control, not a substitute "
        "for Sagear's unpublished visual-QC ledger.",
    ]
    (args.output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print(category_counts.to_string(index=False))
    print(hierarchy.to_string(index=False))
    print(f"WROTE {args.output_dir / 'report.md'}")


if __name__ == "__main__":
    main()
