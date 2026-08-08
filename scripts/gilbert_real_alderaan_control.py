"""Run a source-defined Gilbert small-planet control on recovered ALDERAAN data.

This audit deliberately uses Gilbert et al.'s released analysis catalog for
the stellar and planet cuts.  It then joins those planets to the final Sagear
ALDERAAN posterior inventory and applies the post-fit cuts that can be
reproduced from the nested transit chains.  The result is an intersection
control, not a reconstruction of Gilbert's full 1,646-planet sample, because
the Sagear sample has additional kinematic and chemical-selection boundaries.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from hierarchical_rayleigh import fit_rayleigh, load_population_masses
from hierarchical_table3_order_diagnostic import fit_model, starts_for


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GILBERT_CATALOG = (
    ROOT
    / "tmp"
    / "gjgilbert_kepler_ecc_rp"
    / "Catalogs"
    / "kepler_dr25_gaia_dr2_crossmatch.csv"
)
DEFAULT_SUMMARY = (
    ROOT
    / "tmp"
    / "exact_inventory_berger2020_fixed_dynesty_pairedperiod_rms_periastron_sir1000_20260807"
    / "exact_inventory_summary_density_draw.csv"
)
DEFAULT_POSTFIT = ROOT / "outputs" / "gilbert_postfit_qc_audit_SOURCE_FAITHFUL_RMS_PERIASTRON_20260807.csv"
DEFAULT_VIABILITY = ROOT / "outputs" / "gilbert_importance_viability_SOURCE_FAITHFUL_RMS_PERIASTRON_20K_20260807.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Source-defined Gilbert small-planet control using real ALDERAAN posteriors."
    )
    parser.add_argument("--gilbert-catalog", type=Path, default=DEFAULT_GILBERT_CATALOG)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--postfit-qc", type=Path, default=DEFAULT_POSTFIT)
    parser.add_argument("--viability-audit", type=Path, default=DEFAULT_VIABILITY)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "gilbert_real_alderaan_control")
    parser.add_argument("--grazing-limit", type=float, default=0.05)
    return parser.parse_args()


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns - set(frame.columns)
    if missing:
        raise ValueError(f"{label} is missing required columns: {sorted(missing)}")


def source_prefit_mask(catalog: pd.DataFrame, *, small_planets: bool) -> pd.Series:
    """Reproduce the released Gilbert pre-fit cuts for the published domain."""

    stellar_radius_error = np.hypot(catalog["rstar_err1"], catalog["rstar_err2"]) / np.sqrt(2.0)
    radius_hi = 3.5 if small_planets else 16.0
    return (
        catalog["period"].between(1.0, 100.0, inclusive="both")
        & catalog["rp"].between(0.5, radius_hi, inclusive="both")
        & ((catalog["fpp"] < 0.1) | catalog["disposition"].eq("CONFIRMED"))
        & (catalog["logg"] > 4.0)
        & (((catalog["rcf"] - 1.0) < 0.05) | catalog["rcf"].isna())
        & (catalog["ruwe"] < 1.4)
        & ((stellar_radius_error / catalog["rstar"]) < 0.2)
        & catalog["rstar"].between(0.7, 1.4, inclusive="both")
        & catalog["teff"].between(4700.0, 6500.0, inclusive="both")
        & catalog["feh"].between(-0.6, 0.6, inclusive="both")
        & catalog["age"].between(0.0, 14.0, inclusive="both")
    )


def attrition(catalog: pd.DataFrame, recovered_names: set[str]) -> pd.DataFrame:
    steps: list[tuple[str, pd.Series]] = []
    rerr = np.hypot(catalog["rstar_err1"], catalog["rstar_err2"]) / np.sqrt(2.0)
    steps.extend(
        [
            ("DR25 catalog rows in released Gilbert catalog", pd.Series(True, index=catalog.index)),
            ("period 1-100 d", catalog["period"].between(1.0, 100.0)),
            ("planet radius 0.5-16 R_earth", catalog["rp"].between(0.5, 16.0)),
            ("FPP < 0.1 or confirmed", (catalog["fpp"] < 0.1) | catalog["disposition"].eq("CONFIRMED")),
            ("logg > 4", catalog["logg"] > 4.0),
            ("Furlan radius correction < 5% or unavailable", ((catalog["rcf"] - 1.0) < 0.05) | catalog["rcf"].isna()),
            ("RUWE < 1.4", catalog["ruwe"] < 1.4),
            ("stellar radius precision < 20%", (rerr / catalog["rstar"]) < 0.2),
            ("stellar radius 0.7-1.4 R_sun", catalog["rstar"].between(0.7, 1.4)),
            ("Teff 4700-6500 K", catalog["teff"].between(4700.0, 6500.0)),
            ("[Fe/H] -0.6 to 0.6", catalog["feh"].between(-0.6, 0.6)),
            ("age 0-14 Gyr", catalog["age"].between(0.0, 14.0)),
            ("small planets 0.5-3.5 R_earth", catalog["rp"].between(0.5, 3.5)),
            ("present in final Sagear ALDERAAN inventory", catalog["planet_name"].isin(recovered_names)),
        ]
    )
    keep = pd.Series(True, index=catalog.index)
    rows: list[dict[str, object]] = []
    for label, criterion in steps:
        before = int(keep.sum())
        keep &= criterion.fillna(False)
        after = int(keep.sum())
        rows.append({"cut": label, "n_before": before, "n_after": after, "n_removed": before - after})
    return pd.DataFrame(rows)


def apply_postfit_control(
    selected: pd.DataFrame,
    postfit: pd.DataFrame,
    viability: pd.DataFrame,
    grazing_limit: float,
) -> pd.DataFrame:
    post_cols = {
        "kepoi_name",
        "gilbert_qc_available",
        "nested_grazing_fraction",
        "planet_radius_fractional_uncertainty_approx",
    }
    require_columns(postfit, post_cols, "post-fit QC audit")
    postfit = postfit[list(post_cols)].drop_duplicates("kepoi_name")
    viability_cols = {"kepoi_name", "intended_guard_raises", "status"}
    require_columns(viability, viability_cols, "viability audit")
    viability = viability[list(viability_cols)].drop_duplicates("kepoi_name")

    merged = selected.merge(postfit, on="kepoi_name", how="left", validate="one_to_one")
    merged = merged.merge(viability, on="kepoi_name", how="left", validate="one_to_one")
    density_asymmetry = np.abs(merged["rhostar_err1"].abs() - merged["rhostar_err2"].abs()) / (
        0.5 * (merged["rhostar_err1"].abs() + merged["rhostar_err2"].abs())
    )
    merged["gilbert_density_asymmetry"] = density_asymmetry
    merged["exclude_qc_unavailable"] = ~merged["gilbert_qc_available"].fillna(False).astype(bool)
    merged["exclude_grazing"] = merged["nested_grazing_fraction"].gt(grazing_limit).fillna(True)
    merged["exclude_planet_radius_precision"] = (
        merged["planet_radius_fractional_uncertainty_approx"].gt(0.2).fillna(True)
    )
    merged["exclude_density_asymmetry"] = density_asymmetry.gt(0.30).fillna(True)
    # Gilbert intended to reject cases where fewer than 5% of importance
    # weights were numerically viable.  Use the source-audit interpretation,
    # while retaining missing viability rows as an explicit exclusion.
    merged["exclude_importance_viability"] = (
        merged["intended_guard_raises"].fillna(True).astype(bool)
        | merged["status"].fillna("missing").ne("ok")
    )
    exclusion_columns = [c for c in merged.columns if c.startswith("exclude_")]
    merged["gilbert_control_exclude"] = merged[exclusion_columns].any(axis=1)
    merged["gilbert_control_reasons"] = merged.apply(
        lambda row: ";".join(c.removeprefix("exclude_") for c in exclusion_columns if bool(row[c])),
        axis=1,
    )
    return merged


def fit_groups(selected: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sigmas = np.linspace(1e-4, 1.0, 2000)
    groups = {
        "all_small": selected,
        "observed_singles": selected[selected["npl"] == 1],
        "observed_multis": selected[selected["npl"] > 1],
        "thick_singles": selected[(selected["disk"] == "thick") & (selected["system"] == "single")],
        "thin_singles": selected[(selected["disk"] == "thin") & (selected["system"] == "single")],
        "thick_multis": selected[(selected["disk"] == "thick") & (selected["system"] == "multi")],
        "thin_multis": selected[(selected["disk"] == "thin") & (selected["system"] == "multi")],
    }
    rows: list[dict[str, object]] = []
    shape_rows: list[dict[str, object]] = []
    for group_name, group in groups.items():
        if len(group) < 2:
            continue
        files = group["posterior_file"].tolist()
        for mode in ("manuscript_reciprocal", "none", "legacy_forward_norm"):
            fit = fit_rayleigh(
                files,
                sigmas,
                apply_transit_selection=mode != "none",
                selection_mode=mode,
            )
            rows.append({"group": group_name, "n_planets": len(group), **fit})
        masses, e_grid = load_population_masses(
            group,
            apply_transit_selection=True,
            selection_mode="manuscript_reciprocal",
        )
        for model in ("beta", "monotonic_beta", "half_gaussian"):
            shape_fit = fit_model(masses, e_grid, model, starts_for(model))
            shape_rows.append(
                {
                    "group": group_name,
                    "n_planets": len(group),
                    "model": model,
                    "selection_mode": "manuscript_reciprocal",
                    **shape_fit,
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(shape_rows)


def main() -> None:
    args = parse_args()
    for path in (args.gilbert_catalog, args.summary, args.postfit_qc, args.viability_audit):
        if not path.exists():
            raise FileNotFoundError(path)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    catalog = pd.read_csv(args.gilbert_catalog, index_col=0)
    summary = pd.read_csv(args.summary)
    postfit = pd.read_csv(args.postfit_qc)
    viability = pd.read_csv(args.viability_audit)
    require_columns(
        catalog,
        {
            "planet_name", "koi_id", "npl", "period", "rp", "fpp", "disposition",
            "logg", "rcf", "ruwe", "rstar", "rstar_err1", "rstar_err2",
            "teff", "feh", "age", "rhostar_err1", "rhostar_err2",
        },
        "Gilbert catalog",
    )
    require_columns(summary, {"kepoi_name", "posterior_file"}, "posterior summary")
    if summary["kepoi_name"].duplicated().any():
        raise ValueError("Posterior summary contains duplicate planet names")
    missing_npz = [p for p in summary["posterior_file"] if not Path(p).exists()]
    if missing_npz:
        raise FileNotFoundError(f"{len(missing_npz)} posterior files are missing; first: {missing_npz[0]}")

    recovered_names = set(summary["kepoi_name"].astype(str))
    attr = attrition(catalog, recovered_names)
    prefit = catalog.loc[source_prefit_mask(catalog, small_planets=True)].copy()
    prefit = prefit.merge(
        summary,
        left_on="planet_name",
        right_on="kepoi_name",
        how="inner",
        validate="one_to_one",
        suffixes=("_gilbert", "_sagear"),
    )
    audited = apply_postfit_control(prefit, postfit, viability, args.grazing_limit)
    final = audited.loc[~audited["gilbert_control_exclude"]].copy()
    fits, shape_fits = fit_groups(final)

    attr.to_csv(args.output_dir / "attrition.csv", index=False)
    audited.to_csv(args.output_dir / "planet_qc_ledger.csv", index=False)
    final.to_csv(args.output_dir / "selected_planets.csv", index=False)
    fits.to_csv(args.output_dir / "rayleigh_control_fits.csv", index=False)
    shape_fits.to_csv(args.output_dir / "shape_control_fits.csv", index=False)

    report = {
        "control_definition": "Gilbert released catalog, small planets 0.5-3.5 R_earth, intersected with final Sagear ALDERAAN inventory",
        "gilbert_catalog_rows": int(len(catalog)),
        "gilbert_small_prefit_rows": int(source_prefit_mask(catalog, small_planets=True).sum()),
        "sagear_intersection_prefit_rows": int(len(prefit)),
        "sagear_intersection_postfit_rows": int(len(final)),
        "postfit_exclusions": {
            c: int(audited[c].sum()) for c in audited.columns if c.startswith("exclude_")
        },
        "primary_selection_mode": "manuscript_reciprocal",
        "source_basis": "Gilbert utils/models.py subtracts log(detection_prior) per posterior sample",
        "published_comparators": {
            "Gilbert2025_small_planets_mean_e": 0.05,
            "Gilbert2026_small_planet_observed_singles_mean_e": 0.073,
            "Gilbert2026_small_planet_multis_qualitative": "lower than singles; exact values depend on multiplicity bin",
        },
        "scope_warning": "This is a real-posterior intersection control, not Gilbert's full disk-agnostic sample.",
    }
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(fits[["group", "n_planets", "selection_mode", "expected_e", "expected_e_lo", "expected_e_hi"]].to_string(index=False))
    print(shape_fits[["group", "n_planets", "model", "mean_e", "status"]].to_string(index=False))


if __name__ == "__main__":
    main()
