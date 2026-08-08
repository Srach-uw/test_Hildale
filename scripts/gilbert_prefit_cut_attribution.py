"""Measure which released Gilbert pre-fit cuts change the Sagear hierarchy."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from hierarchical_rayleigh import fit_rayleigh, load_population_masses
from hierarchical_table3_order_diagnostic import fit_model, starts_for


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUMMARY = (
    ROOT
    / "tmp"
    / "exact_inventory_berger2020_fixed_dynesty_pairedperiod_rms_periastron_sir1000_20260807"
    / "exact_inventory_summary_density_draw.csv"
)
DEFAULT_GILBERT = (
    ROOT / "tmp" / "gjgilbert_kepler_ecc_rp" / "Catalogs" / "kepler_dr25_gaia_dr2_crossmatch.csv"
)


def fit_all(stage: str, frame: pd.DataFrame) -> list[dict]:
    files = frame.posterior_file.tolist()
    rayleigh = fit_rayleigh(
        files,
        np.linspace(1e-4, 1.0, 2000),
        apply_transit_selection=True,
        selection_mode="manuscript_reciprocal",
    )
    masses, e_grid = load_population_masses(
        frame, apply_transit_selection=True, selection_mode="manuscript_reciprocal"
    )
    beta = fit_model(masses, e_grid, "beta", starts_for("beta"))
    return [
        {
            "stage": stage,
            "model": "rayleigh",
            "n_planets": len(frame),
            "estimate": rayleigh["expected_e"],
            "lower": rayleigh["expected_e_lo"],
            "upper": rayleigh["expected_e_hi"],
        },
        {
            "stage": stage,
            "model": "beta_map",
            "n_planets": len(frame),
            "estimate": beta["mean_e"],
            "lower": np.nan,
            "upper": np.nan,
        },
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--gilbert-catalog", type=Path, default=DEFAULT_GILBERT)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "outputs" / "gilbert_prefit_cut_attribution"
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(args.summary)
    gilbert = pd.read_csv(args.gilbert_catalog, index_col=0)
    if summary.kepoi_name.duplicated().any() or gilbert.planet_name.duplicated().any():
        raise ValueError("Input planet identifiers must be unique")
    joined = summary.merge(
        gilbert,
        left_on="kepoi_name",
        right_on="planet_name",
        how="inner",
        validate="one_to_one",
        suffixes=("_sagear", "_gilbert"),
    )
    stellar_radius_error = np.hypot(joined.rstar_err1, joined.rstar_err2) / np.sqrt(2.0)
    criteria = [
        ("planet_radius_0p5_3p5", joined.rp.between(0.5, 3.5)),
        (
            "fpp_or_confirmed",
            (joined.fpp < 0.1) | joined.disposition.eq("CONFIRMED"),
        ),
        ("logg_gt_4", joined.logg > 4.0),
        ("furlan_lt_5pct_or_missing", ((joined.rcf - 1.0) < 0.05) | joined.rcf.isna()),
        ("ruwe_lt_1p4", joined.ruwe < 1.4),
        ("stellar_radius_precision_lt_20pct", (stellar_radius_error / joined.rstar) < 0.2),
        ("stellar_radius_0p7_1p4", joined.rstar.between(0.7, 1.4)),
        ("teff_4700_6500", joined.teff.between(4700.0, 6500.0)),
        ("feh_m0p6_p0p6", joined.feh.between(-0.6, 0.6)),
        ("age_0_14", joined.age.between(0.0, 14.0)),
    ]
    rows: list[dict] = []
    ledger_rows: list[dict] = []
    keep = pd.Series(True, index=joined.index)
    rows.extend(fit_all("matched_gilbert_catalog", joined))
    ledger_rows.append({"stage": "matched_gilbert_catalog", "n_planets": len(joined)})
    for name, criterion in criteria:
        keep &= criterion.fillna(False)
        stage = f"cumulative_through_{name}"
        frame = joined[keep]
        rows.extend(fit_all(stage, frame))
        ledger_rows.append({"stage": stage, "n_planets": len(frame)})
    fits = pd.DataFrame(rows)
    counts = pd.DataFrame(ledger_rows)
    fits.to_csv(args.output_dir / "hierarchical_cut_ladder.csv", index=False)
    counts.to_csv(args.output_dir / "attrition.csv", index=False)
    report = [
        "# Gilbert pre-fit cut attribution",
        "",
        "This diagnostic applies the released Gilbert cuts cumulatively to the recovered "
        "Sagear posterior inventory. Cut order follows the released source workflow.",
        "",
        fits.to_markdown(index=False),
    ]
    (args.output_dir / "report.md").write_text("\n".join(report), encoding="utf-8")
    print(fits.to_string(index=False))


if __name__ == "__main__":
    main()
