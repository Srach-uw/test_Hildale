"""Apply reproducible Gilbert-style post-fit controls to the Sagear sample.

Sagear reports visual inspection but does not release the final planet-level
ledger.  This branch therefore tests whether objective transit-fit quality
criteria can account for the population discrepancy; it is not canonical
sample membership.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from sagear_residual_attribution import fit_stage


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUMMARY = (
    ROOT
    / "tmp"
    / "exact_inventory_berger2020_fixed_dynesty_pairedperiod_rms_periastron_sir1000_20260807"
    / "exact_inventory_summary_density_draw.csv"
)
DEFAULT_POSTFIT = ROOT / "outputs" / "gilbert_postfit_qc_audit_SOURCE_FAITHFUL_RMS_PERIASTRON_20260807.csv"
DEFAULT_VIABILITY = ROOT / "outputs" / "gilbert_importance_viability_SOURCE_FAITHFUL_RMS_PERIASTRON_20K_20260807.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--postfit", type=Path, default=DEFAULT_POSTFIT)
    parser.add_argument("--viability", type=Path, default=DEFAULT_VIABILITY)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "outputs" / "sagear_full_postfit_qc"
    )
    return parser.parse_args()


def build_qc_ledger(
    summary: pd.DataFrame, postfit: pd.DataFrame, viability: pd.DataFrame
) -> pd.DataFrame:
    if summary.kepoi_name.duplicated().any():
        raise ValueError("Summary contains duplicate planets")
    post_columns = [
        "kepoi_name",
        "gilbert_qc_available",
        "nested_grazing_fraction",
        "planet_radius_fractional_uncertainty_approx",
    ]
    viability_columns = ["kepoi_name", "intended_guard_raises", "status"]
    merged = summary.merge(
        postfit[post_columns].drop_duplicates("kepoi_name"),
        on="kepoi_name",
        how="left",
        validate="one_to_one",
    ).merge(
        viability[viability_columns].drop_duplicates("kepoi_name"),
        on="kepoi_name",
        how="left",
        validate="one_to_one",
    )
    denominator = 0.5 * (merged.rho_err_hi_solar.abs() + merged.rho_err_lo_solar.abs())
    merged["density_error_asymmetry"] = (
        (merged.rho_err_hi_solar.abs() - merged.rho_err_lo_solar.abs()).abs()
        / denominator.replace(0.0, np.nan)
    )
    merged["exclude_qc_unavailable"] = ~merged.gilbert_qc_available.fillna(False).astype(bool)
    merged["exclude_grazing"] = merged.nested_grazing_fraction.gt(0.05).fillna(True)
    merged["exclude_radius_precision"] = (
        merged.planet_radius_fractional_uncertainty_approx.gt(0.20).fillna(True)
    )
    merged["exclude_density_asymmetry"] = merged.density_error_asymmetry.gt(0.30).fillna(True)
    merged["exclude_importance_viability"] = (
        merged.intended_guard_raises.fillna(True).astype(bool)
        | merged.status.fillna("missing").ne("ok")
    )
    flags = [column for column in merged if column.startswith("exclude_")]
    merged["full_postfit_exclude"] = merged[flags].any(axis=1)
    merged["full_postfit_reasons"] = merged.apply(
        lambda row: ";".join(flag.removeprefix("exclude_") for flag in flags if row[flag]),
        axis=1,
    )
    return merged


def main() -> None:
    args = parse_args()
    for path in (args.summary, args.postfit, args.viability):
        if not path.exists():
            raise FileNotFoundError(path)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ledger = build_qc_ledger(
        pd.read_csv(args.summary), pd.read_csv(args.postfit), pd.read_csv(args.viability)
    )
    stages = {"recovered_all": ledger}
    keep = pd.Series(True, index=ledger.index)
    for flag in (
        "exclude_qc_unavailable",
        "exclude_grazing",
        "exclude_radius_precision",
        "exclude_density_asymmetry",
        "exclude_importance_viability",
    ):
        keep &= ~ledger[flag]
        stages[f"cumulative_through_{flag.removeprefix('exclude_')}"] = ledger[keep]

    fit_rows: list[dict] = []
    count_rows: list[dict] = []
    for stage, frame in stages.items():
        fit_rows.extend(fit_stage(stage, frame))
        for (disk, system), group in frame.groupby(["disk", "system"]):
            count_rows.append(
                {"stage": stage, "disk": disk, "system": system, "n_planets": len(group)}
            )
    fits = pd.DataFrame(fit_rows)
    counts = pd.DataFrame(count_rows)
    ledger.to_csv(args.output_dir / "planet_qc_ledger.csv", index=False)
    counts.to_csv(args.output_dir / "attrition_by_population.csv", index=False)
    fits.to_csv(args.output_dir / "hierarchical_sensitivity.csv", index=False)
    lines = [
        "# Full-sample post-fit QC sensitivity",
        "",
        "These Gilbert-style cuts are reproducible diagnostics. They do not reconstruct "
        "Sagear's unpublished visual-QC decisions.",
        "",
        "## Attrition",
        "",
        counts.to_markdown(index=False),
        "",
        "## Hierarchical results",
        "",
        fits.to_markdown(index=False),
    ]
    (args.output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print(counts.to_string(index=False))
    print(fits.to_string(index=False))


if __name__ == "__main__":
    main()
