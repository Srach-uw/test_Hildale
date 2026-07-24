"""Audit ALDERAAN transit-shape posteriors against the Kepler catalog contract.

This is deliberately upstream of eccentricity inference.  It asks whether the
fitted period, duration, impact parameter, and radius ratio describe the same
transit as the catalog, and whether those shape values imply a plausible
circular stellar density relative to Berger et al. (2020).
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from extract_eccentricity_posteriors_direct import macdougall_rho_star_samp


def read_berger(path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            try:
                rows.append(
                    {
                        "kepid": int(line[0:8]),
                        "rho_log": float(line[116:122]),
                        "rho_log_upper": float(line[123:129]),
                        "rho_log_lower": float(line[130:136]),
                    }
                )
            except (ValueError, IndexError):
                continue
    return pd.DataFrame(rows)


def attach_catalog(fits: pd.DataFrame, catalog: pd.DataFrame) -> pd.DataFrame:
    catalog = catalog.copy()
    for col in ("kepid", "koi_period"):
        catalog[col] = pd.to_numeric(catalog[col], errors="coerce")
    catalog = catalog.dropna(subset=["kepid", "koi_period"])
    values = ["koi_duration", "koi_impact", "koi_ror", "koi_dor", "koi_count"]
    if "koi_model_snr" in catalog.columns:
        values.append("koi_model_snr")
    out = fits.copy()
    for col in values:
        out[col] = np.nan
    for index, row in out.iterrows():
        candidates = catalog[catalog["kepid"].eq(row["kepid"])]
        if candidates.empty:
            continue
        match = (candidates["koi_period"] - row["koi_period"]).abs().idxmin()
        for col in values:
            out.loc[index, col] = candidates.loc[match, col]
    out["catalog_match_period_days"] = out["koi_period"]
    out["catalog_period_relative_difference"] = (
        out["fit_period_days"] - out["koi_period"]
    ).abs() / out["koi_period"]
    return out


def audit(factorial: Path, catalog_path: Path, berger_path: Path) -> pd.DataFrame:
    data = pd.read_csv(factorial)
    data = data[data["arm"].isin(["original_lc", "reference_lc"])].copy()
    data = data.merge(read_berger(berger_path), on="kepid", how="left", validate="many_to_one")
    catalog = pd.read_csv(catalog_path, comment="#", low_memory=False)
    data = attach_catalog(data, catalog)

    period = data["fit_period_days"].to_numpy(float)
    duration_days = data["t14_hr_p50"].to_numpy(float) / 24.0
    impact = data["impact_p50"].to_numpy(float)
    ror = data["rp_over_rs_p50"].to_numpy(float)
    circular_rho = macdougall_rho_star_samp(
        period * 86400.0,
        duration_days * 86400.0,
        ror,
        impact,
        np.zeros(len(data)),
        np.zeros(len(data)),
    )
    data["rho_circular_from_fit"] = circular_rho
    catalog_rho = macdougall_rho_star_samp(
        data["koi_period"].to_numpy(float) * 86400.0,
        data["koi_duration"].to_numpy(float) * 3600.0 / 86400.0 * 86400.0,
        data["koi_ror"].to_numpy(float),
        data["koi_impact"].to_numpy(float),
        np.zeros(len(data)),
        np.zeros(len(data)),
    )
    data["rho_circular_from_catalog"] = catalog_rho
    data["rho_berger"] = np.power(10.0, data["rho_log"])
    data["rho_circular_to_berger"] = circular_rho / data["rho_berger"]
    data["rho_catalog_to_berger"] = catalog_rho / data["rho_berger"]
    data["rho_fit_to_catalog"] = circular_rho / catalog_rho
    data["duration_ratio_fit_to_catalog"] = data["t14_hr_p50"] / data["koi_duration"]
    data["impact_difference_fit_minus_catalog"] = data["impact_p50"] - data["koi_impact"]
    data["ror_ratio_fit_to_catalog"] = data["rp_over_rs_p50"] / data["koi_ror"]
    data["duration_flag"] = data["duration_ratio_fit_to_catalog"].sub(1.0).abs() > 0.15
    data["impact_flag"] = data["impact_difference_fit_minus_catalog"].abs() > 0.20
    data["ror_flag"] = data["ror_ratio_fit_to_catalog"].sub(1.0).abs() > 0.20
    data["density_flag"] = (
        data["rho_circular_to_berger"].lt(0.5) | data["rho_circular_to_berger"].gt(2.0)
    )
    data["contract_flag"] = data[["duration_flag", "impact_flag", "ror_flag", "density_flag"]].any(axis=1)
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--factorial", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--berger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.factorial, args.catalog, args.berger)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"rows={len(result)}")
    print("contract_flag_by_arm")
    print(result.groupby("arm")["contract_flag"].agg(["count", "sum", "mean"]).to_string())
    print("shape_summary_by_arm")
    print(
        result.groupby("arm")
        [[
            "duration_ratio_fit_to_catalog",
            "impact_difference_fit_minus_catalog",
            "ror_ratio_fit_to_catalog",
            "rho_circular_to_berger",
            "rho_catalog_to_berger",
            "rho_fit_to_catalog",
            "e50",
        ]]
        .median()
        .to_string()
    )
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
