from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import output_dir


VARIANTS = {
    "current_no_dwarf_cut": lambda df: pd.Series(True, index=df.index),
    "berger_logg_ge_4p0": lambda df: pd.to_numeric(df["berger_logg"], errors="coerce") >= 4.0,
    "berger_rad_le_2": lambda df: pd.to_numeric(df["berger_rad"], errors="coerce") <= 2.0,
    "rho_linear_ge_0p05": lambda df: np.power(10.0, pd.to_numeric(df["rho_log"], errors="coerce")) >= 0.05,
    "berger_logg_ge_4p0_and_rad_le_2": lambda df: (
        (pd.to_numeric(df["berger_logg"], errors="coerce") >= 4.0)
        & (pd.to_numeric(df["berger_rad"], errors="coerce") <= 2.0)
    ),
    "berger_logg_ge_4p0_and_rho_ge_0p05": lambda df: (
        (pd.to_numeric(df["berger_logg"], errors="coerce") >= 4.0)
        & (np.power(10.0, pd.to_numeric(df["rho_log"], errors="coerce")) >= 0.05)
    ),
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Quantify evolved-star/dwarf-cut sensitivity for the Sagear sample.")
    parser.add_argument("--sample", default=None)
    parser.add_argument("--summary", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--filtered-variant", default="berger_logg_ge_4p0_and_rad_le_2", choices=sorted(VARIANTS))
    parser.add_argument("--filtered-summary-out", default=None)
    parser.add_argument("--filtered-sample-out", default=None)
    args = parser.parse_args()

    out_dir = output_dir()
    sample_path = Path(args.sample) if args.sample else out_dir / "canonical_sample_old_astropy_rawcc.csv"
    summary_path = Path(args.summary) if args.summary else out_dir / "eccentricity_posterior_summary_merged_paired_exact.csv"
    out_path = Path(args.out) if args.out else out_dir / "dwarf_cut_diagnostics.csv"
    filtered_summary_path = (
        Path(args.filtered_summary_out)
        if args.filtered_summary_out
        else out_dir / f"eccentricity_posterior_summary_merged_paired_exact_{args.filtered_variant}.csv"
    )
    filtered_sample_path = (
        Path(args.filtered_sample_out)
        if args.filtered_sample_out
        else out_dir / f"canonical_sample_old_astropy_rawcc_{args.filtered_variant}.csv"
    )

    sample = pd.read_csv(sample_path)
    summary = pd.read_csv(summary_path)
    summary_ids = set(summary["kepoi_name"].astype(str))

    rows = []
    for name, mask_func in VARIANTS.items():
        mask = mask_func(sample).fillna(False)
        sub = sample[mask].copy()
        removed = sample[~mask].copy()
        rows.append(summarize_variant(name, sub, removed, summary_ids))

    diagnostics = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    diagnostics.to_csv(out_path, index=False)

    keep_ids = set(sample[VARIANTS[args.filtered_variant](sample).fillna(False)]["kepoi_name"].astype(str))
    filtered_summary = summary[summary["kepoi_name"].astype(str).isin(keep_ids)].reset_index(drop=True)
    filtered_sample = sample[sample["kepoi_name"].astype(str).isin(keep_ids)].reset_index(drop=True)
    filtered_summary.to_csv(filtered_summary_path, index=False)
    filtered_sample.to_csv(filtered_sample_path, index=False)

    print(f"Wrote: {out_path}")
    print(f"Wrote: {filtered_summary_path}")
    print(f"Wrote: {filtered_sample_path}")
    print("\nDwarf/evolved-star cut diagnostics:")
    print(diagnostics.to_string(index=False))


def summarize_variant(name: str, sub: pd.DataFrame, removed: pd.DataFrame, summary_ids: set[str]) -> dict:
    row: dict[str, object] = {
        "variant": name,
        "sample_planets": int(len(sub)),
        "sample_hosts": int(sub["kepid"].nunique()),
        "removed_planets": int(len(removed)),
        "removed_hosts": int(removed["kepid"].nunique()),
        "posterior_planets": int(sub["kepoi_name"].astype(str).isin(summary_ids).sum()),
        "missing_posterior_planets": int((~sub["kepoi_name"].astype(str).isin(summary_ids)).sum()),
        "removed_median_rho_linear": median_rho(removed),
        "removed_median_berger_logg": float(np.nanmedian(pd.to_numeric(removed.get("berger_logg"), errors="coerce")))
        if len(removed)
        else np.nan,
        "removed_median_berger_rad": float(np.nanmedian(pd.to_numeric(removed.get("berger_rad"), errors="coerce")))
        if len(removed)
        else np.nan,
    }
    for disk in ["thin", "thick"]:
        for system in ["single", "multi"]:
            key = f"{disk}_{system}_planets"
            row[key] = int(((sub["disk"] == disk) & (sub["system"] == system)).sum())
    return row


def median_rho(df: pd.DataFrame) -> float:
    if len(df) == 0 or "rho_log" not in df:
        return np.nan
    rho = np.power(10.0, pd.to_numeric(df["rho_log"], errors="coerce"))
    return float(np.nanmedian(rho))


if __name__ == "__main__":
    main()
