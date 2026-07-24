"""Build explicit Berger-density sensitivity samples for the Sagear audit.

The published Sagear text names Berger et al. (2018) for stellar densities,
but the released Berger 2018 table contains radius, evolution, and binary
columns, not a density column.  This script therefore constructs diagnostic
hybrids rather than pretending that an exact 2018 density prior is available.

The two outputs are:

* ``b18hybrid_fractional``: Berger 2020 mass and Berger 2018 radius, with the
  fractional asymmetric density errors inherited from the published Berger
  2020 density;
* ``b18hybrid_radius_only``: the same central density, with errors propagated
  from the Berger 2018 radius alone while holding mass fixed.

Both retain the original inventory columns and can be passed directly to the
direct eccentricity extractor, whose density merge preserves existing values.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from berger_density_provenance_audit import parse_b18, parse_b20


def build(sample: pd.DataFrame, b20: pd.DataFrame, b18: pd.DataFrame) -> pd.DataFrame:
    out = sample.copy()
    out["kepid"] = pd.to_numeric(out["kepid"], errors="coerce").astype("Int64")
    b20 = b20.copy()
    b18 = b18.copy()
    b20["kepid"] = pd.to_numeric(b20["kepid"], errors="coerce").astype("Int64")
    b18["kepid"] = pd.to_numeric(b18["kepid"], errors="coerce").astype("Int64")
    lookup = b20.merge(b18, on="kepid", how="left", suffixes=("", "_b18"))
    out = out.drop(columns=[c for c in [
        "rho_log", "rho_log_upper", "rho_log_lower",
        "berger2018_rad", "berger2018_rad_err_upper", "berger2018_rad_err_lower",
    ] if c in out])
    out = out.merge(lookup, on="kepid", how="left", suffixes=("", "_lookup"))

    mass = pd.to_numeric(out["mass"], errors="coerce")
    radius = pd.to_numeric(out["radius_b18"], errors="coerce")
    radius_hi = pd.to_numeric(out["radius_b18_err_upper"], errors="coerce")
    radius_lo = pd.to_numeric(out["radius_b18_err_lower"], errors="coerce")
    rho = mass / radius.pow(3)
    out["rho_log_b18hybrid"] = np.log10(rho.where(rho > 0))

    b20_rho = 10.0 ** pd.to_numeric(out["rho_log"], errors="coerce")
    b20_hi = pd.to_numeric(out["rho_log_upper"], errors="coerce")
    b20_lo = pd.to_numeric(out["rho_log_lower"], errors="coerce")
    out["rho_log_upper_b18hybrid_fractional"] = b20_hi
    out["rho_log_lower_b18hybrid_fractional"] = b20_lo

    # Holding B20 mass fixed, the radius-only density envelope is exact for
    # the asymmetric radius errors supplied by Berger 2018.
    rho_hi = mass / (radius - radius_lo).pow(3)
    rho_lo = mass / (radius + radius_hi).pow(3)
    out["rho_log_upper_b18hybrid_radius_only"] = np.log10(rho_hi.where(rho_hi > 0)) - out["rho_log_b18hybrid"]
    out["rho_log_lower_b18hybrid_radius_only"] = out["rho_log_b18hybrid"] - np.log10(rho_lo.where(rho_lo > 0))

    required = ["rho_log_b18hybrid", "rho_log_upper_b18hybrid_fractional", "rho_log_lower_b18hybrid_fractional"]
    out["b18hybrid_available"] = out[required].notna().all(axis=1)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--sample", default="outputs/sagear2026_planet_inventory_pre_visual_qc.csv")
    parser.add_argument("--out-dir", default="outputs")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    b20 = parse_b20(root.parent / "table2.dat.gz")
    b18_path = root / "inputs" / "sagear_2026_published" / "berger2018_table1.dat.gz"
    if not b18_path.exists():
        b18_path = root.parent / "data" / "berger2018_table1_min.tsv"
    b18 = parse_b18(b18_path)
    sample = pd.read_csv(root / args.sample, low_memory=False)
    built = build(sample, b20, b18)

    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    common = [c for c in sample.columns if c in built.columns]
    base = built[common + [
        "rho_log_b18hybrid", "rho_log_upper_b18hybrid_fractional",
        "rho_log_lower_b18hybrid_fractional", "rho_log_upper_b18hybrid_radius_only",
        "rho_log_lower_b18hybrid_radius_only", "b18hybrid_available",
    ]]
    fractional = base.copy()
    fractional["rho_log"] = fractional["rho_log_b18hybrid"]
    fractional["rho_log_upper"] = fractional["rho_log_upper_b18hybrid_fractional"]
    fractional["rho_log_lower"] = fractional["rho_log_lower_b18hybrid_fractional"]
    radius_only = base.copy()
    radius_only["rho_log"] = radius_only["rho_log_b18hybrid"]
    radius_only["rho_log_upper"] = radius_only["rho_log_upper_b18hybrid_radius_only"]
    radius_only["rho_log_lower"] = radius_only["rho_log_lower_b18hybrid_radius_only"]
    fpath = out_dir / "sagear2026_sample_berger2018_hybrid_fractional.csv"
    rpath = out_dir / "sagear2026_sample_berger2018_hybrid_radius_only.csv"
    fractional.to_csv(fpath, index=False)
    radius_only.to_csv(rpath, index=False)
    print(f"hosts={sample['kepid'].nunique()} planets={len(sample)} b18hybrid_rows={int(base['b18hybrid_available'].sum())}")
    print(f"fractional={fpath}")
    print(f"radius_only={rpath}")


if __name__ == "__main__":
    main()
