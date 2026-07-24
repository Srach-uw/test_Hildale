"""Build a Berger-2018 KG-RADII logg/radius density-prior sample.

The official MAST KG-RADII stellar catalog releases logg and radius, but not
mass or density.  In solar units, rho/rho_sun = 10**(logg-logg_sun)/R.
The asymmetric envelope below propagates the released logg uncertainty and
radius uncertainties without importing Berger 2020 values.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table


LOGG_SUN = 4.438


def read_kg(path: Path) -> pd.DataFrame:
    with fits.open(path, memmap=False) as hdul:
        frame = Table(hdul[1].data).to_pandas()
    frame = frame.rename(columns={
        "KIC_ID": "kepid",
        "logg": "berger2018_kg_logg",
        "logg_err": "berger2018_kg_logg_err",
        "Radius": "berger2018_kg_radius",
        "Radius_err_upper": "berger2018_kg_radius_err_upper",
        "Radius_err_lower": "berger2018_kg_radius_err_lower",
    })
    return frame[[
        "kepid", "berger2018_kg_logg", "berger2018_kg_logg_err",
        "berger2018_kg_radius", "berger2018_kg_radius_err_upper",
        "berger2018_kg_radius_err_lower",
    ]].drop_duplicates("kepid")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", default="outputs/sagear2026_planet_inventory_pre_visual_qc.csv")
    parser.add_argument("--kg", default="inputs/berger2018_kg_radii_star_cat.fits")
    parser.add_argument("--out", default="outputs/sagear2026_sample_berger2018_kg_logg_radius.csv")
    args = parser.parse_args()
    sample = pd.read_csv(args.sample, low_memory=False)
    kg = read_kg(Path(args.kg))
    sample["kepid"] = pd.to_numeric(sample["kepid"], errors="coerce").astype("Int64")
    kg["kepid"] = pd.to_numeric(kg["kepid"], errors="coerce").astype("Int64")
    out = sample.drop(columns=[c for c in ["rho_log", "rho_log_upper", "rho_log_lower"] if c in sample]).merge(kg, on="kepid", how="left")

    logg = pd.to_numeric(out["berger2018_kg_logg"], errors="coerce")
    logg_err = pd.to_numeric(out["berger2018_kg_logg_err"], errors="coerce")
    radius = pd.to_numeric(out["berger2018_kg_radius"], errors="coerce")
    radius_hi = pd.to_numeric(out["berger2018_kg_radius_err_upper"], errors="coerce")
    radius_lo = pd.to_numeric(out["berger2018_kg_radius_err_lower"], errors="coerce")
    rho = 10.0 ** (logg - LOGG_SUN) / radius
    rho_hi = 10.0 ** (logg + logg_err - LOGG_SUN) / (radius - radius_lo)
    rho_lo = 10.0 ** (logg - logg_err - LOGG_SUN) / (radius + radius_hi)
    out["rho_log"] = np.log10(rho.where(rho > 0))
    out["rho_log_upper"] = np.log10(rho_hi.where(rho_hi > 0)) - out["rho_log"]
    out["rho_log_lower"] = out["rho_log"] - np.log10(rho_lo.where(rho_lo > 0))
    out["berger2018_density_source"] = "MAST_KG_RADII_logg_radius"
    out["berger2018_density_available"] = out[["rho_log", "rho_log_upper", "rho_log_lower"]].notna().all(axis=1)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"KG rows={len(kg)} hosts={sample['kepid'].nunique()} planets={len(sample)} available={int(out['berger2018_density_available'].sum())}")
    print(f"output={Path(args.out).resolve()}")
    print(out[["rho_log", "rho_log_upper", "rho_log_lower"]].describe().to_string())


if __name__ == "__main__":
    main()
