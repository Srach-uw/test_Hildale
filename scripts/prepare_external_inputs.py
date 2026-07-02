from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from common import load_config, read_gaia_kepler, root_path


FURLAN_TABLE9_URL = "https://cdsarc.cds.unistra.fr/ftp/J/AJ/153/71/table9.dat"
APOGEE_ALLSTARLITE_URL = (
    "https://data.sdss.org/sas/dr17/apogee/spectro/aspcap/dr17/"
    "synspec_rev1/allStarLite-dr17-synspec_rev1.fits"
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Download/crossmatch external Sagear inputs.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--furlan", action="store_true", help="Download Furlan+2017 table9.dat.")
    parser.add_argument(
        "--apogee",
        action="store_true",
        help="Download APOGEE DR17 allStarLite and build Kepler/APOGEE chemical crossmatch.",
    )
    parser.add_argument("--all", action="store_true", help="Run all input-prep steps.")
    parser.add_argument("--max-sep-arcsec", type=float, default=1.0)
    args = parser.parse_args()

    cfg = load_config(args.config)
    if not (args.furlan or args.apogee or args.all):
        parser.error("Choose --furlan, --apogee, or --all.")

    if args.furlan or args.all:
        download_furlan(cfg)
    if args.apogee or args.all:
        prepare_apogee_crossmatch(cfg, args.max_sep_arcsec)


def download_furlan(cfg: dict) -> Path:
    out = root_path(cfg, "furlan_contamination")
    if out is None:
        out = Path(cfg["_root"]) / "data" / "furlan2017_table9.dat"
    download(FURLAN_TABLE9_URL, out)
    return out


def prepare_apogee_crossmatch(cfg: dict, max_sep_arcsec: float) -> Path:
    allstar_path = Path(cfg["_root"]) / "data" / "allStarLite-dr17-synspec_rev1.fits"
    download(APOGEE_ALLSTARLITE_URL, allstar_path)

    out = root_path(cfg, "apogee_dr17_chemical")
    if out is None:
        out = Path(cfg["_root"]) / "data" / "apogee_dr17_kepler_crossmatch.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    apogee = read_apogee_lite(allstar_path)
    kepler_cols = [c for c in ["kepid", "source_id", "ra", "dec"] if c in read_gaia_kepler(cfg).columns]
    kepler = read_gaia_kepler(cfg)[kepler_cols].dropna(subset=["kepid", "ra", "dec"]).drop_duplicates("kepid")

    if "source_id" in kepler.columns and "gaiaedr3_source_id" in apogee.columns:
        matched = match_by_source_id(kepler, apogee)
    else:
        matched = match_by_sky_position(kepler, apogee, max_sep_arcsec)

    matched = matched.sort_values("apogee_sep_arcsec").drop_duplicates("kepid")
    matched.to_csv(out, index=False)
    print(f"Wrote APOGEE crossmatch: {out} ({len(matched)} Kepler stars)")
    return out


def match_by_source_id(kepler: pd.DataFrame, apogee: pd.DataFrame) -> pd.DataFrame:
    kp = kepler.copy()
    ap = apogee.copy()
    kp["gaiaedr3_source_id"] = pd.to_numeric(kp["source_id"], errors="coerce").astype("Int64").astype(str)
    ap["gaiaedr3_source_id"] = pd.to_numeric(ap["gaiaedr3_source_id"], errors="coerce").astype("Int64").astype(str)
    ap = ap[ap["gaiaedr3_source_id"].notna() & (ap["gaiaedr3_source_id"] != "<NA>")]
    ap = ap.sort_values(["aspcapflag", "starflag", "mg_fe_flag", "fe_h_flag"], na_position="last")
    merged = kp.merge(ap, on="gaiaedr3_source_id", how="inner", suffixes=("", "_apogee"))
    merged["apogee_sep_arcsec"] = angular_sep_arcsec(
        merged["ra"], merged["dec"], merged["apogee_ra"], merged["apogee_dec"]
    )
    merged["apogee_match_method"] = "gaiaedr3_source_id"
    return merged


def match_by_sky_position(kepler: pd.DataFrame, apogee: pd.DataFrame, max_sep_arcsec: float) -> pd.DataFrame:
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    ap_coord = SkyCoord(
        ra=apogee["apogee_ra"].to_numpy(float) * u.deg,
        dec=apogee["apogee_dec"].to_numpy(float) * u.deg,
    )
    kp_coord = SkyCoord(
        ra=kepler["ra"].to_numpy(float) * u.deg,
        dec=kepler["dec"].to_numpy(float) * u.deg,
    )
    idx, sep2d, _ = kp_coord.match_to_catalog_sky(ap_coord)
    sep_arcsec = sep2d.to_value(u.arcsec)
    keep = sep_arcsec <= max_sep_arcsec
    matched = kepler.loc[keep].copy().reset_index(drop=True)
    ap_match = apogee.iloc[idx[keep]].reset_index(drop=True)
    matched["apogee_sep_arcsec"] = sep_arcsec[keep]
    for col in ap_match.columns:
        matched[col] = ap_match[col].to_numpy()
    matched["apogee_match_method"] = "sky_position"
    return matched


def angular_sep_arcsec(ra1, dec1, ra2, dec2) -> np.ndarray:
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    c1 = SkyCoord(ra=np.asarray(ra1, dtype=float) * u.deg, dec=np.asarray(dec1, dtype=float) * u.deg)
    c2 = SkyCoord(ra=np.asarray(ra2, dtype=float) * u.deg, dec=np.asarray(dec2, dtype=float) * u.deg)
    return c1.separation(c2).to_value(u.arcsec)


def read_apogee_lite(path: Path) -> pd.DataFrame:
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        names = list(data.names)
        ra_col = find_col(names, ["RA"])
        dec_col = find_col(names, ["DEC"])
        id_col = find_col(names, ["APOGEE_ID", "APSTAR_ID"])
        source_col = find_col(names, ["GAIAEDR3_SOURCE_ID"])
        feh_col = find_col(names, ["FE_H", "M_H"])
        mgfe_col = find_col(names, ["MG_FE"])
        teff_col = find_col(names, ["TEFF"])
        logg_col = find_col(names, ["LOGG"])
        aspcapflag_col = find_col(names, ["ASPCAPFLAG"])
        starflag_col = find_col(names, ["STARFLAG"])
        fehflag_col = find_col(names, ["FE_H_FLAG"])
        mgfeflag_col = find_col(names, ["MG_FE_FLAG"])
        out = pd.DataFrame(
            {
                "apogee_id": decode_fits_strings(data[id_col]),
                "gaiaedr3_source_id": decode_fits_strings(data[source_col]),
                "apogee_ra": np.asarray(data[ra_col], dtype=float),
                "apogee_dec": np.asarray(data[dec_col], dtype=float),
                "feh": np.asarray(data[feh_col], dtype=float),
                "mgfe": np.asarray(data[mgfe_col], dtype=float),
                "apogee_teff": np.asarray(data[teff_col], dtype=float),
                "apogee_logg": np.asarray(data[logg_col], dtype=float),
                "aspcapflag": np.asarray(data[aspcapflag_col]).astype(np.int64),
                "starflag": np.asarray(data[starflag_col]).astype(np.int64),
                "fe_h_flag": np.asarray(data[fehflag_col]).astype(np.int64),
                "mg_fe_flag": np.asarray(data[mgfeflag_col]).astype(np.int64),
            }
        )
    finite = np.isfinite(out[["apogee_ra", "apogee_dec", "feh", "mgfe"]]).all(axis=1)
    return out[finite].copy()


def find_col(names: list[str], options: list[str]) -> str:
    by_upper = {name.upper(): name for name in names}
    for option in options:
        if option.upper() in by_upper:
            return by_upper[option.upper()]
    raise KeyError(f"Could not find any of {options} in APOGEE columns.")


def decode_fits_strings(values) -> np.ndarray:
    arr = np.asarray(values)
    if arr.dtype.kind == "S":
        return np.char.decode(arr, "utf-8")
    return arr.astype(str)


def download(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        print(f"Already exists, skipping download: {path}")
        return
    print(f"Downloading {url}")
    print(f" -> {path}")
    urllib.request.urlretrieve(url, path)


if __name__ == "__main__":
    main()
