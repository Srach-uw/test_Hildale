from __future__ import annotations

import gzip
import json
import math
from io import StringIO
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


PIPELINE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = PIPELINE_DIR / "config.json"


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    cfg_path = Path(path) if path else DEFAULT_CONFIG
    with cfg_path.open("r", encoding="utf-8") as f:
        cfg = json.load(f)
    root = Path(cfg["paths"]["research_root"]).expanduser()
    cfg["_config_path"] = str(cfg_path)
    cfg["_root"] = str(root)
    return cfg


def root_path(cfg: dict[str, Any], key: str) -> Path | None:
    value = cfg["paths"].get(key)
    if value in (None, ""):
        return None
    p = Path(value)
    if not p.is_absolute():
        p = Path(cfg["_root"]) / p
    return p


def output_dir() -> Path:
    out = PIPELINE_DIR / "outputs"
    out.mkdir(parents=True, exist_ok=True)
    return out


def read_koi(cfg: dict[str, Any]) -> pd.DataFrame:
    path = root_path(cfg, "koi_catalog")
    if path is None or not path.exists():
        raise FileNotFoundError(f"KOI catalog not found: {path}")
    return pd.read_csv(path, comment="#")


def read_gaia_kepler(cfg: dict[str, Any]) -> pd.DataFrame:
    from astropy.table import Table

    path = root_path(cfg, "gaia_kepler_fits")
    if path is None or not path.exists():
        raise FileNotFoundError(f"Gaia/Kepler FITS not found: {path}")
    return Table.read(path, format="fits").to_pandas()


def read_angus(cfg: dict[str, Any]) -> pd.DataFrame:
    path = root_path(cfg, "angus_velocities")
    if path is None or not path.exists():
        raise FileNotFoundError(f"Angus velocity table not found: {path}")
    raw_path = path.with_name("angus_velocities.dat.gz")
    if raw_path.exists():
        rows = []
        with gzip.open(raw_path, "rt") as handle:
            for line in handle:
                try:
                    kepid = int(line[7:15])
                    ra = float(line[36:43])
                    dec = float(line[50:56])
                    plx = float(line[63:69])
                    dist_text = line[76:83].strip()
                    dist_pc = float(dist_text) if dist_text else np.nan
                    vx_calc = line[175:182].strip()
                    vy_calc = line[195:202].strip()
                    vz_calc = line[214:221].strip()
                    vx = float(vx_calc) if vx_calc else float(line[183:189])
                    vy = float(vy_calc) if vy_calc else float(line[203:208])
                    vz = float(vz_calc) if vz_calc else float(line[222:228])
                    rows.append((kepid, ra, dec, plx, dist_pc, vx, vy, vz, bool(vx_calc), "Angus2022_raw"))
                except (TypeError, ValueError):
                    continue
        raw = pd.DataFrame(
            rows,
            columns=[
                "kepid",
                "angus_ra",
                "angus_dec",
                "angus_parallax",
                "angus_dist_pc",
                "vx",
                "vy",
                "vz",
                "has_rv",
                "angus_source",
            ],
        )
        missing_distance = ~np.isfinite(raw["angus_dist_pc"])
        usable_parallax = missing_distance & np.isfinite(raw["angus_parallax"]) & (raw["angus_parallax"] > 0)
        raw.loc[usable_parallax, "angus_dist_pc"] = 1000.0 / raw.loc[usable_parallax, "angus_parallax"]
        if not raw.empty:
            return raw
    df = pd.read_csv(path)
    required = {"kepid", "vx", "vy", "vz"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Angus table missing required columns: {sorted(missing)}")
    df["angus_source"] = "parsed_csv_without_positions"
    return df


def read_berger_table2(cfg: dict[str, Any]) -> pd.DataFrame:
    path = root_path(cfg, "berger_table2")
    if path is None or not path.exists():
        raise FileNotFoundError(f"Berger table2 not found: {path}")
    rows: list[dict[str, Any]] = []
    with gzip.open(path, "rt") as f:
        for line in f:
            try:
                age_text = line[158:163].strip()
                age_flag = line[164:165].strip()
                rows.append(
                    {
                        "kepid": int(line[0:8].strip()),
                        "berger_mass": _float(line[9:14]),
                        "berger_teff": _float(line[28:35]),
                        "berger_logg": _float(line[51:57]),
                        "berger_feh": _float(line[71:77]),
                        "berger_rad": _float(line[91:98]),
                        "rho_log": _float(line[116:122]),
                        "rho_log_upper": _float(line[123:129]),
                        "rho_log_lower": _float(line[130:136]),
                        "berger_age": float(age_text)
                        if age_text and age_flag != "*"
                        else np.nan,
                    }
                )
            except Exception:
                continue
    return pd.DataFrame(rows)


def read_berger2018_stellar_table(path: str | Path) -> pd.DataFrame:
    """Read the VizieR TSV export for Berger+2018 radii/evolution flags."""
    path = Path(path)
    lines = path.read_text(errors="replace").splitlines()
    header_idx = next((i for i, line in enumerate(lines) if line.startswith("KIC\t")), None)
    if header_idx is None:
        raise ValueError(f"Could not find Berger+2018 VizieR header in {path}")
    data = "\n".join([lines[header_idx]] + lines[header_idx + 3 :])
    tab = pd.read_csv(StringIO(data), sep="\t")
    out = tab.rename(
        columns={
            "KIC": "kepid",
            "Teff": "berger2018_teff",
            "R*": "berger2018_rad",
            "E_R*": "berger2018_rad_err_upper",
            "e_R*": "berger2018_rad_err_lower",
            "Evol": "berger2018_evol",
            "Bin": "berger2018_bin",
        }
    )
    out["kepid"] = pd.to_numeric(out["kepid"], errors="coerce")
    out = out.dropna(subset=["kepid"]).copy()
    out["kepid"] = out["kepid"].astype(int)
    for col in [
        "berger2018_teff",
        "berger2018_rad",
        "berger2018_rad_err_upper",
        "berger2018_rad_err_lower",
        "berger2018_evol",
        "berger2018_bin",
    ]:
        if col in out:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out.drop_duplicates("kepid")


def read_sagear2026_kinematic_hosts(cfg: dict[str, Any]) -> pd.DataFrame:
    """Read Sagear et al. (2026) AJ 172, 42 machine-readable Table 1."""
    path = root_path(cfg, "sagear2026_kinematic_hosts")
    if path is None or not path.exists():
        raise FileNotFoundError(f"Published Sagear host table not found: {path}")
    names = [
        "kepid", "gaia_dr3_source_id", "ra_deg", "dec_deg", "parallax_mas",
        "z_galcen_kpc", "r_galcen_kpc", "phi_galcen_rad", "vr_infer_kms",
        "vphi_infer_kms", "vz_infer_kms", "vr_calc_kms", "vphi_calc_kms",
        "vz_calc_kms", "p_thick_published", "disk_published",
    ]
    colspecs = [
        (0, 8), (9, 28), (29, 37), (38, 45), (46, 53), (54, 60), (61, 67),
        (68, 74), (75, 84), (85, 94), (95, 100), (101, 110), (111, 120),
        (121, 128), (129, 135), (136, 141),
    ]
    table = pd.read_fwf(path, colspecs=colspecs, names=names, skiprows=28)
    table["kepid"] = pd.to_numeric(table["kepid"], errors="coerce")
    table = table.dropna(subset=["kepid"]).copy()
    table["kepid"] = table["kepid"].astype(int)
    table["gaia_dr3_source_id"] = table["gaia_dr3_source_id"].astype(str).str.strip()
    numeric = [c for c in names if c not in {"kepid", "gaia_dr3_source_id", "disk_published"}]
    for column in numeric:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table["disk_published"] = table["disk_published"].astype(str).str.strip().str.lower()
    table["has_measured_velocity"] = table[["vr_calc_kms", "vphi_calc_kms", "vz_calc_kms"]].ne(0).any(axis=1)
    table["vr_used_kms"] = table["vr_calc_kms"].where(table["has_measured_velocity"], table["vr_infer_kms"])
    table["vphi_used_kms"] = table["vphi_calc_kms"].where(table["has_measured_velocity"], table["vphi_infer_kms"])
    table["vz_used_kms"] = table["vz_calc_kms"].where(table["has_measured_velocity"], table["vz_infer_kms"])
    if table["kepid"].duplicated().any():
        raise ValueError("Published Sagear host table contains duplicate KIC identifiers")
    if not table["disk_published"].isin({"thin", "thick"}).all():
        raise ValueError("Published Sagear host table contains invalid disk labels")
    return table.reset_index(drop=True)


def _float(text: str) -> float:
    text = text.strip()
    return float(text) if text else np.nan


def koi_target(kepoi_name: Any) -> str | None:
    if not isinstance(kepoi_name, str):
        return None
    return kepoi_name.split(".")[0]


def add_target_and_system(
    df: pd.DataFrame,
    host_multiplicity: pd.Series | None = None,
) -> pd.DataFrame:
    """Attach KOI target and architecture labels.

    Pass a host-indexed multiplicity series derived before downstream quality
    cuts to preserve known system architecture. Without it, multiplicity is
    computed from the supplied frame and is therefore suitable only for
    explicitly requested filtered-sample sensitivities.
    """
    out = df.copy()
    out["koi_target"] = out["kepoi_name"].map(koi_target)
    if host_multiplicity is None:
        host_counts = out.groupby("kepid")["kepoi_name"].transform("count")
    else:
        if host_multiplicity.index.has_duplicates:
            raise ValueError("Host multiplicity index contains duplicate KIC identifiers")
        host_counts = out["kepid"].map(host_multiplicity)
        if host_counts.isna().any():
            missing = sorted(out.loc[host_counts.isna(), "kepid"].astype(int).unique())
            raise ValueError(f"Raw host multiplicity is missing KIC identifiers: {missing[:10]}")
    out["system"] = np.where(host_counts == 1, "single", "multi")
    return out


def add_audit_row(rows: list[dict[str, Any]], stage: str, df: pd.DataFrame, note: str = "") -> None:
    rows.append(
        {
            "stage": stage,
            "planets": int(len(df)),
            "hosts": int(df["kepid"].nunique()) if "kepid" in df else np.nan,
            "thin_single_planets": _count(df, disk="thin", system="single"),
            "thick_single_planets": _count(df, disk="thick", system="single"),
            "thin_multi_planets": _count(df, disk="thin", system="multi"),
            "thick_multi_planets": _count(df, disk="thick", system="multi"),
            "note": note,
        }
    )


def _count(df: pd.DataFrame, disk: str, system: str) -> int | float:
    if "disk" not in df or "system" not in df:
        return np.nan
    return int(((df["disk"] == disk) & (df["system"] == system)).sum())


def normalize_dynesty_weights(logwt: np.ndarray) -> np.ndarray:
    finite = np.isfinite(logwt)
    if not finite.any():
        return np.full(len(logwt), 1.0 / len(logwt))
    m = np.nanmax(logwt[finite])
    w = np.zeros(len(logwt), dtype=float)
    w[finite] = np.exp(logwt[finite] - m)
    total = w.sum()
    if total <= 0:
        return np.full(len(logwt), 1.0 / len(logwt))
    return w / total


def trapezoid(y: np.ndarray, x: np.ndarray | None = None, axis: int = -1) -> np.ndarray:
    """NumPy 2 removed np.trapz; keep one compatibility point for all scripts."""
    fn = getattr(np, "trapezoid", None)
    if fn is None:
        fn = np.trapz
    return fn(y, x=x, axis=axis)


def stellar_a_over_r(rho_log: np.ndarray, period_days: np.ndarray) -> np.ndarray:
    rho_sun_cgs = 1.408
    g_cgs = 6.674e-8
    rho_cgs = np.power(10.0, rho_log) * rho_sun_cgs
    p_sec = period_days * 86400.0
    return np.power(g_cgs * rho_cgs * p_sec * p_sec / (3.0 * np.pi), 1.0 / 3.0)


def circular_duration_days(
    period_days: np.ndarray,
    impact: np.ndarray,
    ror: np.ndarray,
    a_over_r: np.ndarray,
) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        numerator_sq = (1.0 + ror) ** 2 - impact**2
        denominator_sq = a_over_r**2 - impact**2
        valid = (numerator_sq > 0.0) & (denominator_sq > 0.0)
        arg = np.full(np.broadcast(period_days, impact, ror, a_over_r).shape, np.nan, dtype=float)
        arg[valid] = np.sqrt(numerator_sq[valid] / denominator_sq[valid])
        arg = np.clip(arg, -1.0, 1.0)
        return (period_days / np.pi) * np.arcsin(arg)


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True)


def finite_columns(df: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    mask = pd.Series(True, index=df.index)
    for col in columns:
        mask &= np.isfinite(pd.to_numeric(df[col], errors="coerce"))
    return mask
