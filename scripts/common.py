from __future__ import annotations

import gzip
import json
import math
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
    df = pd.read_csv(path)
    required = {"kepid", "vx", "vy", "vz"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Angus table missing required columns: {sorted(missing)}")
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


def _float(text: str) -> float:
    text = text.strip()
    return float(text) if text else np.nan


def koi_target(kepoi_name: Any) -> str | None:
    if not isinstance(kepoi_name, str):
        return None
    return kepoi_name.split(".")[0]


def add_target_and_system(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["koi_target"] = out["kepoi_name"].map(koi_target)
    host_counts = out.groupby("kepid")["kepoi_name"].transform("count")
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
        chord = np.sqrt(np.maximum((1.0 + ror) ** 2 - impact**2, 0.0))
        arg = np.clip(chord / a_over_r, -1.0, 1.0)
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

