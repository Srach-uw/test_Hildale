"""Prepare a provenance-preserving input for the standalone injection diagnostic.

This tool reads local Kepler PDCSAP FITS and one ALDERAAN result FITS. It does
not download data, detrend flux, subtract fitted transits, run ALDERAAN, or
change any canonical result. The output contains the three fields required by
``diagnostic_transit_injection.py`` plus cadence provenance needed to audit it.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
from pathlib import Path
import platform
import re
import sys
from typing import Any, Iterable, Sequence

import numpy as np
from astropy import __version__ as ASTROPY_VERSION
from astropy.io import fits


PINNED_LIGHTKURVE_VERSION = "2.1.1"
OUTPUT_NAME = "diagnostic_injection_input.npz"
MANIFEST_NAME = "diagnostic_injection_input_manifest.json"
TTIMES_PATTERN = re.compile(r"^TTIMES_(\d+)$")
TIMING_COLUMNS = ("INDEX", "TTIME", "MODEL", "OUT_PROB", "OUT_FLAG")
SYSTEM_CATALOG_COLUMNS = ("koi_id", "kic_id", "npl", "period", "epoch", "duration")
BKJD_REFERENCE = 2454833.0
EXPOSURE_MATCH_TOLERANCE_MINUTES = 0.05


class PreparationInputError(ValueError):
    """Raised when a local input cannot support an auditable preparation."""


class LightkurveDependencyError(RuntimeError):
    """Raised when the pinned ALDERAAN quality-mask implementation is absent."""


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of one immutable local source file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_array(values: np.ndarray) -> str:
    """Return a content hash after making array byte order and layout explicit."""
    array = np.ascontiguousarray(values)
    return hashlib.sha256(array.tobytes()).hexdigest()


def _output_directory(output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise PreparationInputError("output directory must be absent or empty")
    output_dir.mkdir(parents=True, exist_ok=True)


def _integer_header(header: fits.Header, names: Iterable[str], label: str) -> int:
    for name in names:
        raw = header.get(name)
        try:
            value = int(raw)
        except (TypeError, ValueError, OverflowError):
            continue
        if value > 0:
            return value
    raise PreparationInputError(f"primary header is missing a positive {label}")


def _quarter_header(primary: fits.Header, extension: fits.Header, path: Path) -> int:
    raw = _shared_header_value(primary, extension, "QUARTER", path)
    try:
        quarter = int(raw)
    except (TypeError, ValueError, OverflowError) as error:
        raise PreparationInputError(f"{path.name} has a non-integer QUARTER") from error
    if quarter < 0:
        raise PreparationInputError(f"{path.name} QUARTER must be non-negative")
    return quarter


def _shared_header_value(
    primary: fits.Header, extension: fits.Header, name: str, path: Path
) -> Any:
    primary_value = primary.get(name)
    extension_value = extension.get(name)
    if primary_value is not None and extension_value is not None and primary_value != extension_value:
        raise PreparationInputError(f"{path.name} has contradictory primary and LIGHTCURVE {name}")
    value = extension_value if extension_value is not None else primary_value
    if value is None:
        raise PreparationInputError(f"{path.name} is missing standard {name}")
    return value


def _shared_normalized_header_value(
    primary: fits.Header,
    extension: fits.Header,
    name: str,
    path: Path,
    normalizer: Any,
    equivalent: Any | None = None,
) -> Any:
    """Read one standard keyword while rejecting contradictory header values."""
    primary_raw = primary.get(name)
    extension_raw = extension.get(name)
    primary_value = normalizer(primary_raw, path) if primary_raw is not None else None
    extension_value = normalizer(extension_raw, path) if extension_raw is not None else None
    same = equivalent or (lambda left, right: left == right)
    if primary_value is not None and extension_value is not None and not same(primary_value, extension_value):
        raise PreparationInputError(f"{path.name} has contradictory primary and LIGHTCURVE {name}")
    value = extension_value if extension_value is not None else primary_value
    if value is None:
        raise PreparationInputError(f"{path.name} is missing standard {name}")
    return value


def _normalize_obsmode(raw: Any, path: Path) -> str:
    value = str(raw).strip().lower().replace("_", " ").replace("-", " ")
    value = " ".join(value.split())
    aliases = {
        "long": "long",
        "long cadence": "long",
        "lc": "long",
        "short": "short",
        "short cadence": "short",
        "sc": "short",
    }
    try:
        return aliases[value]
    except KeyError as error:
        raise PreparationInputError(f"{path.name} has unsupported OBSMODE {raw!r}") from error


def _normalize_timedelta(raw: Any, path: Path) -> float:
    try:
        value = float(raw)
    except (TypeError, ValueError, OverflowError) as error:
        raise PreparationInputError(f"{path.name} has an invalid TIMEDEL") from error
    if not np.isfinite(value) or value <= 0.0:
        raise PreparationInputError(f"{path.name} TIMEDEL must be finite and positive")
    return value


def _same_timedelta(left: float, right: float) -> bool:
    """Allow only floating-point representation noise between paired headers."""
    return bool(np.isclose(left, right, rtol=0.0, atol=1e-12))


def _cadence_headers(primary: fits.Header, extension: fits.Header, path: Path) -> dict[str, Any]:
    observed_mode = _shared_normalized_header_value(
        primary, extension, "OBSMODE", path, _normalize_obsmode
    )
    timedelta_days = _shared_normalized_header_value(
        primary, extension, "TIMEDEL", path, _normalize_timedelta, _same_timedelta
    )
    return {
        "observed_mode": observed_mode,
        "timedel_days": float(timedelta_days),
        "observed_exposure_minutes": float(timedelta_days * 1440.0),
    }


def _bjd_reference(header: fits.Header, path: Path) -> float | None:
    direct = header.get("BJDREF")
    integer = header.get("BJDREFI")
    fractional = header.get("BJDREFF")
    if direct is None and integer is None and fractional is None:
        return None
    try:
        split = None if integer is None else float(integer) + float(fractional or 0.0)
        value = float(direct) if direct is not None else split
    except (TypeError, ValueError, OverflowError) as error:
        raise PreparationInputError(f"{path.name} has an invalid BJDREF value") from error
    if value is None or not np.isfinite(value):
        raise PreparationInputError(f"{path.name} has a non-finite BJDREF value")
    if split is not None and direct is not None and not np.isclose(value, split, rtol=0.0, atol=1e-9):
        raise PreparationInputError(f"{path.name} has contradictory BJDREF and BJDREFI/BJDREFF")
    return value


def _kepler_time_convention(primary: fits.Header, extension: fits.Header, path: Path) -> dict[str, Any]:
    """Validate the Kepler BKJD convention instead of inferring a time basis."""
    raw_unit = str(_shared_header_value(primary, extension, "TIMEUNIT", path)).strip().lower()
    if raw_unit not in {"d", "day", "days"}:
        raise PreparationInputError(f"{path.name} has unsupported TIMEUNIT {raw_unit!r}; expected days")
    raw_system = str(_shared_header_value(primary, extension, "TIMESYS", path)).strip().upper()
    if raw_system != "TDB":
        raise PreparationInputError(f"{path.name} has unsupported TIMESYS {raw_system!r}; expected TDB")
    primary_reference = _bjd_reference(primary, path)
    extension_reference = _bjd_reference(extension, path)
    if primary_reference is None and extension_reference is None:
        raise PreparationInputError(f"{path.name} is missing standard BJDREF or BJDREFI/BJDREFF")
    if primary_reference is not None and extension_reference is not None and not np.isclose(
        primary_reference, extension_reference, rtol=0.0, atol=1e-9
    ):
        raise PreparationInputError(f"{path.name} has contradictory primary and LIGHTCURVE BJDREF")
    reference = extension_reference if extension_reference is not None else primary_reference
    assert reference is not None
    if not np.isclose(reference, BKJD_REFERENCE, rtol=0.0, atol=1e-9):
        raise PreparationInputError(
            f"{path.name} has unsupported BJDREF {reference}; expected Kepler BKJD reference {BKJD_REFERENCE}"
        )
    for header in (primary, extension):
        raw_zero = header.get("TIMEZERO", 0.0)
        try:
            zero = float(raw_zero)
        except (TypeError, ValueError, OverflowError) as error:
            raise PreparationInputError(f"{path.name} has an invalid TIMEZERO") from error
        if not np.isfinite(zero) or not np.isclose(zero, 0.0, rtol=0.0, atol=1e-12):
            raise PreparationInputError(f"{path.name} has unsupported nonzero TIMEZERO")
    return {
        "time_reference": "BKJD",
        "timeunit": "d",
        "timesys": raw_system,
        "bjdref": float(reference),
        "timezero": 0.0,
    }


def _load_pinned_lightkurve() -> Any:
    try:
        lightkurve = importlib.import_module("lightkurve")
    except ModuleNotFoundError as error:
        raise LightkurveDependencyError(
            "lightkurve 2.1.1 from the pinned ALDERAAN environment is required "
            "for the default Kepler quality mask"
        ) from error
    version = str(getattr(lightkurve, "__version__", "unknown"))
    if version != PINNED_LIGHTKURVE_VERSION:
        raise LightkurveDependencyError(
            f"expected pinned lightkurve {PINNED_LIGHTKURVE_VERSION}, found {version}"
        )
    try:
        lightkurve.KeplerQualityFlags.create_quality_mask
    except AttributeError as error:
        raise LightkurveDependencyError(
            "pinned lightkurve lacks KeplerQualityFlags.create_quality_mask"
        ) from error
    return lightkurve


def default_kepler_quality_mask(quality: np.ndarray) -> tuple[np.ndarray, str]:
    """Call the same default Kepler quality-mask API used by pinned ALDERAAN."""
    lightkurve = _load_pinned_lightkurve()
    result = lightkurve.KeplerQualityFlags.create_quality_mask(quality, bitmask="default")
    mask = np.asarray(result, dtype=bool)
    if mask.shape != quality.shape:
        raise PreparationInputError("lightkurve returned a quality mask with the wrong shape")
    return mask, str(lightkurve.__version__)


def _table_column(table: fits.FITS_rec, name: str, path: Path) -> np.ndarray:
    names = {str(item).upper(): str(item) for item in (table.names or ())}
    if name not in names:
        raise PreparationInputError(f"{path.name} is missing LIGHTCURVE column {name}")
    return np.asarray(table[names[name]])


def _number_array(values: np.ndarray, name: str, path: Path) -> np.ndarray:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise PreparationInputError(f"{path.name} {name} is not numeric: {error}") from error
    if array.ndim != 1:
        raise PreparationInputError(f"{path.name} {name} must be one-dimensional")
    return array


def _integer_array(values: np.ndarray, name: str, path: Path) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 1 or array.dtype.kind not in "iu":
        raise PreparationInputError(f"{path.name} {name} must be a one-dimensional integer array")
    return array.astype(np.int64, copy=False)


def _catalog_integer(value: str | None, name: str, row_number: int) -> int:
    try:
        numeric = float(value) if value is not None else np.nan
        integer = int(numeric)
    except (TypeError, ValueError, OverflowError) as error:
        raise PreparationInputError(f"catalog row {row_number} has invalid {name}") from error
    if not np.isfinite(numeric) or numeric != integer or integer <= 0:
        raise PreparationInputError(f"catalog row {row_number} has invalid positive integer {name}")
    return integer


def _catalog_number(value: str | None, name: str, row_number: int) -> float:
    try:
        numeric = float(value) if value is not None else np.nan
    except (TypeError, ValueError) as error:
        raise PreparationInputError(f"catalog row {row_number} has invalid {name}") from error
    if not np.isfinite(numeric):
        raise PreparationInputError(f"catalog row {row_number} has non-finite {name}")
    return numeric


def read_system_catalog(path: Path, target: str, kic: int, npl: int) -> dict[str, Any]:
    """Validate one explicit ALDERAAN input catalog and its full target system."""
    if not path.is_file():
        raise PreparationInputError(f"ALDERAAN system catalog does not exist: {path}")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None:
                raise PreparationInputError("ALDERAAN system catalog has no header")
            names = {name.strip().lower(): name for name in reader.fieldnames if name is not None}
            missing = [name for name in SYSTEM_CATALOG_COLUMNS if name not in names]
            if missing:
                raise PreparationInputError(
                    "ALDERAAN system catalog is missing required column(s): " + ", ".join(missing)
                )
            rows = list(reader)
    except OSError as error:
        raise PreparationInputError(f"cannot read ALDERAAN system catalog: {error}") from error
    target = target.strip().upper()
    matching: list[dict[str, Any]] = []
    kics_for_target: set[int] = set()
    targets_for_kic: set[str] = set()
    for row_number, raw_row in enumerate(rows, start=2):
        normalized = {key: raw_row[value] for key, value in names.items()}
        row_target = str(normalized["koi_id"] or "").strip().upper()
        row_kic = _catalog_integer(normalized["kic_id"], "kic_id", row_number)
        row_npl = _catalog_integer(normalized["npl"], "npl", row_number)
        if row_target == target:
            kics_for_target.add(row_kic)
            matching.append(
                {
                    "row_number": row_number,
                    "koi_id": row_target,
                    "kic_id": row_kic,
                    "npl": row_npl,
                    "period_days": _catalog_number(normalized["period"], "period", row_number),
                    "epoch_bkjd": _catalog_number(normalized["epoch"], "epoch", row_number),
                    "duration_hours": _catalog_number(normalized["duration"], "duration", row_number),
                }
            )
        if row_kic == kic:
            targets_for_kic.add(row_target)
    if not matching:
        raise PreparationInputError(f"ALDERAAN system catalog has no rows for {target}")
    if kics_for_target != {kic}:
        raise PreparationInputError(
            f"ALDERAAN system catalog maps {target} to KIC values {sorted(kics_for_target)}, expected {kic}"
        )
    if targets_for_kic != {target}:
        raise PreparationInputError(
            f"ALDERAAN system catalog maps KIC {kic} to KOI identifiers {sorted(targets_for_kic)}, expected {target}"
        )
    if len(matching) != npl:
        raise PreparationInputError(
            f"ALDERAAN system catalog has {len(matching)} rows for {target}, but result NPL is {npl}"
        )
    row_npl_values = {row["npl"] for row in matching}
    if row_npl_values != {npl}:
        raise PreparationInputError(
            f"ALDERAAN system catalog NPL values {sorted(row_npl_values)} do not match result NPL {npl}"
        )
    if len({row["period_days"] for row in matching}) != len(matching):
        raise PreparationInputError("ALDERAAN system catalog has duplicate planet periods for one target")
    if any(row["duration_hours"] <= 0.0 for row in matching):
        raise PreparationInputError("ALDERAAN system catalog duration must be positive for every companion")
    return {
        "path": path,
        "sha256": sha256_file(path),
        "bytes": int(path.stat().st_size),
        "target": target,
        "kic": kic,
        "npl": npl,
        "system_rows": matching,
        "maximum_duration_days": max(row["duration_hours"] for row in matching) / 24.0,
    }


def read_pdcsap_file(path: Path, expected_kic: int, source_file_index: int) -> dict[str, Any]:
    """Read one local PDCSAP file and retain only mask-passing finite cadences."""
    if not path.is_file():
        raise PreparationInputError(f"PDCSAP FITS does not exist: {path}")
    try:
        # MAST Kepler PDCSAP files can carry extension checksum cards that do
        # not validate even for a fresh byte-identical download. The immutable
        # source-file SHA-256 in the manifest is the provenance check here.
        with fits.open(path, memmap=False, checksum=False) as hdul:
            kic = _integer_header(hdul[0].header, ("KEPLERID", "KICID"), "KIC identifier")
            if kic != expected_kic:
                raise PreparationInputError(
                    f"{path.name} has KIC {kic}, expected {expected_kic}"
                )
            if len(hdul) < 2 or not isinstance(hdul[1], fits.BinTableHDU) or hdul[1].data is None:
                raise PreparationInputError(f"{path.name} has no populated LIGHTCURVE binary table")
            table = hdul[1].data
            ext_header = hdul[1].header
            time_convention = _kepler_time_convention(hdul[0].header, ext_header, path)
            cadence_headers = _cadence_headers(hdul[0].header, ext_header, path)
            time = _number_array(_table_column(table, "TIME", path), "TIME", path)
            flux = _number_array(_table_column(table, "PDCSAP_FLUX", path), "PDCSAP_FLUX", path)
            flux_err = _number_array(
                _table_column(table, "PDCSAP_FLUX_ERR", path), "PDCSAP_FLUX_ERR", path
            )
            cadence = _integer_array(_table_column(table, "CADENCENO", path), "CADENCENO", path)
            quality = _integer_array(_table_column(table, "SAP_QUALITY", path), "SAP_QUALITY", path)
            if not (len(time) == len(flux) == len(flux_err) == len(cadence) == len(quality)):
                raise PreparationInputError(f"{path.name} LIGHTCURVE columns have unequal lengths")
            quarter = _quarter_header(hdul[0].header, ext_header, path)
            quality_pass, lightkurve_version = default_kepler_quality_mask(quality)
            finite = np.isfinite(time) & np.isfinite(flux) & np.isfinite(flux_err)
            keep = quality_pass & finite
            rows = np.flatnonzero(keep).astype(np.int64)
            return {
                "path": path,
                "sha256": sha256_file(path),
                "bytes": int(path.stat().st_size),
                "kic": kic,
                "quarter": quarter,
                "time_convention": time_convention,
                **cadence_headers,
                "raw_count": int(len(time)),
                "quality_pass_count": int(np.sum(quality_pass)),
                "finite_count": int(np.sum(finite)),
                "retained_count": int(len(rows)),
                "lightkurve_version": lightkurve_version,
                "time": time[rows],
                "flux": flux[rows],
                "flux_err": flux_err[rows],
                "cadence": cadence[rows],
                "quality": quality[rows],
                "quarter_array": np.full(len(rows), quarter, dtype=np.int16),
                "source_file_index": np.full(len(rows), source_file_index, dtype=np.int32),
                "source_row_index": rows,
            }
    except (OSError, ValueError, TypeError) as error:
        if isinstance(error, PreparationInputError):
            raise
        raise PreparationInputError(f"cannot read {path}: {error}") from error


def _required_timing_column(table: fits.FITS_rec, name: str, extension: str) -> np.ndarray:
    names = {str(item).upper(): str(item) for item in (table.names or ())}
    if name not in names:
        raise PreparationInputError(f"{extension} is missing {name}")
    try:
        values = np.asarray(table[names[name]], dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise PreparationInputError(f"{extension} {name} is not numeric: {error}") from error
    if values.ndim != 1 or not np.all(np.isfinite(values)):
        raise PreparationInputError(f"{extension} {name} must be finite and one-dimensional")
    return values


def _validate_timing_table(timing_data: dict[str, np.ndarray], extension: str) -> None:
    index = timing_data["index"]
    ttime = timing_data["ttime"]
    model = timing_data["model"]
    out_prob = timing_data["out_prob"]
    out_flag = timing_data["out_flag"]
    if len(index) < 2:
        raise PreparationInputError(f"{extension} must contain at least two timing rows")
    if np.any(index < 0) or not np.allclose(index, np.round(index), rtol=0.0, atol=1e-9):
        raise PreparationInputError(f"{extension} INDEX must contain non-negative integers")
    if not np.all(np.diff(index) > 0):
        raise PreparationInputError(f"{extension} INDEX must be strictly increasing")
    if not np.all(np.diff(ttime) > 0):
        raise PreparationInputError(f"{extension} TTIME must be strictly increasing")
    if not np.all(np.diff(model) > 0):
        raise PreparationInputError(f"{extension} MODEL must be strictly increasing")
    if np.any((out_prob < 0.0) | (out_prob > 1.0)):
        raise PreparationInputError(f"{extension} OUT_PROB must lie in [0, 1]")
    if not np.all(np.isin(out_flag, (0.0, 1.0))):
        raise PreparationInputError(f"{extension} OUT_FLAG must contain only 0 or 1")


def _timing_convention(
    time_reference: str, timeunit: str, timesys: str
) -> dict[str, str]:
    reference = time_reference.strip().upper()
    unit = timeunit.strip().lower()
    system = timesys.strip().upper()
    if reference != "BKJD" or unit != "d" or system != "TDB":
        raise PreparationInputError(
            "ALDERAAN timing convention must be declared exactly as BKJD, d, TDB"
        )
    return {"time_reference": reference, "timeunit": unit, "timesys": system}


def read_result_timings(
    path: Path, target: str, timing_time_reference: str, timing_timeunit: str, timing_timesys: str
) -> dict[str, Any]:
    """Read all declared ALDERAAN model transit times without inferring new ones."""
    if not path.is_file():
        raise PreparationInputError(f"ALDERAAN result FITS does not exist: {path}")
    target = target.strip().upper()
    if not target:
        raise PreparationInputError("target is required")
    convention = _timing_convention(timing_time_reference, timing_timeunit, timing_timesys)
    try:
        with fits.open(path, memmap=False, checksum=False) as hdul:
            fits_target = str(hdul[0].header.get("TARGET", "")).strip().upper()
            if fits_target != target:
                raise PreparationInputError(
                    f"ALDERAAN result TARGET={fits_target or '<missing>'}, expected {target}"
                )
            npl = _integer_header(hdul[0].header, ("NPL",), "NPL")
            by_name: dict[str, list[fits.hdu.base.ExtensionHDU]] = {}
            for hdu in hdul[1:]:
                name = str(hdu.name).upper()
                if TTIMES_PATTERN.fullmatch(name):
                    by_name.setdefault(name, []).append(hdu)
            expected = {f"TTIMES_{index:02d}" for index in range(npl)}
            unexpected = sorted(set(by_name) - expected)
            if unexpected:
                raise PreparationInputError("unexpected timing extension(s): " + ", ".join(unexpected))
            planets: list[dict[str, Any]] = []
            all_models: list[np.ndarray] = []
            for index in range(npl):
                name = f"TTIMES_{index:02d}"
                matches = by_name.get(name, [])
                if len(matches) != 1 or not isinstance(matches[0], fits.BinTableHDU) or matches[0].data is None:
                    raise PreparationInputError(f"{name} must occur once as a populated binary table")
                table = matches[0].data
                timing_data = {
                    column.lower(): _required_timing_column(table, column, name)
                    for column in TIMING_COLUMNS
                }
                _validate_timing_table(timing_data, name)
                model = timing_data["model"]
                if len(model) == 0:
                    raise PreparationInputError(f"{name} has no model transit times")
                record: dict[str, Any] = {
                    "planet_index": index,
                    "extension": name,
                    **{column: values.tolist() for column, values in timing_data.items()},
                }
                planets.append(record)
                all_models.append(model)
            return {
                "path": path,
                "sha256": sha256_file(path),
                "bytes": int(path.stat().st_size),
                "target": fits_target,
                "npl": npl,
                "time_convention": convention,
                "planets": planets,
                "model_times": np.sort(np.concatenate(all_models)),
            }
    except (OSError, ValueError, TypeError) as error:
        if isinstance(error, PreparationInputError):
            raise
        raise PreparationInputError(f"cannot read ALDERAAN result {path}: {error}") from error


def conservative_transit_mask(time: np.ndarray, model_times: np.ndarray, window_days: float) -> np.ndarray:
    """Mark cadences within a declared symmetric interval of any model transit."""
    if not np.isfinite(window_days) or window_days <= 0:
        raise PreparationInputError("transit_window_days must be finite and positive")
    if time.ndim != 1 or not np.all(np.isfinite(time)) or np.any(np.diff(time) <= 0):
        raise PreparationInputError("time must be finite and strictly increasing before transit masking")
    if model_times.ndim != 1 or len(model_times) == 0 or not np.all(np.isfinite(model_times)):
        raise PreparationInputError("model transit times must be a non-empty finite array")
    mask = np.zeros(len(time), dtype=bool)
    for transit_time in model_times:
        left = int(np.searchsorted(time, transit_time - window_days, side="left"))
        right = int(np.searchsorted(time, transit_time + window_days, side="right"))
        mask[left:right] = True
    return mask


def transit_window_half_width(
    duration_days: float, multiple: float, timing_uncertainty_margin_days: float
) -> float:
    """Compute the declared conservative mask half-width without hidden defaults."""
    values = {
        "duration_days": duration_days,
        "multiple": multiple,
        "timing_uncertainty_margin_days": timing_uncertainty_margin_days,
    }
    for name, value in values.items():
        if not np.isfinite(value) or value < 0:
            raise PreparationInputError(f"{name} must be finite and non-negative")
    if duration_days <= 0.0 or multiple < 1.0:
        raise PreparationInputError("duration_days must be positive and multiple must be at least one")
    half_width = duration_days * multiple + timing_uncertainty_margin_days
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise PreparationInputError("computed transit window half-width must be positive")
    return float(half_width)


def _concatenate(records: Sequence[dict[str, Any]], field: str, dtype: Any) -> np.ndarray:
    return np.concatenate([np.asarray(record[field], dtype=dtype) for record in records])


def prepare_input(
    pdcsap_paths: Sequence[Path],
    result_path: Path,
    system_catalog_path: Path,
    target: str,
    kic: int,
    duration_days: float,
    duration_multiple: float,
    timing_uncertainty_margin_days: float,
    cadence_mode: str,
    exposure_minutes: float,
    timing_time_reference: str,
    timing_timeunit: str,
    timing_timesys: str,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Build the immutable NPZ arrays and its JSON-serializable provenance."""
    if not pdcsap_paths:
        raise PreparationInputError("at least one local PDCSAP FITS file is required")
    if kic <= 0:
        raise PreparationInputError("kic must be a positive integer")
    if cadence_mode not in {"long", "short"}:
        raise PreparationInputError("cadence_mode must be 'long' or 'short'")
    if not np.isfinite(exposure_minutes) or exposure_minutes <= 0.0:
        raise PreparationInputError("exposure_minutes must be finite and positive")
    result = read_result_timings(
        result_path, target, timing_time_reference, timing_timeunit, timing_timesys
    )
    system_catalog = read_system_catalog(system_catalog_path, result["target"], kic, result["npl"])
    half_width = transit_window_half_width(
        duration_days, duration_multiple, timing_uncertainty_margin_days
    )
    if duration_days < system_catalog["maximum_duration_days"]:
        raise PreparationInputError(
            "duration_days must be at least the maximum duration in the declared full system catalog"
        )
    records = [read_pdcsap_file(Path(path), kic, index) for index, path in enumerate(pdcsap_paths)]
    versions = {record["lightkurve_version"] for record in records}
    if versions != {PINNED_LIGHTKURVE_VERSION}:
        raise PreparationInputError("PDCSAP files were not masked with one pinned lightkurve version")
    observed_modes = {record["observed_mode"] for record in records}
    if observed_modes != {cadence_mode}:
        raise PreparationInputError(
            "PDCSAP OBSMODE values "
            f"{sorted(observed_modes)} do not match declared cadence_mode {cadence_mode!r}"
        )
    observed_timedel_days = [float(record["timedel_days"]) for record in records]
    reference_timedel_days = observed_timedel_days[0]
    if any(not _same_timedelta(reference_timedel_days, value) for value in observed_timedel_days[1:]):
        raise PreparationInputError("PDCSAP files have inconsistent TIMEDEL values")
    observed_exposure_minutes = reference_timedel_days * 1440.0
    if not np.isclose(
        exposure_minutes,
        observed_exposure_minutes,
        rtol=0.0,
        atol=EXPOSURE_MATCH_TOLERANCE_MINUTES,
    ):
        raise PreparationInputError(
            "declared exposure_minutes does not match PDCSAP TIMEDEL: "
            f"declared {exposure_minutes:.10g}, observed {observed_exposure_minutes:.10g}"
        )
    conventions = {tuple(sorted(record["time_convention"].items())) for record in records}
    if len(conventions) != 1:
        raise PreparationInputError("PDCSAP files have contradictory time conventions")
    input_convention = records[0]["time_convention"]
    shared_convention = {
        key: input_convention[key] for key in ("time_reference", "timeunit", "timesys")
    }
    if shared_convention != result["time_convention"]:
        raise PreparationInputError("PDCSAP and declared ALDERAAN timing conventions do not match")
    arrays = {
        name: _concatenate(records, name, dtype)
        for name, dtype in {
            "time": np.float64,
            "flux": np.float64,
            "flux_err": np.float64,
            "cadence": np.int64,
            "quality": np.int64,
            "quarter_array": np.int16,
            "source_file_index": np.int32,
            "source_row_index": np.int64,
        }.items()
    }
    arrays["quarter"] = arrays.pop("quarter_array")
    order = np.argsort(arrays["time"], kind="stable")
    arrays = {name: values[order] for name, values in arrays.items()}
    if len(arrays["time"]) == 0 or np.any(np.diff(arrays["time"]) <= 0):
        raise PreparationInputError("retained PDCSAP cadence times must be non-empty and strictly increasing")
    known_transit_mask = conservative_transit_mask(
        arrays["time"], result["model_times"], half_width
    )
    mask_index = np.flatnonzero(~known_transit_mask).astype(np.int64)
    if len(mask_index) == 0:
        raise PreparationInputError("known-transit window excludes every retained cadence")
    arrays["known_transit_mask"] = known_transit_mask
    arrays["mask_index"] = mask_index
    arrays["source_files"] = np.asarray([str(record["path"].resolve()) for record in records], dtype="U")
    manifest = {
        "schema_version": 1,
        "scope": "Standalone injection-input preparation. No detrending, fitted-model subtraction, transit injection, fitting, or ALDERAAN result changes were performed.",
        "target": result["target"],
        "kic": int(kic),
        "cadence": {
            "mode": cadence_mode,
            "exposure_minutes": float(exposure_minutes),
            "observed_mode": records[0]["observed_mode"],
            "timedel_days": reference_timedel_days,
            "observed_exposure_minutes": observed_exposure_minutes,
            "exposure_match_tolerance_minutes": EXPOSURE_MATCH_TOLERANCE_MINUTES,
        },
        "time_convention": input_convention,
        "system_catalog": {
            "filename": system_catalog["path"].name,
            "sha256": system_catalog["sha256"],
            "bytes": system_catalog["bytes"],
            "target": system_catalog["target"],
            "kic": system_catalog["kic"],
            "npl": system_catalog["npl"],
            "maximum_duration_days": system_catalog["maximum_duration_days"],
            "system_rows": system_catalog["system_rows"],
        },
        "quality_mask": {
            "implementation": "lightkurve.KeplerQualityFlags.create_quality_mask",
            "bitmask": "default",
            "pinned_lightkurve_version": PINNED_LIGHTKURVE_VERSION,
        },
        "known_transit_exclusion": {
            "duration_days": float(duration_days),
            "multiple": float(duration_multiple),
            "timing_uncertainty_margin_days": float(timing_uncertainty_margin_days),
            "half_width_days": half_width,
            "rule": "Cadences within or on the declared symmetric window around every TTIMES_nn MODEL time are excluded from mask_index.",
            "result_fits": {
                "filename": result["path"].name,
                "sha256": result["sha256"],
                "bytes": result["bytes"],
                "target": result["target"],
                "npl": result["npl"],
                "timing_tables": result["planets"],
            },
        },
        "source_files": [
            {
                key: str(record[key]) if key == "path" else record[key]
                for key in (
                    "path", "sha256", "bytes", "kic", "quarter", "observed_mode",
                    "timedel_days", "observed_exposure_minutes", "raw_count",
                    "quality_pass_count", "finite_count", "retained_count",
                )
            }
            for record in records
        ],
        "cadence_counts": {
            "raw_total": int(sum(record["raw_count"] for record in records)),
            "quality_pass_total": int(sum(record["quality_pass_count"] for record in records)),
            "finite_total": int(sum(record["finite_count"] for record in records)),
            "retained_total": int(len(arrays["time"])),
            "known_transit_excluded_total": int(np.sum(known_transit_mask)),
            "injection_eligible_total": int(len(mask_index)),
        },
        "arrays": {
            name: {"count": int(len(values)), "sha256": sha256_array(values)}
            for name, values in arrays.items()
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "astropy": ASTROPY_VERSION,
            "lightkurve": PINNED_LIGHTKURVE_VERSION,
        },
    }
    return arrays, manifest


def _write_json_once(path: Path, payload: dict[str, Any]) -> None:
    try:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    except FileExistsError as error:
        raise PreparationInputError(f"refusing to replace immutable manifest: {path}") from error


def write_prepared_input(
    arrays: dict[str, np.ndarray], manifest: dict[str, Any], output_dir: Path
) -> tuple[Path, Path]:
    """Write a new NPZ and one manifest without replacing prior provenance."""
    _output_directory(output_dir)
    output_path = output_dir / OUTPUT_NAME
    np.savez_compressed(output_path, **arrays)
    manifest = dict(manifest)
    manifest["output"] = {"filename": output_path.name, "sha256": sha256_file(output_path)}
    manifest_path = output_dir / MANIFEST_NAME
    _write_json_once(manifest_path, manifest)
    return output_path, manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdcsap-fits", type=Path, nargs="+", required=True)
    parser.add_argument("--alderaan-result-fits", type=Path, required=True)
    parser.add_argument("--alderaan-system-catalog", type=Path, required=True)
    parser.add_argument("--target", required=True, help="Expected ALDERAAN TARGET header, such as K00367.")
    parser.add_argument("--kic", type=int, required=True, help="Expected KEPLERID/KICID in every PDCSAP primary header.")
    parser.add_argument("--duration-days", type=float, required=True, help="Declared conservative full-system duration in days.")
    parser.add_argument("--duration-multiple", type=float, required=True, help="Positive multiplier applied to duration-days for each side of the mask.")
    parser.add_argument("--timing-uncertainty-margin-days", type=float, required=True)
    parser.add_argument("--cadence-mode", choices=("long", "short"), required=True)
    parser.add_argument("--exposure-minutes", type=float, required=True)
    parser.add_argument("--timing-time-reference", required=True, help="Explicit convention for headerless ALDERAAN TTIMES, currently BKJD only.")
    parser.add_argument("--timing-timeunit", required=True, help="Explicit convention for headerless ALDERAAN TTIMES, currently d only.")
    parser.add_argument("--timing-timesys", required=True, help="Explicit convention for headerless ALDERAAN TTIMES, currently TDB only.")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        arrays, manifest = prepare_input(
            args.pdcsap_fits,
            args.alderaan_result_fits,
            args.alderaan_system_catalog,
            args.target,
            args.kic,
            args.duration_days,
            args.duration_multiple,
            args.timing_uncertainty_margin_days,
            args.cadence_mode,
            args.exposure_minutes,
            args.timing_time_reference,
            args.timing_timeunit,
            args.timing_timesys,
        )
        output_path, manifest_path = write_prepared_input(arrays, manifest, args.output_dir)
    except (PreparationInputError, LightkurveDependencyError) as error:
        build_parser().error(str(error))
    print(f"Input NPZ: {output_path}")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
