"""Audit what pinned ALDERAAN retains from a circular-injection fixture.

Run this only in the pinned ALDERAAN Linux environment. It does not fit a
transit, write a result FITS, or modify a fixture.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


FIELD_RTOL = 2.0e-7
FIELD_ATOL = 1.0e-10


class ReaderAuditError(ValueError):
    """Raised when a fixture cannot be checked against the ALDERAAN reader."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def classify_reader_times(
    prepared_time: np.ndarray, known_transit_mask: np.ndarray, reader_time: np.ndarray
) -> dict[str, int]:
    """Count whether the reader retained expected injected and masked cadences."""
    prepared_time = np.asarray(prepared_time, dtype=np.float64)
    known_transit_mask = np.asarray(known_transit_mask, dtype=bool)
    reader_time = np.asarray(reader_time, dtype=np.float64)
    if prepared_time.ndim != 1 or known_transit_mask.shape != prepared_time.shape or reader_time.ndim != 1:
        raise ReaderAuditError("time arrays and known_transit_mask must be aligned one-dimensional arrays")
    if not np.all(np.isfinite(prepared_time)) or not np.all(np.isfinite(reader_time)):
        raise ReaderAuditError("reader comparison requires finite times")
    retained = np.isin(prepared_time, reader_time)
    return {
        "prepared_total": int(len(prepared_time)),
        "eligible_total": int(np.sum(~known_transit_mask)),
        "eligible_retained": int(np.sum(retained & ~known_transit_mask)),
        "eligible_missing": int(np.sum(~retained & ~known_transit_mask)),
        "known_transit_total": int(np.sum(known_transit_mask)),
        "known_transit_retained": int(np.sum(retained & known_transit_mask)),
        "reader_total": int(len(reader_time)),
        "reader_extra": int(np.sum(~np.isin(reader_time, prepared_time))),
    }


def _numeric_values(value: object, dtype: type[np.generic]) -> np.ndarray:
    """Return a numeric Lightkurve value array without changing its order."""
    raw = getattr(value, "value", value)
    return np.asarray(raw, dtype=dtype)


def _quarter_values(value: object, count: int) -> np.ndarray:
    """Expand a scalar quarter or retain a per-cadence quarter array."""
    quarter = _numeric_values(value, np.int64)
    if quarter.ndim == 0:
        return np.full(count, int(quarter), dtype=np.int64)
    if quarter.ndim == 1 and len(quarter) == count:
        return quarter
    raise ReaderAuditError("reader quarter values must be scalar or align with reader times")


def _reader_arrays(collection: object) -> dict[str, np.ndarray]:
    """Collect reader fields before any detrending or transit-model fitting."""
    chunks: dict[str, list[np.ndarray]] = {"time": [], "flux": []}
    optional = {"flux_err": [], "cadence": [], "quarter": []}
    seen_optional = {name: True for name in optional}
    for lightcurve in collection:
        time = _numeric_values(getattr(lightcurve, "time"), np.float64)
        flux = _numeric_values(getattr(lightcurve, "flux"), np.float64)
        if time.ndim != 1 or flux.ndim != 1 or len(time) != len(flux):
            raise ReaderAuditError("reader time and flux must be aligned one-dimensional arrays")
        chunks["time"].append(time)
        chunks["flux"].append(flux)
        attributes = {
            "flux_err": "flux_err",
            "cadence": "cadenceno",
            "quarter": "quarter",
        }
        for field, attribute in attributes.items():
            if not hasattr(lightcurve, attribute):
                seen_optional[field] = False
                continue
            if field == "quarter":
                values = _quarter_values(getattr(lightcurve, attribute), len(time))
            else:
                values = _numeric_values(getattr(lightcurve, attribute), np.float64 if field == "flux_err" else np.int64)
                if values.ndim != 1 or len(values) != len(time):
                    raise ReaderAuditError(f"reader {field} must align with reader times")
            optional[field].append(values)
    if not chunks["time"]:
        raise ReaderAuditError("pinned ALDERAAN reader returned no light curves")
    payload = {name: np.concatenate(values) for name, values in chunks.items()}
    for field, values in optional.items():
        if seen_optional[field] and len(values) == len(chunks["time"]):
            payload[field] = np.concatenate(values)
    return payload


def _matched_indices(prepared_time: np.ndarray, reader_time: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pair exact BKJD values while preserving reader output order."""
    lookup = {float(value): index for index, value in enumerate(prepared_time)}
    reader_index = np.array([index for index, value in enumerate(reader_time) if float(value) in lookup], dtype=np.int64)
    prepared_index = np.array([lookup[float(reader_time[index])] for index in reader_index], dtype=np.int64)
    if len(np.unique(prepared_index)) != len(prepared_index):
        raise ReaderAuditError("reader duplicated a fixture cadence time")
    return prepared_index, reader_index


def _check_exact_field(
    name: str, expected: np.ndarray, actual: np.ndarray, prepared_index: np.ndarray, reader_index: np.ndarray
) -> dict[str, object]:
    """Check a discrete field only at fixture cadences the reader retained."""
    if len(prepared_index) == 0:
        return {"applicable": True, "matched": 0, "pass": True}
    equal = np.array_equal(expected[prepared_index], actual[reader_index])
    return {"applicable": True, "matched": int(len(prepared_index)), "pass": bool(equal)}


def _quarter_scales(
    expected_flux: np.ndarray, actual_flux: np.ndarray, quarter: np.ndarray
) -> tuple[dict[str, float], np.ndarray]:
    """Infer stitch normalization scales from finite, nonzero fixture fluxes."""
    scales: dict[str, float] = {}
    model = np.empty_like(actual_flux, dtype=np.float64)
    for value in np.unique(quarter):
        use = quarter == value
        reference = expected_flux[use]
        observed = actual_flux[use]
        support = np.isfinite(reference) & np.isfinite(observed) & (np.abs(reference) > FIELD_ATOL)
        if not np.any(support):
            scale = 1.0
        else:
            scale = float(np.median(observed[support] / reference[support]))
        if not np.isfinite(scale) or scale <= 0.0:
            raise ReaderAuditError("reader flux cannot support a positive stitch-normalization scale")
        scales[str(int(value))] = scale
        model[use] = reference * scale
    return scales, model


def _check_scaled_field(
    name: str,
    expected: np.ndarray,
    actual: np.ndarray,
    prepared_index: np.ndarray,
    reader_index: np.ndarray,
    quarter: np.ndarray,
    scales: dict[str, float],
) -> dict[str, object]:
    """Check continuous reader values after the same per-quarter stitch scale."""
    if len(prepared_index) == 0:
        return {"applicable": True, "matched": 0, "pass": True, "max_abs_residual": 0.0}
    scale = np.array([scales[str(int(value))] for value in quarter], dtype=np.float64)
    model = expected[prepared_index] * scale
    observed = actual[reader_index]
    finite = np.isfinite(model) & np.isfinite(observed)
    passed = bool(np.all(finite) and np.allclose(observed, model, rtol=FIELD_RTOL, atol=FIELD_ATOL))
    residual = np.abs(observed - model)
    return {
        "applicable": True,
        "matched": int(len(prepared_index)),
        "pass": passed,
        "max_abs_residual": float(np.max(residual)) if len(residual) else 0.0,
    }


def compare_reader_provenance(
    prepared: dict[str, np.ndarray], reader: dict[str, np.ndarray]
) -> dict[str, object]:
    """Compare reader cadences with the injected fixture, not a fitted model."""
    required = {"time", "flux", "known_transit_mask"}
    missing = required - set(prepared)
    if missing:
        raise ReaderAuditError("fixture provenance is missing " + ", ".join(sorted(missing)))
    if "time" not in reader or "flux" not in reader:
        raise ReaderAuditError("reader output is missing time or flux")
    time = np.asarray(prepared["time"], dtype=np.float64)
    known = np.asarray(prepared["known_transit_mask"], dtype=bool)
    if time.ndim != 1 or known.shape != time.shape or np.any(~np.isfinite(time)):
        raise ReaderAuditError("fixture time and known_transit_mask must be aligned finite arrays")
    if len(np.unique(time)) != len(time):
        raise ReaderAuditError("fixture provenance repeats a cadence time")
    reader_time = np.asarray(reader["time"], dtype=np.float64)
    counts = classify_reader_times(time, known, reader_time)
    prepared_index, reader_index = _matched_indices(time, reader_time)
    eligible = ~known[prepared_index]
    prepared_index = prepared_index[eligible]
    reader_index = reader_index[eligible]
    checks: dict[str, dict[str, object]] = {}
    expected_quarter = prepared.get("quarter")
    reader_quarter = reader.get("quarter")
    if expected_quarter is None:
        checks["quarter"] = {"applicable": False, "reason": "fixture provenance has no quarter array"}
        comparison_quarter = np.zeros(len(prepared_index), dtype=np.int64)
    elif reader_quarter is None:
        checks["quarter"] = {"applicable": True, "matched": int(len(prepared_index)), "pass": False, "reason": "reader has no quarter values"}
        comparison_quarter = np.zeros(len(prepared_index), dtype=np.int64)
    else:
        expected_quarter = np.asarray(expected_quarter, dtype=np.int64)
        reader_quarter = np.asarray(reader_quarter, dtype=np.int64)
        if expected_quarter.shape != time.shape or reader_quarter.shape != reader_time.shape:
            raise ReaderAuditError("quarter arrays must align with their time arrays")
        checks["quarter"] = _check_exact_field("quarter", expected_quarter, reader_quarter, prepared_index, reader_index)
        comparison_quarter = expected_quarter[prepared_index]
    for field in ("cadence",):
        expected = prepared.get(field)
        actual = reader.get(field)
        if expected is None:
            checks[field] = {"applicable": False, "reason": f"fixture provenance has no {field} array"}
        elif actual is None:
            checks[field] = {"applicable": True, "matched": int(len(prepared_index)), "pass": False, "reason": f"reader has no {field} values"}
        else:
            expected = np.asarray(expected, dtype=np.int64)
            actual = np.asarray(actual, dtype=np.int64)
            if expected.shape != time.shape or actual.shape != reader_time.shape:
                raise ReaderAuditError(f"{field} arrays must align with their time arrays")
            checks[field] = _check_exact_field(field, expected, actual, prepared_index, reader_index)
    expected_flux = np.asarray(prepared["flux"], dtype=np.float64)
    actual_flux = np.asarray(reader["flux"], dtype=np.float64)
    if expected_flux.shape != time.shape or actual_flux.shape != reader_time.shape:
        raise ReaderAuditError("flux arrays must align with their time arrays")
    scales, _ = _quarter_scales(expected_flux[prepared_index], actual_flux[reader_index], comparison_quarter)
    checks["flux"] = _check_scaled_field("flux", expected_flux, actual_flux, prepared_index, reader_index, comparison_quarter, scales)
    expected_error = prepared.get("flux_err")
    actual_error = reader.get("flux_err")
    if expected_error is None:
        checks["flux_err"] = {"applicable": False, "reason": "fixture provenance has no flux_err array"}
    elif actual_error is None:
        checks["flux_err"] = {"applicable": True, "matched": int(len(prepared_index)), "pass": False, "reason": "reader has no flux_err values"}
    else:
        expected_error = np.asarray(expected_error, dtype=np.float64)
        actual_error = np.asarray(actual_error, dtype=np.float64)
        if expected_error.shape != time.shape or actual_error.shape != reader_time.shape:
            raise ReaderAuditError("flux_err arrays must align with their time arrays")
        checks["flux_err"] = _check_scaled_field("flux_err", expected_error, actual_error, prepared_index, reader_index, comparison_quarter, scales)
    checked = [item["pass"] for item in checks.values() if item.get("applicable")]
    return {
        "counts": counts,
        "field_checks": checks,
        "quarter_stitch_scales": scales,
        "pass": counts["eligible_missing"] == 0 and counts["known_transit_retained"] == 0 and all(checked),
    }


def audit_fixture(fixture_dir: Path, injection_path: Path, kic: int, output: Path) -> Path:
    """Read a fixture through ALDERAAN and write a non-overwriting JSON audit."""
    if output.exists():
        raise ReaderAuditError(f"refusing to replace existing audit: {output}")
    manifest = fixture_dir / "diagnostic_fixture_manifest.json"
    data_dir = fixture_dir / "Data"
    if not manifest.is_file() or not data_dir.is_dir() or not injection_path.is_file():
        raise ReaderAuditError("fixture manifest, Data directory, or injection NPZ is missing")
    try:
        fixture_record = json.loads(manifest.read_text(encoding="utf-8"))
        expected_injection_hash = fixture_record["injection"]["sha256"]
        if expected_injection_hash != _sha256(injection_path):
            raise ReaderAuditError("fixture manifest does not identify this injection NPZ")
        with np.load(injection_path, allow_pickle=False) as payload:
            prepared = {name: np.asarray(payload[name]) for name in payload.files}
    except (KeyError, OSError, ValueError, TypeError, json.JSONDecodeError) as error:
        raise ReaderAuditError("fixture manifest or injection NPZ is unreadable") from error
    try:
        from alderaan import io
    except ModuleNotFoundError as error:
        raise ReaderAuditError("pinned ALDERAAN must be on PYTHONPATH") from error
    files = sorted(data_dir.glob(f"kplr{kic:09d}*.fits"))
    if not files:
        raise ReaderAuditError(f"fixture has no MAST FITS for KIC {kic}")
    collection = io.read_mast_files([str(path) for path in files], kic, "long cadence")
    comparison = compare_reader_provenance(prepared, _reader_arrays(collection))
    payload = {
        "schema_version": 2,
        "scope": "Pinned-ALDERAAN reader audit for a diagnostic circular-injection fixture. It reads and compares fixture provenance only; it does not detrend data or fit a transit model.",
        "fixture_manifest_sha256": _sha256(manifest),
        "injection_sha256": _sha256(injection_path),
        "kic": int(kic),
        "source_file_count": len(files),
        **comparison,
        "interpretation": "Extra reader cadences are reported, not ignored. Flux and flux_err are checked after a common positive per-quarter scale because Lightkurve stitch may normalize each quarter; cadence number and quarter remain exact checks.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-dir", type=Path, required=True)
    parser.add_argument("--injection", type=Path, required=True)
    parser.add_argument("--kic", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    result = audit_fixture(args.fixture_dir, args.injection, args.kic, args.output)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
