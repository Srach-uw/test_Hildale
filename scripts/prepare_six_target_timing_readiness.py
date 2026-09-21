"""Validate embedded ALDERAAN transit timings for the six-target pilot.

This preparation step reads result FITS files only. It does not download light
curves, reconstruct missing timings, or perform transit injection or fitting.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from astropy.io import fits


PILOT_TARGETS = (
    "K00367",
    "K01868",
    "K01474",
    "K02204",
    "K00366",
    "K01852",
)
REQUIRED_COLUMNS = ("INDEX", "TTIME", "MODEL", "OUT_PROB", "OUT_FLAG")
ENV_VAR = "SAGEAR_INJECTION_RESULT_FITS_JSON"
SCHEMA_VERSION = "1.0"
TTIMES_PATTERN = re.compile(r"^TTIMES_(\d+)$")


class ReadinessInputError(ValueError):
    """Raised when CLI inputs cannot define the frozen six-target audit."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _issue(code: str, message: str, **context: Any) -> dict[str, Any]:
    return {"code": code, "message": message, **context}


def parse_result_specifications(specifications: Sequence[str]) -> dict[str, Path]:
    mapping: dict[str, Path] = {}
    for specification in specifications:
        if "=" not in specification:
            raise ReadinessInputError(
                f"Invalid --result value {specification!r}; expected TARGET=PATH"
            )
        target, raw_path = specification.split("=", 1)
        target = target.strip().upper()
        raw_path = raw_path.strip()
        if not target or not raw_path:
            raise ReadinessInputError(
                f"Invalid --result value {specification!r}; target and path are required"
            )
        if target in mapping:
            raise ReadinessInputError(f"Duplicate result path for {target}")
        mapping[target] = Path(raw_path).expanduser()
    return validate_target_mapping(mapping)


def parse_environment_mapping(raw_json: str) -> dict[str, Path]:
    try:
        payload = json.loads(raw_json)
    except json.JSONDecodeError as error:
        raise ReadinessInputError(f"{ENV_VAR} is not valid JSON: {error}") from error
    if not isinstance(payload, dict):
        raise ReadinessInputError(f"{ENV_VAR} must contain a JSON object")
    mapping: dict[str, Path] = {}
    for raw_target, raw_path in payload.items():
        target = str(raw_target).strip().upper()
        if not isinstance(raw_path, str) or not raw_path.strip():
            raise ReadinessInputError(f"{ENV_VAR} path for {target!r} must be a string")
        if target in mapping:
            raise ReadinessInputError(f"Duplicate result path for {target}")
        mapping[target] = Path(raw_path).expanduser()
    return validate_target_mapping(mapping)


def validate_target_mapping(mapping: Mapping[str, Path]) -> dict[str, Path]:
    expected = set(PILOT_TARGETS)
    supplied = set(mapping)
    missing = sorted(expected - supplied)
    unexpected = sorted(supplied - expected)
    if missing or unexpected:
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if unexpected:
            details.append("unexpected " + ", ".join(unexpected))
        raise ReadinessInputError("Result mapping must contain exactly the pilot targets: " + "; ".join(details))
    return {target: Path(mapping[target]) for target in PILOT_TARGETS}


def _numeric_column(
    table: fits.FITS_rec,
    column: str,
    extension: str,
    planet_index: int,
    issues: list[dict[str, Any]],
) -> np.ndarray | None:
    try:
        values = np.asarray(table[column], dtype=float)
    except (TypeError, ValueError) as error:
        issues.append(
            _issue(
                "non_numeric_column",
                f"{extension} column {column} is not numeric: {error}",
                planet_index=planet_index,
                extension=extension,
                column=column,
            )
        )
        return None
    if values.ndim != 1:
        issues.append(
            _issue(
                "non_scalar_column",
                f"{extension} column {column} is not one-dimensional",
                planet_index=planet_index,
                extension=extension,
                column=column,
            )
        )
        return None
    if not np.all(np.isfinite(values)):
        issues.append(
            _issue(
                "non_finite_column",
                f"{extension} column {column} contains non-finite values",
                planet_index=planet_index,
                extension=extension,
                column=column,
            )
        )
        return None
    return values


def validate_ttimes_extension(
    hdu: fits.hdu.base.ExtensionHDU,
    planet_index: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    extension = f"TTIMES_{planet_index:02d}"
    summary: dict[str, Any] = {
        "planet_index": planet_index,
        "extension": extension,
        "present": True,
        "ready": False,
        "rows": 0,
    }
    issues: list[dict[str, Any]] = []
    if not isinstance(hdu, fits.BinTableHDU) or hdu.data is None:
        issues.append(
            _issue(
                "invalid_ttimes_hdu",
                f"{extension} is not a populated binary table",
                planet_index=planet_index,
                extension=extension,
            )
        )
        return summary, issues

    table = hdu.data
    summary["rows"] = int(len(table))
    names = tuple(table.names or ())
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in names]
    if missing_columns:
        issues.append(
            _issue(
                "missing_ttimes_columns",
                f"{extension} is missing required columns: {', '.join(missing_columns)}",
                planet_index=planet_index,
                extension=extension,
                columns=missing_columns,
            )
        )
        return summary, issues
    if len(table) < 2:
        issues.append(
            _issue(
                "insufficient_timing_rows",
                f"{extension} has {len(table)} row(s); at least two are required",
                planet_index=planet_index,
                extension=extension,
            )
        )
        return summary, issues

    arrays = {
        column: _numeric_column(table, column, extension, planet_index, issues)
        for column in REQUIRED_COLUMNS
    }
    if issues:
        return summary, issues

    index = arrays["INDEX"]
    ttime = arrays["TTIME"]
    model = arrays["MODEL"]
    out_prob = arrays["OUT_PROB"]
    out_flag = arrays["OUT_FLAG"]
    assert index is not None and ttime is not None and model is not None
    assert out_prob is not None and out_flag is not None

    if np.any(index < 0) or not np.allclose(index, np.round(index), rtol=0.0, atol=1e-9):
        issues.append(
            _issue(
                "invalid_transit_index",
                f"{extension} INDEX must contain non-negative integers",
                planet_index=planet_index,
                extension=extension,
            )
        )
    if not np.all(np.diff(index) > 0):
        issues.append(
            _issue(
                "nonmonotonic_transit_index",
                f"{extension} INDEX must be strictly increasing; gaps are allowed",
                planet_index=planet_index,
                extension=extension,
            )
        )
    if not np.all(np.diff(ttime) > 0):
        issues.append(
            _issue(
                "nonmonotonic_ttime",
                f"{extension} TTIME must be strictly increasing",
                planet_index=planet_index,
                extension=extension,
            )
        )
    if not np.all(np.diff(model) > 0):
        issues.append(
            _issue(
                "nonmonotonic_model",
                f"{extension} MODEL must be strictly increasing",
                planet_index=planet_index,
                extension=extension,
            )
        )
    if np.any((out_prob < 0.0) | (out_prob > 1.0)):
        issues.append(
            _issue(
                "invalid_out_probability",
                f"{extension} OUT_PROB must lie in [0, 1]",
                planet_index=planet_index,
                extension=extension,
            )
        )
    if not np.all(np.isin(out_flag, (0.0, 1.0))):
        issues.append(
            _issue(
                "invalid_out_flag",
                f"{extension} OUT_FLAG must contain only 0 or 1",
                planet_index=planet_index,
                extension=extension,
            )
        )

    summary.update(
        {
            "index_min": int(np.min(index)),
            "index_max": int(np.max(index)),
            "index_gap_count": int(np.sum(np.diff(index) > 1)),
            "ttime_min": float(np.min(ttime)),
            "ttime_max": float(np.max(ttime)),
            "model_min": float(np.min(model)),
            "model_max": float(np.max(model)),
            "out_flagged_rows": int(np.sum(out_flag == 1.0)),
            "ready": not issues,
        }
    )
    return summary, issues


def inspect_result(target: str, path: Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "target": target,
        "input_path": str(path),
        "input_sha256": None,
        "input_bytes": None,
        "fits_target": None,
        "npl": None,
        "ready": False,
        "reasons": [],
        "extensions": [],
    }
    issues: list[dict[str, Any]] = record["reasons"]
    if not path.is_file():
        issues.append(_issue("missing_input_file", f"Result FITS does not exist: {path}"))
        return record

    try:
        record["input_sha256"] = sha256_file(path)
        record["input_bytes"] = path.stat().st_size
    except OSError as error:
        issues.append(_issue("unreadable_input_file", f"Cannot hash result FITS: {error}"))
        return record

    try:
        with fits.open(path, memmap=False, checksum=True) as hdul:
            fits_target = str(hdul[0].header.get("TARGET", "")).strip().upper()
            record["fits_target"] = fits_target or None
            if not fits_target:
                issues.append(_issue("missing_target_header", "Primary header is missing TARGET"))
            elif fits_target != target:
                issues.append(
                    _issue(
                        "target_mismatch",
                        f"Expected TARGET={target}, found {fits_target}",
                    )
                )

            raw_npl = hdul[0].header.get("NPL")
            try:
                numeric_npl = float(raw_npl)
                npl = int(numeric_npl)
                valid_npl = np.isfinite(numeric_npl) and numeric_npl == npl and npl > 0
            except (TypeError, ValueError, OverflowError):
                npl = 0
                valid_npl = False
            if not valid_npl:
                issues.append(_issue("invalid_npl", f"Primary header NPL is not a positive integer: {raw_npl!r}"))
                return record
            record["npl"] = npl

            timing_hdus: dict[str, list[fits.hdu.base.ExtensionHDU]] = {}
            for hdu in hdul[1:]:
                name = str(hdu.name).upper()
                if TTIMES_PATTERN.fullmatch(name):
                    timing_hdus.setdefault(name, []).append(hdu)

            expected_names = {f"TTIMES_{index:02d}" for index in range(npl)}
            unexpected = sorted(set(timing_hdus) - expected_names)
            for name in unexpected:
                issues.append(
                    _issue(
                        "unexpected_ttimes_extension",
                        f"{name} does not correspond to a planet index in NPL={npl}",
                        extension=name,
                    )
                )

            for planet_index in range(npl):
                name = f"TTIMES_{planet_index:02d}"
                matches = timing_hdus.get(name, [])
                if not matches:
                    issues.append(
                        _issue(
                            "missing_ttimes_extension",
                            f"Required extension {name} is absent; no timing data were inferred",
                            planet_index=planet_index,
                            extension=name,
                        )
                    )
                    record["extensions"].append(
                        {
                            "planet_index": planet_index,
                            "extension": name,
                            "present": False,
                            "ready": False,
                            "rows": 0,
                        }
                    )
                    continue
                if len(matches) > 1:
                    issues.append(
                        _issue(
                            "duplicate_ttimes_extension",
                            f"Required extension {name} appears {len(matches)} times",
                            planet_index=planet_index,
                            extension=name,
                        )
                    )
                    record["extensions"].append(
                        {
                            "planet_index": planet_index,
                            "extension": name,
                            "present": True,
                            "ready": False,
                            "rows": None,
                        }
                    )
                    continue
                extension_summary, extension_issues = validate_ttimes_extension(
                    matches[0], planet_index
                )
                record["extensions"].append(extension_summary)
                issues.extend(extension_issues)
    except (OSError, ValueError, TypeError) as error:
        issues.append(_issue("invalid_fits", f"Cannot parse result FITS: {error}"))
        return record

    record["ready"] = not issues and all(item["ready"] for item in record["extensions"])
    return record


def build_readiness(result_paths: Mapping[str, Path]) -> dict[str, Any]:
    validated = validate_target_mapping(result_paths)
    targets = [inspect_result(target, validated[target]) for target in PILOT_TARGETS]
    return {
        "schema_version": SCHEMA_VERSION,
        "stage": "timing_readiness_only",
        "pilot_targets": list(PILOT_TARGETS),
        "required_columns": list(REQUIRED_COLUMNS),
        "overall_ready": all(record["ready"] for record in targets),
        "ready_target_count": sum(bool(record["ready"]) for record in targets),
        "target_count": len(targets),
        "timing_inference_performed": False,
        "light_curves_downloaded": False,
        "transits_injected_or_fit": False,
        "targets": targets,
    }


def _target_csv_rows(payload: Mapping[str, Any]) -> Iterable[dict[str, Any]]:
    for record in payload["targets"]:
        reasons = record["reasons"]
        yield {
            "target": record["target"],
            "ready": record["ready"],
            "input_path": record["input_path"],
            "input_sha256": record["input_sha256"] or "",
            "input_bytes": record["input_bytes"] if record["input_bytes"] is not None else "",
            "fits_target": record["fits_target"] or "",
            "npl": record["npl"] if record["npl"] is not None else "",
            "validated_extension_count": sum(bool(item["ready"]) for item in record["extensions"]),
            "reason_codes": ";".join(item["code"] for item in reasons),
            "reason_messages": " | ".join(item["message"] for item in reasons),
        }


def write_readiness(payload: Mapping[str, Any], output_dir: Path) -> tuple[Path, Path]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ReadinessInputError(f"Output directory must be absent or empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "timing_readiness.json"
    csv_path = output_dir / "timing_readiness.csv"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rows = list(_target_csv_rows(payload))
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return json_path, csv_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate embedded ALDERAAN timings for the frozen six-target pilot."
    )
    parser.add_argument(
        "--result",
        action="append",
        default=[],
        metavar="TARGET=PATH",
        help=(
            "Result FITS mapping; repeat exactly once for each pilot target. "
            f"If omitted, supply a JSON object through {ENV_VAR}."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="New or empty local directory for readiness JSON and CSV only.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.result:
            if os.environ.get(ENV_VAR):
                raise ReadinessInputError(f"Use either --result or {ENV_VAR}, not both")
            result_paths = parse_result_specifications(args.result)
        else:
            raw_environment = os.environ.get(ENV_VAR)
            if not raw_environment:
                raise ReadinessInputError(
                    f"Provide six --result TARGET=PATH arguments or set {ENV_VAR}"
                )
            result_paths = parse_environment_mapping(raw_environment)
        payload = build_readiness(result_paths)
        json_path, csv_path = write_readiness(payload, args.output_dir)
    except ReadinessInputError as error:
        parser.error(str(error))

    print(f"Readiness: {payload['ready_target_count']}/{payload['target_count']} targets")
    print(f"JSON: {json_path}")
    print(f"CSV: {csv_path}")
    return 0 if payload["overall_ready"] else 2


if __name__ == "__main__":
    sys.exit(main())
