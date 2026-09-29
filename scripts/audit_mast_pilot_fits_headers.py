"""Audit MAST pilot FITS headers with bounded HTTP range requests only.

This utility does not download a light-curve payload, fit a transit, or judge
whether a product is suitable for a likelihood replay.  It only attests the
identity and basic cadence metadata that are available in the first 64 KiB.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import requests
from astropy.io import fits


SUMMARY_SHA256 = "020008a34eec7993affb8bb67ff553e32b9e4c5de413800e79b152d7273b9bae"
PILOT_TARGET_COUNT = 24
CADENCES = ("long", "short")
MAX_RANGE_BYTES = 65_536
RANGE_HEADER = "bytes=0-65535"
MAST_ROOT = "https://archive.stsci.edu/pub/kepler/lightcurves"
CONTENT_RANGE_RE = re.compile(r"^bytes\s+(\d+)-(\d+)/(\d+)$", re.IGNORECASE)
FILENAME_RE = re.compile(r"^kplr(\d{9})-[A-Za-z0-9]+_(llc|slc)\.fits$")


class AuditInputError(ValueError):
    """Raised when a bounded audit input does not satisfy its frozen contract."""


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def require_columns(rows: list[dict[str, str]], names: set[str], label: str) -> None:
    if not rows:
        raise AuditInputError(f"{label} has no rows")
    missing = names - set(rows[0])
    if missing:
        raise AuditInputError(f"{label} is missing columns: {sorted(missing)}")


def as_int(value: str, label: str) -> int:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise AuditInputError(f"Invalid integer for {label}: {value!r}") from exc


def is_true(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def read_and_validate_summary(path: Path, expected_hash: str = SUMMARY_SHA256) -> dict[str, int]:
    if sha256_path(path) != expected_hash:
        raise AuditInputError("Summary SHA-256 does not match the fixed source-faithful cohort")
    rows = read_csv(path)
    require_columns(
        rows,
        {
            "kepoi_name", "koi_target", "kepid", "qc_primary_exclude", "qc_reasons",
            "impact_mode", "nested_weight_mode", "include_transit_prior", "posterior_source",
        },
        "summary",
    )
    names = [row["kepoi_name"] for row in rows]
    if len(rows) != 2123 or len(set(names)) != 2123:
        raise AuditInputError("Fixed summary must contain exactly 2,123 unique planets")
    excluded = [row for row in rows if is_true(row["qc_primary_exclude"])]
    if len(excluded) != 4 or any("importance_ess_below_threshold" not in row["qc_reasons"] for row in excluded):
        raise AuditInputError("Fixed summary must retain exactly four documented ESS exclusions")
    kept = [row for row in rows if not is_true(row["qc_primary_exclude"])]
    if len(kept) != 2119 or len({row["koi_target"] for row in kept}) != 1583:
        raise AuditInputError("Fixed summary membership must be 2,119 planets in 1,583 targets")
    for row in kept:
        if row["impact_mode"] != "alderaan":
            raise AuditInputError(f"{row['kepoi_name']} does not use paired ALDERAAN impact samples")
        # The summary records the weighting method as ``dynesty``; the source
        # FITS weight column used by that method is LN_WT.
        if row["nested_weight_mode"] != "dynesty":
            raise AuditInputError(f"{row['kepoi_name']} does not retain Dynesty weighting")
        if is_true(row["include_transit_prior"]):
            raise AuditInputError(f"{row['kepoi_name']} has a transit prior already included")
    mapping: dict[str, int] = {}
    for row in kept:
        target, kepid = row["koi_target"], as_int(row["kepid"], f"kepid for {row['kepoi_name']}")
        if target in mapping and mapping[target] != kepid:
            raise AuditInputError(f"Summary has ambiguous KIC mapping for {target}")
        mapping[target] = kepid
    return mapping


def read_and_validate_pilot(path: Path, summary_kics: dict[str, int]) -> dict[str, dict[str, str]]:
    rows = read_csv(path)
    require_columns(rows, {"koi_target", "source_class", "candidate_role", "priority_planet"}, "pilot CSV")
    targets = [row["koi_target"] for row in rows]
    if len(rows) != PILOT_TARGET_COUNT or len(set(targets)) != PILOT_TARGET_COUNT:
        raise AuditInputError("Pilot CSV must contain exactly 24 unique predeclared targets")
    missing = set(targets) - set(summary_kics)
    if missing:
        raise AuditInputError(f"Pilot targets are absent from the fixed summary: {sorted(missing)}")
    for row in rows:
        if not row["source_class"] or not row["candidate_role"] or not row["priority_planet"]:
            raise AuditInputError(f"Pilot linkage is incomplete for {row['koi_target']}")
    return {row["koi_target"]: row for row in rows}


def read_and_validate_availability(
    path: Path, pilot: dict[str, dict[str, str]], summary_kics: dict[str, int]
) -> list[dict[str, str]]:
    rows = read_csv(path)
    require_columns(rows, {"koi_target", "kepid", "cadence", "query_status", "filenames"}, "availability CSV")
    pairs = [(row["koi_target"], row["cadence"]) for row in rows]
    expected = {(target, cadence) for target in pilot for cadence in CADENCES}
    if len(rows) != 48 or len(set(pairs)) != 48 or set(pairs) != expected:
        raise AuditInputError("Availability CSV must contain each of 24 pilot targets once for long and short cadence")
    for row in rows:
        target = row["koi_target"]
        if as_int(row["kepid"], f"availability kepid for {target}") != summary_kics[target]:
            raise AuditInputError(f"Availability KIC disagrees with fixed summary for {target}")
    return rows


def parse_content_range(value: str) -> tuple[int, int, int]:
    match = CONTENT_RANGE_RE.fullmatch((value or "").strip())
    if not match:
        raise ValueError("missing_or_malformed_content_range")
    start, end, total = (int(part) for part in match.groups())
    if start != 0 or end < start or end >= total or end >= MAX_RANGE_BYTES:
        raise ValueError("invalid_content_range_bounds")
    return start, end, total


def mast_url(kepid: int, filename: str, cadence: str) -> str:
    match = FILENAME_RE.fullmatch(filename)
    if not match:
        raise AuditInputError(f"Unsafe or invalid listed filename: {filename!r}")
    file_kepid, suffix = int(match.group(1)), match.group(2)
    if file_kepid != kepid:
        raise AuditInputError(f"Listed filename KIC does not match row KIC: {filename}")
    if suffix != ("llc" if cadence == "long" else "slc"):
        raise AuditInputError(f"Listed filename cadence does not match row cadence: {filename}")
    kid = f"{kepid:09d}"
    return f"{MAST_ROOT}/{kid[:4]}/{kid}/{filename}"


def base_record(row: dict[str, str], pilot_row: dict[str, str], filename: str = "") -> dict[str, Any]:
    return {
        "koi_target": row["koi_target"],
        "kepid": as_int(row["kepid"], f"availability kepid for {row['koi_target']}"),
        "source_class": pilot_row["source_class"],
        "candidate_role": pilot_row["candidate_role"],
        "priority_planet": pilot_row["priority_planet"],
        "cadence": row["cadence"],
        "filename": filename,
        "url": "",
        "http_status": None,
        "range_requested": RANGE_HEADER,
        "range_bytes_received": 0,
        "full_size_bytes": None,
        "primary_keplerid": None,
        "primary_obsmode": None,
        "primary_quarter": None,
        "hdu1_extname": None,
        "hdu1_tstart": None,
        "hdu1_tstop": None,
        "hdu1_timedel": None,
        "hdu1_int_time": None,
        "hdu1_num_frm": None,
        "hdu1_ttype_columns": "",
        "identity_consistent": False,
        "cadence_consistent": False,
        "status": "unstarted",
        "error": "",
    }


def safe_header_value(header: Any, key: str) -> Any:
    value = header.get(key)
    return value.item() if hasattr(value, "item") else value


def inspect_fits_headers(data: bytes, record: dict[str, Any]) -> dict[str, Any]:
    try:
        with fits.open(io.BytesIO(data), memmap=False, lazy_load_hdus=True, ignore_missing_end=True) as hdul:
            primary = hdul[0].header
            record["primary_keplerid"] = safe_header_value(primary, "KEPLERID")
            record["primary_obsmode"] = safe_header_value(primary, "OBSMODE")
            record["primary_quarter"] = safe_header_value(primary, "QUARTER")
            if len(hdul) > 1:
                header = hdul[1].header
                record["hdu1_extname"] = safe_header_value(header, "EXTNAME")
                record["hdu1_tstart"] = safe_header_value(header, "TSTART")
                record["hdu1_tstop"] = safe_header_value(header, "TSTOP")
                record["hdu1_timedel"] = safe_header_value(header, "TIMEDEL")
                record["hdu1_int_time"] = safe_header_value(header, "INT_TIME")
                record["hdu1_num_frm"] = safe_header_value(header, "NUM_FRM")
                names = getattr(getattr(hdul[1], "columns", None), "names", None) or []
                record["hdu1_ttype_columns"] = "|".join(str(name) for name in names)
    except Exception as exc:  # astropy errors vary by truncated header arrangement.
        record["status"] = "fits_header_parse_failed"
        record["error"] = f"{type(exc).__name__}:{exc}"
        return record

    errors: list[str] = []
    try:
        record["identity_consistent"] = int(record["primary_keplerid"]) == record["kepid"]
    except (TypeError, ValueError):
        errors.append("missing_or_invalid_primary_KEPLERID")
    if not record["identity_consistent"] and not errors:
        errors.append("primary_KEPLERID_mismatch")
    mode = str(record["primary_obsmode"] or "").lower()
    record["cadence_consistent"] = record["cadence"] in mode
    if not mode:
        errors.append("missing_primary_OBSMODE")
    elif not record["cadence_consistent"]:
        errors.append("primary_OBSMODE_cadence_mismatch")
    if record["hdu1_extname"] is None:
        errors.append("missing_HDU1_header")
    if errors:
        record["status"] = "header_incomplete"
        record["error"] = ";".join(errors)
    else:
        record["status"] = "ok"
    return record


def fetch_product(record: dict[str, Any], timeout: float, requester=requests.get) -> dict[str, Any]:
    try:
        response = requester(record["url"], headers={"Range": RANGE_HEADER}, stream=True, timeout=timeout)
    except Exception as exc:
        record["status"] = "request_failed"
        record["error"] = f"{type(exc).__name__}:{exc}"
        return record
    try:
        record["http_status"] = response.status_code
        if response.status_code != 206:
            record["status"] = "range_not_honored" if response.status_code == 200 else "unexpected_http_status"
            record["error"] = f"HTTP {response.status_code}; expected 206 Partial Content"
            return record
        try:
            start, end, total = parse_content_range(response.headers.get("Content-Range", ""))
        except ValueError as exc:
            record["status"] = "invalid_content_range"
            record["error"] = str(exc)
            return record
        expected_length = end - start + 1
        content_length = response.headers.get("Content-Length")
        if content_length is not None and int(content_length) != expected_length:
            record["status"] = "invalid_content_range"
            record["error"] = "Content-Length disagrees with Content-Range"
            return record
        chunks: list[bytes] = []
        received = 0
        for chunk in response.iter_content(chunk_size=8192):
            if not chunk:
                continue
            if received + len(chunk) > MAX_RANGE_BYTES:
                record["status"] = "range_byte_cap_exceeded"
                record["error"] = "response exceeded 65,536-byte cap"
                return record
            chunks.append(chunk)
            received += len(chunk)
        record["range_bytes_received"] = received
        record["full_size_bytes"] = total
        if received != expected_length:
            record["status"] = "range_length_mismatch"
            record["error"] = f"received {received} bytes; Content-Range declared {expected_length}"
            return record
        return inspect_fits_headers(b"".join(chunks), record)
    except Exception as exc:
        record["status"] = "request_read_failed"
        record["error"] = f"{type(exc).__name__}:{exc}"
        return record
    finally:
        response.close()


def product_jobs(
    availability: list[dict[str, str]], pilot: dict[str, dict[str, str]]
) -> list[dict[str, Any]]:
    jobs: list[dict[str, Any]] = []
    for row in availability:
        base = base_record(row, pilot[row["koi_target"]])
        if row["query_status"] != "success":
            base["status"] = "listing_not_success"
            base["error"] = row.get("query_error", "") or f"listing status {row['query_status']}"
            jobs.append(base)
            continue
        names = [name for name in row["filenames"].split("|") if name]
        if not names:
            base["status"] = "no_listed_products"
            jobs.append(base)
            continue
        for name in names:
            record = base_record(row, pilot[row["koi_target"]], name)
            record["url"] = mast_url(record["kepid"], name, record["cadence"])
            jobs.append(record)
    return jobs


def audit_products(jobs: list[dict[str, Any]], workers: int, timeout: float) -> list[dict[str, Any]]:
    no_request = [job for job in jobs if job["status"] != "unstarted"]
    requested = [job for job in jobs if job["status"] == "unstarted"]
    if workers == 1:
        completed = [fetch_product(job, timeout) for job in requested]
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            completed = list(executor.map(lambda job: fetch_product(job, timeout), requested))
    return no_request + completed


def write_outputs(records: list[dict[str, Any]], args: argparse.Namespace) -> None:
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    if args.output_csv.exists() or args.output_json.exists():
        raise AuditInputError("Refusing to overwrite an existing explicit audit output")
    fields = list(records[0]) if records else list(base_record({"koi_target":"", "kepid":"0", "cadence":"long"}, {"source_class":"", "candidate_role":"", "priority_planet":""}))
    with args.output_csv.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    by_cadence: dict[str, int] = {}
    by_source: dict[str, int] = {}
    by_status: dict[str, int] = {}
    for record in records:
        by_cadence[record["cadence"]] = by_cadence.get(record["cadence"], 0) + 1
        by_source[record["source_class"]] = by_source.get(record["source_class"], 0) + 1
        by_status[record["status"]] = by_status.get(record["status"], 0) + 1
    payload = {
        "audit": "MAST FITS header audit with Range bytes=0-65535 only",
        "availability_csv": str(args.availability_csv.resolve()),
        "summary": str(args.summary.resolve()),
        "summary_sha256": sha256_path(args.summary),
        "pilot_csv": str(args.pilot_csv.resolve()),
        "workers": args.workers,
        "timeout_seconds": args.timeout,
        "product_records": len(records),
        "successful_headers": by_status.get("ok", 0),
        "transfer_bytes": sum(record["range_bytes_received"] for record in records),
        "counts_by_cadence": by_cadence,
        "counts_by_source": by_source,
        "counts_by_status": by_status,
        "explicit_limits": [
            "Header presence does not verify finite flux or flux uncertainty.",
            "Header presence does not verify quality flags, transit coverage, or LC/SC overlap.",
            "This audit does not attest transit timings, native masks, likelihood conventions, or replay readiness.",
            "No FITS payload beyond the requested first 65,536 bytes was intentionally read.",
            "This audit makes no fit requests and does not alter canonical data.",
        ],
    }
    with args.output_json.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--availability-csv", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--pilot-csv", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=30.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not 1 <= args.workers <= 4:
        raise AuditInputError("--workers must be between 1 and 4")
    if args.timeout <= 0 or args.timeout > 120:
        raise AuditInputError("--timeout must be greater than zero and no more than 120 seconds")
    summary_kics = read_and_validate_summary(args.summary)
    pilot = read_and_validate_pilot(args.pilot_csv, summary_kics)
    availability = read_and_validate_availability(args.availability_csv, pilot, summary_kics)
    records = audit_products(product_jobs(availability, pilot), args.workers, args.timeout)
    write_outputs(records, args)
    print(f"Wrote {len(records)} header audit records to {args.output_csv}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AuditInputError as exc:
        print(f"input error: {exc}", file=sys.stderr)
        raise SystemExit(2)
