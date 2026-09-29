"""Build an isolated MAST-FITS fixture for one circular ALDERAAN diagnostic.

This is a new, single-planet control.  It writes no result FITS and does not
alter source PDCSAP files, catalog files, or canonical eccentricity products.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits

from diagnostic_transit_injection import circular_a_over_rstar


class FixtureInputError(ValueError):
    """Raised when an injection cannot become a trustworthy ALDERAAN fixture."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise FixtureInputError(f"could not read JSON manifest: {path}") from error


def _single_planet_catalog(catalog_path: Path, target: str) -> pd.DataFrame:
    try:
        catalog = pd.read_csv(catalog_path, index_col=0)
    except (OSError, pd.errors.ParserError) as error:
        raise FixtureInputError(f"could not read catalog: {catalog_path}") from error
    required = {"koi_id", "kic_id", "npl", "period", "epoch", "depth", "duration", "impact", "limbdark_1", "limbdark_2"}
    missing = required - set(catalog.columns)
    if missing:
        raise FixtureInputError("catalog is missing " + ", ".join(sorted(missing)))
    rows = catalog.loc[catalog.koi_id.astype(str) == target].copy()
    if len(rows) != 1 or int(rows.iloc[0].npl) != 1:
        raise FixtureInputError("diagnostic fixture currently requires exactly one catalog row and npl=1")
    return rows


def _duration_hours(period_days: float, rho_star_solar: float, radius_ratio: float, impact: float) -> float:
    a_over_rstar = circular_a_over_rstar(period_days, rho_star_solar)
    if impact >= a_over_rstar:
        raise FixtureInputError("impact parameter is outside the circular orbit")
    sin_inclination = np.sqrt(1.0 - (impact / a_over_rstar) ** 2)
    argument = np.sqrt((1.0 + radius_ratio) ** 2 - impact**2) / (a_over_rstar * sin_inclination)
    if not 0.0 < argument < 1.0:
        raise FixtureInputError("circular geometry does not yield a finite first-to-fourth duration")
    return float(period_days * 24.0 / np.pi * np.arcsin(argument))


def _load_contract(prepared_path: Path, injection_path: Path) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    try:
        with np.load(prepared_path, allow_pickle=False) as prepared_payload:
            prepared = {name: np.asarray(prepared_payload[name]) for name in prepared_payload.files}
        with np.load(injection_path, allow_pickle=False) as injection_payload:
            injection = {name: np.asarray(injection_payload[name]) for name in injection_payload.files}
    except (OSError, ValueError) as error:
        raise FixtureInputError("could not read prepared or injected NPZ") from error
    required = {"time", "flux", "source_file_index", "source_row_index", "known_transit_mask", "source_files"}
    missing = required - set(prepared)
    if missing:
        raise FixtureInputError("prepared input is missing " + ", ".join(sorted(missing)))
    if set(("time", "flux")) - set(injection):
        raise FixtureInputError("injection output must include time and flux")
    n = len(prepared["time"])
    if n == 0 or any(np.asarray(prepared[name]).ndim != 1 or len(prepared[name]) != n for name in required - {"source_files"}):
        raise FixtureInputError("prepared point arrays must be non-empty one-dimensional arrays of equal length")
    if prepared["source_files"].ndim != 1 or len(prepared["source_files"]) == 0:
        raise FixtureInputError("prepared source_files must be a non-empty one-dimensional mapping")
    if len(injection["time"]) != n or len(injection["flux"]) != n:
        raise FixtureInputError("injection arrays do not match prepared input length")
    if not np.array_equal(np.asarray(prepared["time"]), np.asarray(injection["time"])):
        raise FixtureInputError("injection time array does not match prepared input")
    if not np.all(np.isfinite(injection["flux"])):
        raise FixtureInputError("injected flux must be finite before known-transit masking")
    return prepared, injection


def _validate_manifests(
    prepared_manifest: Path, injection_manifest: Path, prepared_path: Path, injection_path: Path
) -> tuple[dict[str, Any], dict[str, Any]]:
    prepared = _read_json(prepared_manifest)
    injection = _read_json(injection_manifest)
    if prepared.get("output", {}).get("sha256") != _sha256(prepared_path):
        raise FixtureInputError("prepared NPZ SHA-256 does not match its manifest")
    if injection.get("output", {}).get("sha256") != _sha256(injection_path):
        raise FixtureInputError("injection NPZ SHA-256 does not match its manifest")
    if injection.get("input", {}).get("sha256") != _sha256(prepared_path):
        raise FixtureInputError("injection manifest does not identify this prepared NPZ as its input")
    specification = injection.get("specification")
    if not isinstance(specification, dict):
        raise FixtureInputError("injection manifest lacks a specification")
    cadence = prepared.get("cadence")
    if not isinstance(cadence, dict):
        raise FixtureInputError("prepared manifest lacks cadence provenance")
    try:
        tolerance = float(cadence["exposure_match_tolerance_minutes"])
        prepared_exposure = float(cadence["exposure_minutes"])
        injected_exposure = float(specification["exposure_minutes"])
    except (KeyError, TypeError, ValueError) as error:
        raise FixtureInputError("prepared or injection manifest lacks a valid exposure declaration") from error
    if abs(prepared_exposure - injected_exposure) > tolerance:
        raise FixtureInputError("injection exposure does not match prepared cadence exposure")
    return specification, prepared


def build_fixture(
    *, prepared_path: Path, prepared_manifest: Path, injection_path: Path, injection_manifest: Path,
    catalog_path: Path, target: str, output_dir: Path,
) -> Path:
    """Write cloned PDCSAP files and a matching single-planet ALDERAAN catalog."""
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FixtureInputError("output directory must be absent or empty")
    prepared, injection = _load_contract(prepared_path, injection_path)
    spec, prepared_record = _validate_manifests(prepared_manifest, injection_manifest, prepared_path, injection_path)
    catalog_rows = _single_planet_catalog(catalog_path, target)
    if prepared_record.get("target") != target or int(prepared_record.get("kic", -1)) != int(catalog_rows.iloc[0].kic_id):
        raise FixtureInputError("prepared manifest target or KIC does not agree with the fixture catalog")
    source_paths = [Path(str(item)) for item in prepared["source_files"]]
    source_index = np.asarray(prepared["source_file_index"], dtype=np.int64)
    source_row = np.asarray(prepared["source_row_index"], dtype=np.int64)
    if np.any(source_index < 0) or np.any(source_index >= len(source_paths)) or np.any(source_row < 0):
        raise FixtureInputError("prepared source mapping contains an invalid row")
    source_pairs = np.column_stack((source_index, source_row))
    if len(np.unique(source_pairs, axis=0)) != len(source_pairs):
        raise FixtureInputError("prepared source mapping repeats a source file row")
    source_records = prepared_record.get("source_files")
    if not isinstance(source_records, list) or len(source_records) != len(source_paths):
        raise FixtureInputError("prepared manifest source-file records do not match the input mapping")
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = output_dir / "Data"
    data_dir.mkdir()
    output_files: list[dict[str, Any]] = []
    for file_index, source_path in enumerate(source_paths):
        if not source_path.is_file():
            raise FixtureInputError(f"source PDCSAP file does not exist: {source_path}")
        source_record = source_records[file_index]
        if Path(str(source_record.get("path", ""))).name != source_path.name:
            raise FixtureInputError(f"prepared manifest source filename does not match {source_path.name}")
        if source_record.get("sha256") != _sha256(source_path):
            raise FixtureInputError(f"source PDCSAP SHA-256 does not match the prepared manifest: {source_path.name}")
        use = np.flatnonzero(source_index == file_index)
        destination = data_dir / source_path.name
        with fits.open(source_path, memmap=False, checksum=False) as source_hdul:
            table = source_hdul[1].data
            for column in ("TIME", "PDCSAP_FLUX"):
                if column not in table.names:
                    raise FixtureInputError(f"{source_path.name} lacks {column}")
            if len(use) and int(np.max(source_row[use])) >= len(table):
                raise FixtureInputError(f"prepared source rows exceed {source_path.name}")
            original_flux = np.asarray(table["PDCSAP_FLUX"], dtype=np.float64)
            if len(use) and not np.allclose(original_flux[source_row[use]], prepared["flux"][use], rtol=0.0, atol=0.0, equal_nan=True):
                raise FixtureInputError(f"prepared flux does not reproduce {source_path.name}")
            if len(use) and not np.allclose(np.asarray(table["TIME"], dtype=np.float64)[source_row[use]], prepared["time"][use], rtol=0.0, atol=0.0):
                raise FixtureInputError(f"prepared times do not reproduce {source_path.name}")
            for row in use:
                if bool(prepared["known_transit_mask"][row]):
                    table["PDCSAP_FLUX"][source_row[row]] = np.nan
                else:
                    table["PDCSAP_FLUX"][source_row[row]] = injection["flux"][row]
            source_hdul[0].header["HILDDIAG"] = (True, "Synthetic circular-transit diagnostic")
            source_hdul.writeto(destination, overwrite=False, checksum=True)
        output_files.append({"source": str(source_path), "source_sha256": _sha256(source_path), "output": destination.name, "output_sha256": _sha256(destination), "prepared_rows": int(len(use)), "known_transit_rows_masked": int(np.sum(prepared["known_transit_mask"][use]))})
    row = catalog_rows.iloc[0].copy()
    row["period"] = float(spec["period_days"])
    row["epoch"] = float(spec["t0_days"])
    row["depth"] = float(spec["radius_ratio"]) ** 2 * 1e6
    row["duration"] = _duration_hours(float(spec["period_days"]), float(spec["rho_star_solar"]), float(spec["radius_ratio"]), float(spec["impact_parameter"]))
    row["impact"] = float(spec["impact_parameter"])
    row["limbdark_1"] = float(spec["limb_darkening_u1"])
    row["limbdark_2"] = float(spec["limb_darkening_u2"])
    fixture_catalog = output_dir / "diagnostic_catalog.csv"
    pd.DataFrame([row]).to_csv(fixture_catalog, index=True)
    manifest = {
        "schema_version": 1,
        "scope": "New single-planet ALDERAAN diagnostic fixture, not a Sagear result.",
        "target": target,
        "prepared_input": {"path": prepared_path.name, "sha256": _sha256(prepared_path)},
        "injection": {"path": injection_path.name, "sha256": _sha256(injection_path), "specification": spec},
        "catalog": {"source": str(catalog_path), "source_sha256": _sha256(catalog_path), "output": fixture_catalog.name, "output_sha256": _sha256(fixture_catalog)},
        "files": output_files,
        "known_transits_are_nan_masked": True,
        "results_written": False,
    }
    manifest_path = output_dir / "diagnostic_fixture_manifest.json"
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest_path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-input", type=Path, required=True)
    parser.add_argument("--prepared-manifest", type=Path, required=True)
    parser.add_argument("--injection", type=Path, required=True)
    parser.add_argument("--injection-manifest", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    print(build_fixture(prepared_path=args.prepared_input, prepared_manifest=args.prepared_manifest, injection_path=args.injection, injection_manifest=args.injection_manifest, catalog_path=args.catalog, target=args.target, output_dir=args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
