"""Create a bounded, provenance-preserving circular transit injection.

This is a new diagnostic implementation.  It is not Sagear's unavailable
``alderaan.validate`` module and it does not alter ALDERAAN results.

The input NPZ must contain finite one-dimensional ``time`` and ``flux`` arrays
of equal length plus an ordered integer ``mask_index`` array.  A transit is
applied only at ``mask_index``; all other flux values are retained exactly.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import importlib
import json
from pathlib import Path
import platform
import sys
from typing import Any

import numpy as np


G_SI = 6.67430e-11
RHO_SUN_KG_M3 = 1408.0
MANIFEST_NAME = "diagnostic_injection_manifest.json"
OUTPUT_NAME = "diagnostic_injection.npz"
POINT_ALIGNED_OPTIONAL_ARRAYS = (
    "flux_err",
    "cadence",
    "quality",
    "quarter",
    "source_file_index",
    "source_row_index",
    "known_transit_mask",
)
SOURCE_MAPPING_ARRAY = "source_files"
IN_TRANSIT_TOLERANCE = 1.0e-12


class InjectionInputError(ValueError):
    """Raised when an injection specification or input contract is invalid."""


class BatmanDependencyError(RuntimeError):
    """Raised when rendering is requested without the optional batman package."""


@dataclass(frozen=True)
class InjectionSpec:
    """Parameters for a circular, exposure-integrated diagnostic transit."""

    period_days: float
    rho_star_solar: float
    radius_ratio: float
    impact_parameter: float
    t0_days: float
    limb_darkening_u1: float
    limb_darkening_u2: float
    exposure_minutes: float
    supersample_factor: int
    seed: int

    def validate(self) -> None:
        finite_positive = {
            "period_days": self.period_days,
            "rho_star_solar": self.rho_star_solar,
            "radius_ratio": self.radius_ratio,
            "exposure_minutes": self.exposure_minutes,
        }
        for name, value in finite_positive.items():
            if not np.isfinite(value) or value <= 0:
                raise InjectionInputError(f"{name} must be finite and positive")
        if self.radius_ratio >= 1.0:
            raise InjectionInputError("radius_ratio must be smaller than one")
        if not np.isfinite(self.impact_parameter) or self.impact_parameter < 0:
            raise InjectionInputError("impact_parameter must be finite and non-negative")
        if self.impact_parameter > 1.0 + self.radius_ratio:
            raise InjectionInputError("impact_parameter exceeds the transiting limit 1 + radius_ratio")
        if not np.isfinite(self.t0_days):
            raise InjectionInputError("t0_days must be finite")
        if not np.isfinite(self.limb_darkening_u1) or not np.isfinite(self.limb_darkening_u2):
            raise InjectionInputError("limb-darkening coefficients must be finite")
        if self.limb_darkening_u1 < 0 or self.limb_darkening_u1 + self.limb_darkening_u2 > 1:
            raise InjectionInputError("limb-darkening coefficients violate the quadratic intensity bounds")
        if self.limb_darkening_u1 + 2.0 * self.limb_darkening_u2 < 0:
            raise InjectionInputError("limb-darkening coefficients violate the quadratic intensity bounds")
        if int(self.supersample_factor) != self.supersample_factor or self.supersample_factor < 1:
            raise InjectionInputError("supersample_factor must be an integer of at least one")
        if int(self.seed) != self.seed or self.seed < 0:
            raise InjectionInputError("seed must be a non-negative integer")
        if self.impact_parameter >= circular_a_over_rstar(self.period_days, self.rho_star_solar):
            raise InjectionInputError("impact_parameter must be smaller than circular a_over_rstar")


def circular_a_over_rstar(period_days: float, rho_star_solar: float) -> float:
    """Return circular a/Rstar from period and stellar density in solar units."""
    if not np.isfinite(period_days) or period_days <= 0:
        raise InjectionInputError("period_days must be finite and positive")
    if not np.isfinite(rho_star_solar) or rho_star_solar <= 0:
        raise InjectionInputError("rho_star_solar must be finite and positive")
    period_seconds = period_days * 86400.0
    return float((G_SI * RHO_SUN_KG_M3 * rho_star_solar * period_seconds**2 / (3.0 * np.pi)) ** (1.0 / 3.0))


def circular_rho_star_solar(period_days: float, a_over_rstar: float) -> float:
    """Return density in solar units from a circular period and a/Rstar."""
    if not np.isfinite(period_days) or period_days <= 0:
        raise InjectionInputError("period_days must be finite and positive")
    if not np.isfinite(a_over_rstar) or a_over_rstar <= 0:
        raise InjectionInputError("a_over_rstar must be finite and positive")
    period_seconds = period_days * 86400.0
    return float(3.0 * np.pi * a_over_rstar**3 / (G_SI * period_seconds**2 * RHO_SUN_KG_M3))


def circular_first_to_fourth_duration_days(
    period_days: float,
    rho_star_solar: float,
    radius_ratio: float,
    impact_parameter: float,
) -> float:
    """Return the exact circular first-to-fourth-contact duration in days."""
    a_over_rstar = circular_a_over_rstar(period_days, rho_star_solar)
    if not np.isfinite(radius_ratio) or radius_ratio <= 0.0 or radius_ratio >= 1.0:
        raise InjectionInputError("radius_ratio must be finite, positive, and smaller than one")
    if not np.isfinite(impact_parameter) or impact_parameter < 0.0:
        raise InjectionInputError("impact_parameter must be finite and non-negative")
    if impact_parameter >= a_over_rstar:
        raise InjectionInputError("impact_parameter must be smaller than circular a_over_rstar")
    chord_squared = (1.0 + radius_ratio) ** 2 - impact_parameter**2
    if chord_squared <= 0.0:
        raise InjectionInputError("circular geometry does not yield a transiting first-to-fourth chord")
    sin_inclination = np.sqrt(1.0 - (impact_parameter / a_over_rstar) ** 2)
    argument = np.sqrt(chord_squared) / (a_over_rstar * sin_inclination)
    if not np.isfinite(argument) or not 0.0 < argument < 1.0:
        raise InjectionInputError("circular geometry does not yield a finite first-to-fourth duration")
    return float(period_days / np.pi * np.arcsin(argument))


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_array(values: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(values)
    return hashlib.sha256(contiguous.tobytes()).hexdigest()


def _array_report(values: np.ndarray) -> dict[str, Any]:
    """Describe an array without converting its dtype or values."""
    return {
        "dtype": values.dtype.str,
        "shape": [int(size) for size in values.shape],
        "sha256": _sha256_array(values),
    }


def load_input_contract(path: Path) -> dict[str, np.ndarray]:
    """Load and validate an NPZ light curve without changing its selected rows."""
    if not path.is_file():
        raise InjectionInputError(f"input NPZ does not exist: {path}")
    try:
        with np.load(path, allow_pickle=False) as payload:
            missing = {"time", "flux", "mask_index"} - set(payload.files)
            if missing:
                raise InjectionInputError("input NPZ is missing " + ", ".join(sorted(missing)))
            time = np.asarray(payload["time"], dtype=np.float64)
            flux = np.asarray(payload["flux"], dtype=np.float64)
            raw_index = np.asarray(payload["mask_index"])
            optional = {
                name: np.asarray(payload[name])
                for name in (*POINT_ALIGNED_OPTIONAL_ARRAYS, SOURCE_MAPPING_ARRAY)
                if name in payload.files
            }
    except (OSError, ValueError) as exc:
        if isinstance(exc, InjectionInputError):
            raise
        raise InjectionInputError(f"could not read input NPZ: {exc}") from exc
    if time.ndim != 1 or flux.ndim != 1 or len(time) != len(flux) or len(time) == 0:
        raise InjectionInputError("time and flux must be non-empty one-dimensional arrays of equal length")
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(flux)):
        raise InjectionInputError("time and flux must be finite")
    if np.any(np.diff(time) <= 0):
        raise InjectionInputError("time must be strictly increasing")
    if raw_index.ndim != 1 or len(raw_index) == 0:
        raise InjectionInputError("mask_index must be a non-empty one-dimensional array")
    if raw_index.dtype.kind not in "iu":
        raise InjectionInputError("mask_index must have an integer dtype")
    mask_index = raw_index.astype(np.int64, copy=False)
    if mask_index[0] < 0 or mask_index[-1] >= len(time):
        raise InjectionInputError("mask_index is outside the time array")
    if np.any(np.diff(mask_index) <= 0):
        raise InjectionInputError("mask_index must be strictly increasing and unique")
    for name in POINT_ALIGNED_OPTIONAL_ARRAYS:
        if name in optional and (optional[name].ndim == 0 or len(optional[name]) != len(time)):
            raise InjectionInputError(f"optional array {name} must have length matching time")
    if SOURCE_MAPPING_ARRAY in optional and optional[SOURCE_MAPPING_ARRAY].ndim == 0:
        raise InjectionInputError("optional array source_files must be an array source mapping")
    return {"time": time, "flux": flux, "mask_index": mask_index, **optional}


def _batman_module() -> Any:
    try:
        return importlib.import_module("batman")
    except ModuleNotFoundError as exc:
        raise BatmanDependencyError(
            "batman-package is required to render an injection. Install batman-package "
            "or use --dry-run for contract validation only."
        ) from exc


def _supports_parameter_contract(params: Any, *names: str) -> bool:
    """Identify declared batman parameters without assigning unsupported fields."""
    slots = getattr(type(params), "__slots__", ())
    if isinstance(slots, str):
        slots = (slots,)
    declared = set(slots)
    return all(hasattr(params, name) or name in declared for name in names)


def exposure_integrated_circular_model(spec: InjectionSpec, time: np.ndarray) -> tuple[np.ndarray, str]:
    """Render a circular transit at selected times with batman's exposure integration."""
    spec.validate()
    batman = _batman_module()
    params = batman.TransitParams()
    a_over_rstar = circular_a_over_rstar(spec.period_days, spec.rho_star_solar)
    params.t0 = spec.t0_days
    params.per = spec.period_days
    params.rp = spec.radius_ratio
    if _supports_parameter_contract(params, "a", "inc"):
        params.a = a_over_rstar
        params.inc = float(np.degrees(np.arccos(spec.impact_parameter / a_over_rstar)))
        params.ecc = 0.0
        params.w = 90.0
    elif _supports_parameter_contract(params, "b", "T14"):
        params.b = spec.impact_parameter
        params.T14 = circular_first_to_fourth_duration_days(
            spec.period_days,
            spec.rho_star_solar,
            spec.radius_ratio,
            spec.impact_parameter,
        )
    else:
        raise BatmanDependencyError(
            "batman's TransitParams supports neither the public a/inc contract nor "
            "the ALDERAAN b/T14 contract"
        )
    params.u = [spec.limb_darkening_u1, spec.limb_darkening_u2]
    params.limb_dark = "quadratic"
    model = batman.TransitModel(
        params,
        np.asarray(time, dtype=np.float64),
        supersample_factor=int(spec.supersample_factor),
        exp_time=spec.exposure_minutes / 1440.0,
    )
    return np.asarray(model.light_curve(params), dtype=np.float64), str(getattr(batman, "__version__", "unknown"))


def evaluate_injection_coverage(full_model: np.ndarray, mask_index: np.ndarray) -> dict[str, Any]:
    """Verify that the preserved injection mask retains every modeled transit cadence."""
    model = np.asarray(full_model, dtype=np.float64)
    if model.ndim != 1 or not np.all(np.isfinite(model)):
        raise RuntimeError("batman returned an invalid model for the full time grid")
    eligible = np.zeros(len(model), dtype=bool)
    eligible[np.asarray(mask_index, dtype=np.int64)] = True
    in_transit = model < 1.0 - IN_TRANSIT_TOLERANCE
    total = int(np.count_nonzero(in_transit))
    retained = int(np.count_nonzero(in_transit & eligible))
    excluded = total - retained
    coverage = {
        "full_grid_count": int(len(model)),
        "eligible_grid_count": int(np.count_nonzero(eligible)),
        "in_transit_full_grid_count": total,
        "eligible_in_transit_count": retained,
        "masked_in_transit_count": excluded,
        "in_transit_tolerance": IN_TRANSIT_TOLERANCE,
        "all_in_transit_cadences_eligible": excluded == 0,
    }
    if total == 0:
        raise InjectionInputError("the supplied time grid does not contain an injected transit cadence")
    if excluded:
        raise InjectionInputError(
            "the preserved mask would remove "
            f"{excluded} of {total} modeled in-transit cadences; choose a non-overlapping injection epoch"
        )
    return coverage


def build_manifest(
    spec: InjectionSpec,
    input_path: Path,
    contract: dict[str, np.ndarray],
    *,
    dry_run: bool,
    batman_version: str | None,
    coverage: dict[str, Any] | None = None,
    output_path: Path | None = None,
) -> dict[str, Any]:
    """Build deterministic provenance for the exact declared input and geometry."""
    a_over_rstar = circular_a_over_rstar(spec.period_days, spec.rho_star_solar)
    optional_reports = {
        name: _array_report(contract[name])
        for name in (*POINT_ALIGNED_OPTIONAL_ARRAYS, SOURCE_MAPPING_ARRAY)
        if name in contract
    }
    output = None if output_path is None else {
        "filename": output_path.name,
        "sha256": _sha256_file(output_path),
        "optional_arrays": optional_reports,
    }
    return {
        "schema_version": 1,
        "scope": "New standalone diagnostic injection; not Sagear's missing alderaan.validate module and not an ALDERAAN result.",
        "mode": "preflight" if dry_run else "rendered_injection",
        "specification": asdict(spec),
        "derived_circular_geometry": {
            "eccentricity": 0.0,
            "a_over_rstar": a_over_rstar,
            "rho_star_solar_round_trip": circular_rho_star_solar(spec.period_days, a_over_rstar),
        },
        "input": {
            "filename": input_path.name,
            "sha256": _sha256_file(input_path),
            "point_count": int(len(contract["time"])),
            "mask_count": int(len(contract["mask_index"])),
            "time_sha256": _sha256_array(contract["time"]),
            "flux_sha256": _sha256_array(contract["flux"]),
            "mask_index_sha256": _sha256_array(contract["mask_index"]),
            "optional_arrays": optional_reports,
        },
        "output": output,
        "injected_transit_coverage": coverage,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "batman": batman_version if batman_version is not None else "not_loaded_preflight",
        },
        "seed": int(spec.seed),
        "random_noise_added": False,
    }


def _write_json_once(path: Path, payload: dict[str, Any]) -> None:
    content = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    try:
        with path.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
    except FileExistsError as exc:
        raise InjectionInputError(f"refusing to replace existing immutable manifest: {path}") from exc


def run_injection(spec: InjectionSpec, input_path: Path, output_dir: Path, *, dry_run: bool = False) -> Path:
    """Validate, render when requested, and write a manifest into an empty directory."""
    spec.validate()
    contract = load_input_contract(input_path)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise InjectionInputError("output_dir must be absent or empty to preserve immutable provenance")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / MANIFEST_NAME
    if dry_run:
        _write_json_once(manifest_path, build_manifest(spec, input_path, contract, dry_run=True, batman_version=None))
        return manifest_path
    full_model, version = exposure_integrated_circular_model(spec, contract["time"])
    if full_model.shape != contract["time"].shape:
        raise RuntimeError("batman returned an invalid model for the full time grid")
    coverage = evaluate_injection_coverage(full_model, contract["mask_index"])
    model = full_model[contract["mask_index"]]
    injected_flux = contract["flux"].copy()
    injected_flux[contract["mask_index"]] *= model
    output_path = output_dir / OUTPUT_NAME
    output_payload = {
        "time": contract["time"],
        "flux": injected_flux,
        "mask_index": contract["mask_index"],
        "transit_model": model,
    }
    output_payload.update({
        name: contract[name]
        for name in (*POINT_ALIGNED_OPTIONAL_ARRAYS, SOURCE_MAPPING_ARRAY)
        if name in contract
    })
    np.savez_compressed(output_path, **output_payload)
    _write_json_once(
        manifest_path,
        build_manifest(
            spec,
            input_path,
            contract,
            dry_run=False,
            batman_version=version,
            coverage=coverage,
            output_path=output_path,
        ),
    )
    return manifest_path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-npz", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--period-days", type=float, required=True)
    parser.add_argument("--rho-star-solar", type=float, required=True)
    parser.add_argument("--radius-ratio", type=float, required=True)
    parser.add_argument("--impact-parameter", type=float, required=True)
    parser.add_argument("--t0-days", type=float, required=True)
    parser.add_argument("--limb-darkening", type=float, nargs=2, metavar=("U1", "U2"), required=True)
    parser.add_argument("--exposure-minutes", type=float, required=True)
    parser.add_argument("--supersample-factor", type=int, default=15)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs and write only immutable provenance.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    spec = InjectionSpec(
        period_days=args.period_days,
        rho_star_solar=args.rho_star_solar,
        radius_ratio=args.radius_ratio,
        impact_parameter=args.impact_parameter,
        t0_days=args.t0_days,
        limb_darkening_u1=args.limb_darkening[0],
        limb_darkening_u2=args.limb_darkening[1],
        exposure_minutes=args.exposure_minutes,
        supersample_factor=args.supersample_factor,
        seed=args.seed,
    )
    manifest = run_injection(spec, args.input_npz, args.output_dir, dry_run=args.dry_run)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
