"""Create an immutable baseline for a controlled transit-injection diagnostic.

The input is a prepared raw-PDCSAP contract from
``prepare_diagnostic_injection_input.py``. The output retains the cadence
contract and can be supplied to ``diagnostic_transit_injection.py``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
from typing import Any

import numpy as np

from diagnostic_transit_injection import InjectionInputError, load_input_contract


MANIFEST_NAME = "diagnostic_background_manifest.json"
OUTPUT_NAME = "diagnostic_background.npz"
MODES = ("observed_pdcsap", "quoted_error_gaussian", "variance_matched_iid")


class BackgroundInputError(ValueError):
    """Raised when a requested diagnostic background is not well-defined."""


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _array_sha256(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()


def _write_json_once(path: Path, payload: dict[str, Any]) -> None:
    try:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    except FileExistsError as error:
        raise BackgroundInputError(f"refusing to replace immutable manifest: {path}") from error


def _source_groups(contract: dict[str, np.ndarray]) -> np.ndarray:
    """Use source files as independent raw-flux scale groups when available."""
    n_points = len(contract["time"])
    source = contract.get("source_file_index")
    if source is None:
        return np.zeros(n_points, dtype=np.int64)
    source = np.asarray(source)
    if source.ndim != 1 or len(source) != n_points or source.dtype.kind not in "iu":
        raise BackgroundInputError("source_file_index must be an integer array matching time")
    return source.astype(np.int64, copy=False)


def _eligible_mask(contract: dict[str, np.ndarray]) -> np.ndarray:
    eligible = np.zeros(len(contract["time"]), dtype=bool)
    eligible[contract["mask_index"]] = True
    return eligible


def _robust_scale(values: np.ndarray) -> float:
    median = float(np.median(values))
    scale = float(1.4826 * np.median(np.abs(values - median)))
    if not np.isfinite(scale) or scale <= 0.0:
        raise BackgroundInputError("eligible PDCSAP residuals have no positive robust scale")
    return scale


def generate_background(
    contract: dict[str, np.ndarray], *, mode: str, seed: int
) -> tuple[np.ndarray, dict[str, Any]]:
    """Return one raw-flux baseline and an auditable noise-model summary."""
    if mode not in MODES:
        raise BackgroundInputError(f"unsupported background mode {mode!r}")
    if int(seed) != seed or seed < 0:
        raise BackgroundInputError("seed must be a non-negative integer")
    observed = np.asarray(contract["flux"], dtype=np.float64)
    if mode == "observed_pdcsap":
        return observed.copy(), {"mode": mode, "random_noise_added": False}

    eligible = _eligible_mask(contract)
    groups = _source_groups(contract)
    if not np.any(eligible):
        raise BackgroundInputError("input has no injection-eligible cadence")
    global_center = float(np.median(observed[eligible]))
    global_mad_scale = _robust_scale(observed[eligible]) if mode == "variance_matched_iid" else None
    generated = np.empty_like(observed)
    rng = np.random.default_rng(int(seed))
    group_reports: list[dict[str, Any]] = []
    errors = None
    if mode == "quoted_error_gaussian":
        if "flux_err" not in contract:
            raise BackgroundInputError("quoted_error_gaussian requires flux_err")
        errors = np.asarray(contract["flux_err"], dtype=np.float64)
        if not np.all(np.isfinite(errors)) or np.any(errors <= 0.0):
            raise BackgroundInputError("flux_err must be finite and positive")

    for group in np.unique(groups):
        group_mask = groups == group
        fitting_mask = group_mask & eligible
        has_eligible = bool(np.any(fitting_mask))
        center = float(np.median(observed[fitting_mask])) if has_eligible else global_center
        if not np.isfinite(center):
            raise BackgroundInputError(f"source group {int(group)} has an invalid raw-flux center")
        if mode == "quoted_error_gaussian":
            sigma = errors[group_mask]
            summary_sigma = float(np.median(errors[fitting_mask])) if has_eligible else float(np.median(sigma))
        else:
            summary_sigma = _robust_scale(observed[fitting_mask]) if has_eligible else float(global_mad_scale)
            sigma = np.full(int(np.count_nonzero(group_mask)), summary_sigma, dtype=np.float64)
        generated[group_mask] = center + rng.normal(loc=0.0, scale=sigma)
        group_reports.append(
            {
                "source_file_index": int(group),
                "point_count": int(np.count_nonzero(group_mask)),
                "eligible_point_count": int(np.count_nonzero(fitting_mask)),
                "fallback_to_global_eligible_statistics": not has_eligible,
                "raw_flux_center": center,
                "noise_scale_median": summary_sigma,
            }
        )
    if not np.all(np.isfinite(generated)):
        raise RuntimeError("background generator produced non-finite flux")
    return generated, {
        "mode": mode,
        "random_noise_added": True,
        "generator": "numpy.random.Generator(PCG64)",
        "seed": int(seed),
        "centering": "Per-source-file median over injection-eligible raw PDCSAP cadences.",
        "scale": (
            "Per-cadence PDCSAP_FLUX_ERR."
            if mode == "quoted_error_gaussian"
            else "Per-source-file 1.4826 times MAD of injection-eligible raw PDCSAP flux."
        ),
        "groups": group_reports,
    }


def build_manifest(
    input_path: Path, contract: dict[str, np.ndarray], output_path: Path, details: dict[str, Any]
) -> dict[str, Any]:
    """Record the exact baseline contract without claiming a recovery result."""
    return {
        "schema_version": 1,
        "scope": "Standalone injection-background control; not an ALDERAAN result or Sagear validation product.",
        "input": {
            "filename": input_path.name,
            "sha256": _sha256_file(input_path),
            "time_sha256": _array_sha256(contract["time"]),
            "mask_index_sha256": _array_sha256(contract["mask_index"]),
            "point_count": int(len(contract["time"])),
            "injection_eligible_count": int(len(contract["mask_index"])),
        },
        "background": details,
        "output": {"filename": output_path.name, "sha256": _sha256_file(output_path)},
        "software": {"python": platform.python_version(), "numpy": np.__version__},
    }


def run_background(input_path: Path, output_dir: Path, *, mode: str, seed: int) -> Path:
    """Generate one immutable background while retaining the full input contract."""
    try:
        contract = load_input_contract(input_path)
    except InjectionInputError as error:
        raise BackgroundInputError(str(error)) from error
    if output_dir.exists() and any(output_dir.iterdir()):
        raise BackgroundInputError("output_dir must be absent or empty to preserve immutable provenance")
    output_dir.mkdir(parents=True, exist_ok=True)
    flux, details = generate_background(contract, mode=mode, seed=seed)
    output_path = output_dir / OUTPUT_NAME
    payload = dict(contract)
    payload["flux"] = flux
    np.savez_compressed(output_path, **payload)
    manifest_path = output_dir / MANIFEST_NAME
    _write_json_once(manifest_path, build_manifest(input_path, contract, output_path, details))
    return manifest_path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-npz", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mode", choices=MODES, required=True)
    parser.add_argument("--seed", type=int, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        manifest = run_background(args.input_npz, args.output_dir, mode=args.mode, seed=args.seed)
    except BackgroundInputError as error:
        _parser().error(str(error))
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
