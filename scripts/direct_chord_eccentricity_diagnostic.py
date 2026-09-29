"""Run a bounded direct ``(e, omega)`` check on one ALDERAAN fixture.

This diagnostic preserves ALDERAAN's public detrended-light-curve chord
likelihood.  It replaces the sampled duration with the exact finite-duration
mapping from a declared stellar density, eccentricity, argument of periastron,
radius ratio, and impact parameter.  It is intentionally limited to a
single-planet fixture.  It is not an exact eccentric-orbit light-curve model
and it must never be used to update the canonical population archive.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.special import erfinv


RHO_SUN_KG_M3 = 1408.0
G_SI = 6.67430e-11
ECCENTRICITY_MAX = 0.95


class DirectChordError(ValueError):
    """Raised when a direct-chord diagnostic is outside its declared scope."""


def sha256(path: Path) -> str:
    """Return a content hash for a declared diagnostic input."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def a_over_rstar(period_days: float, rho_star_solar: float) -> float:
    """Compute ``a/Rstar`` from period and a fixed mean stellar density."""
    if not np.isfinite(period_days) or period_days <= 0:
        raise DirectChordError("period_days must be finite and positive")
    if not np.isfinite(rho_star_solar) or rho_star_solar <= 0:
        raise DirectChordError("rho_star_solar must be finite and positive")
    period_seconds = period_days * 86400.0
    return float((G_SI * rho_star_solar * RHO_SUN_KG_M3 * period_seconds**2 / (3.0 * np.pi)) ** (1.0 / 3.0))


def finite_duration_days(
    period_days: float,
    radius_ratio: float,
    impact_parameter: float,
    eccentricity: float,
    omega_radians: float,
    rho_star_solar: float,
) -> float:
    """Map physical parameters to ALDERAAN's finite chord duration ``T14``.

    The duration follows the Sagear/MacDougall finite-duration expression.
    It is deliberately evaluated separately from the native ALDERAAN chord
    renderer so the diagnostic changes only the parameterization, not that
    renderer.
    """
    values = (period_days, radius_ratio, impact_parameter, eccentricity, omega_radians, rho_star_solar)
    if not all(np.isfinite(value) for value in values):
        raise DirectChordError("duration inputs must be finite")
    if radius_ratio <= 0 or radius_ratio >= 1:
        raise DirectChordError("radius_ratio must lie in (0, 1)")
    if eccentricity < 0 or eccentricity >= 1:
        raise DirectChordError("eccentricity must lie in [0, 1)")
    if impact_parameter < 0 or impact_parameter > 1.0 + radius_ratio:
        raise DirectChordError("impact_parameter falls outside the transit chord support")
    aor = a_over_rstar(period_days, rho_star_solar)
    numerator = (1.0 + radius_ratio) ** 2 - impact_parameter**2
    denominator = aor**2 - impact_parameter**2
    if numerator <= 0 or denominator <= 0:
        raise DirectChordError("duration geometry is not physical")
    arcsin_argument = math.sqrt(numerator / denominator)
    if not 0 < arcsin_argument < 1:
        raise DirectChordError("duration geometry has no valid finite chord")
    g = (1.0 + eccentricity * math.sin(omega_radians)) / math.sqrt(1.0 - eccentricity**2)
    if not np.isfinite(g) or g <= 0:
        raise DirectChordError("photoeccentric velocity factor is invalid")
    return float(period_days / math.pi * math.asin(arcsin_argument) / g)


def direct_prior_transform(unit_cube: np.ndarray) -> np.ndarray:
    """Transform one unit-cube row into the fixed eight-parameter basis."""
    u = np.asarray(unit_cube, dtype=np.float64)
    if u.shape != (8,) or np.any(~np.isfinite(u)) or np.any((u < 0) | (u > 1)):
        raise DirectChordError("unit_cube must contain eight values on [0, 1]")
    # These are the same Normal(0, 0.1) ephemeris priors ALDERAAN uses.
    c0 = 0.1 * np.sqrt(2.0) * erfinv((2.0 * u[0] - 1.0) / (1.0 + 1.0e-12))
    c1 = 0.1 * np.sqrt(2.0) * erfinv((2.0 * u[1] - 1.0) / (1.0 + 1.0e-12))
    radius_ratio = float(np.exp(np.log(1.0e-5) + u[2] * np.log(0.99 / 1.0e-5)))
    impact_parameter = (1.0 + radius_ratio) * u[3]
    eccentricity = ECCENTRICITY_MAX * u[4]
    omega = -0.5 * np.pi + 2.0 * np.pi * u[5]
    q1, q2 = u[6], u[7]
    return np.array([c0, c1, radius_ratio, impact_parameter, eccentricity, omega, q1, q2])


def native_alderaan_vector(direct_parameters: np.ndarray, period_days: float, rho_star_solar: float) -> np.ndarray:
    """Map direct parameters into ALDERAAN's native ``C0,C1,r,b,T14,q1,q2`` basis."""
    direct = np.asarray(direct_parameters, dtype=np.float64)
    if direct.shape != (8,):
        raise DirectChordError("direct parameter vector must have length eight")
    c0, c1, radius_ratio, impact, eccentricity, omega, q1, q2 = direct
    duration = finite_duration_days(period_days, radius_ratio, impact, eccentricity, omega, rho_star_solar)
    return np.array([c0, c1, radius_ratio, impact, duration, q1, q2], dtype=np.float64)


def fixture_preflight(project_dir: Path, run_id: str, target: str, rho_star_solar: float) -> dict[str, Any]:
    """Record and validate files needed to reconstruct a one-planet fixture."""
    results_dir = project_dir / "Results" / run_id / target
    required = [
        results_dir / f"{target}_transit_parameters.csv",
        results_dir / f"{target}_lc_filtered.fits",
        results_dir / f"{target}_00_quick.ttvs",
    ]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise DirectChordError("missing required fixture inputs: " + "; ".join(missing))
    with required[0].open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1 or rows[0].get("npl") != "1":
        raise DirectChordError("the direct diagnostic requires exactly one catalog planet")
    if rows[0].get("koi_id") != target:
        raise DirectChordError("catalog target does not match the requested target")
    return {
        "schema_version": 1,
        "scope": "One-planet direct chord diagnostic; not a canonical population result.",
        "project_dir": str(project_dir.resolve()),
        "run_id": run_id,
        "target": target,
        "rho_star_solar": float(rho_star_solar),
        "catalog_planets": 1,
        "required_inputs": [
            {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in required
        ],
        "model_boundary": (
            "Uses ALDERAAN dynesty_helpers.lnlike after reconstructing its public "
            "detrended chord inputs. It does not render an eccentric orbit directly."
        ),
    }


def _build_native_context(project_dir: Path, run_id: str, target: str) -> tuple[Any, tuple[Any, ...], float, tuple[float, float]]:
    """Rebuild the one-planet data structures consumed by ``dynhelp.lnlike``.

    This follows the public ALDERAAN transit-shape script for long cadence
    fixture data.  It deliberately excludes short cadence and multi-planet
    systems from this diagnostic rather than silently approximating either.
    """
    import batman
    import alderaan.dynesty_helpers as dynhelp
    import alderaan.io as io
    from alderaan.astro import make_transit_mask, set_oversample_factor
    from alderaan.constants import lcit
    from alderaan.Ephemeris import Ephemeris

    results_dir = project_dir / "Results" / run_id / target
    catalog = io.parse_catalog(results_dir / f"{target}_transit_parameters.csv", "Kepler", target)
    n_planets = int(catalog.npl.to_numpy()[0])
    if n_planets != 1:
        raise DirectChordError("this diagnostic is intentionally limited to one-planet fixtures")
    period = float(catalog.period.to_numpy()[0])
    depth = float(catalog.depth.to_numpy()[0]) * 1.0e-6
    catalog_duration = float(catalog.duration.to_numpy()[0]) / 24.0
    u1 = float(catalog.limbdark_1.to_numpy()[0])
    u2 = float(catalog.limbdark_2.to_numpy()[0])
    lc = io.load_detrended_lightcurve(results_dir / f"{target}_lc_filtered.fits")
    data = np.genfromtxt(results_dir / f"{target}_00_quick.ttvs")
    if data.ndim != 2 or data.shape[1] < 3:
        raise DirectChordError("quick-TTV table must have epoch, fitted, and regular times")
    transit_indices = np.asarray(data[:, 0], dtype=int)
    regular_times = np.asarray(data[:, 2], dtype=float)
    if len(transit_indices) < 2:
        raise DirectChordError("at least two timing rows are required")
    ephemeris = Ephemeris(transit_indices, regular_times)
    mask = make_transit_mask(lc.time, regular_times, masksize=max(1.0 / 24.0, 1.5 * catalog_duration))
    count_expected = max(1, int(np.floor(catalog_duration / lcit)))
    selected_quarters = []
    times: list[np.ndarray] = []
    fluxes: list[np.ndarray] = []
    errors: list[np.ndarray] = []
    warped_times: list[list[np.ndarray]] = [[]]
    warped_indices: list[list[np.ndarray]] = [[]]
    oversample = np.zeros(18, dtype=int)
    exposure = np.zeros(18, dtype=float)
    offsets: list[np.ndarray | None] = [None] * 18
    oversample_lc = int(set_oversample_factor(np.array([period]), np.array([depth]), np.array([catalog_duration]), lc.flux, lc.error, lcit))
    for quarter in range(18):
        in_quarter = lc.quarter == quarter
        use = mask & in_quarter
        if int(np.sum(use)) <= count_expected:
            continue
        selected_quarters.append(quarter)
        t = np.asarray(lc.time[use], dtype=float)
        times.append(t)
        fluxes.append(np.asarray(lc.flux[use], dtype=float))
        errors.append(np.asarray(lc.error[use], dtype=float))
        warped_time, warped_index = ephemeris._warp_times(t, return_inds=True)
        warped_times[0].append(warped_time)
        warped_indices[0].append((warped_index - transit_indices[-1] // 2) / (transit_indices[-1] / 2))
        oversample[quarter] = oversample_lc
        exposure[quarter] = lcit
        offsets[quarter] = np.linspace(-lcit / 2.0, lcit / 2.0, oversample_lc)
    if not selected_quarters:
        raise DirectChordError("the fixture has no retained long-cadence transit windows")
    theta = batman.TransitParams()
    theta.per = ephemeris.period
    theta.t0 = 0.0
    theta.rp = math.sqrt(depth)
    theta.b = 0.5
    theta.T14 = catalog_duration
    theta.u = [u1, u2]
    theta.limb_dark = "quadratic"
    phot_args = {
        "time": times,
        "flux": fluxes,
        "error": errors,
        "quarters": np.asarray(selected_quarters, dtype=int),
        "warped_t": warped_times,
        "warped_x": warped_indices,
        "exptime": exposure,
        "oversample": oversample,
        "texp_offsets": offsets,
    }
    centered_indices = transit_indices - transit_indices[-1] // 2
    ephem_args = {
        "transit_inds": [centered_indices],
        "transit_times": [regular_times],
        "transit_legx": [centered_indices / (transit_indices[-1] / 2)],
    }
    # Match the native log-uniform T14 support. Dividing the interim prior in
    # post-model importance sampling does not restore likelihood outside this
    # finite support, so a matched direct comparison must retain the boundary.
    duration_support = (float(lcit), float(3.0 * catalog_duration))
    return dynhelp.lnlike, (1, [theta], ephem_args, phot_args, [u1, u2], None), period, duration_support


def run_direct_nested(
    project_dir: Path,
    run_id: str,
    target: str,
    rho_star_solar: float,
    *,
    seed: int,
    maxcall: int,
) -> dict[str, np.ndarray]:
    """Run a bounded direct nested sampler against ALDERAAN's public chord likelihood."""
    import dynesty

    native_loglike, native_args, period, duration_support = _build_native_context(project_dir, run_id, target)

    def loglike(direct: np.ndarray) -> float:
        try:
            native = native_alderaan_vector(direct, period, rho_star_solar)
        except DirectChordError:
            return -1.0e300
        if not duration_support[0] <= native[4] <= duration_support[1]:
            return -1.0e300
        return float(native_loglike(native, *native_args))

    sampler = dynesty.DynamicNestedSampler(loglike, direct_prior_transform, 8, rstate=np.random.default_rng(seed), bound="multi", sample="rwalk")
    sampler.run_nested(maxcall=maxcall, print_progress=False)
    return {
        "samples": np.asarray(sampler.results.samples),
        "logwt": np.asarray(sampler.results.logwt),
        "logz": np.asarray(sampler.results.logz),
        "logl": np.asarray(sampler.results.logl),
        "period_days": np.asarray([period]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--rho-star-solar", type=float, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--execute", action="store_true", help="Run bounded direct nested sampling after preflight.")
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--maxcall", type=int, default=20000)
    parser.add_argument("--output", type=Path, help="Output NPZ required with --execute.")
    args = parser.parse_args()
    manifest = fixture_preflight(args.project_dir, args.run_id, args.target, args.rho_star_solar)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if args.execute:
        if args.output is None:
            raise DirectChordError("--output is required with --execute")
        results = run_direct_nested(args.project_dir, args.run_id, args.target, args.rho_star_solar, seed=args.seed, maxcall=args.maxcall)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.output, **results)
        manifest["execution"] = {
            "seed": args.seed,
            "maxcall": args.maxcall,
            "output": str(args.output.resolve()),
            "output_sha256": sha256(args.output),
            "nested_rows": int(len(results["samples"])),
        }
        args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
