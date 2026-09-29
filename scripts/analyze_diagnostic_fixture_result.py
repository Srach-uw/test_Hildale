"""Summarize one provenance-preserving circular-injection ALDERAAN result."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits

from common import normalize_dynesty_weights
from diagnostic_transit_injection import circular_first_to_fourth_duration_days
from extract_eccentricity_posteriors_direct import macdougall_rho_star_samp, weighted_quantile


class DiagnosticResultError(ValueError):
    """Raised when a diagnostic result cannot support a density comparison."""


def sha256(path: Path) -> str:
    """Return the SHA-256 digest for one immutable input file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_specification(path: Path) -> dict[str, float]:
    """Load and validate the circular geometry declared for one fixture."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        # The standalone injection manifest uses ``specification`` at the root,
        # while the prepared ALDERAAN fixture records it beneath ``injection``.
        specification = payload.get("specification")
        if specification is None:
            specification = payload["injection"]["specification"]
        period_days = float(specification["period_days"])
        rho_star_solar = float(specification["rho_star_solar"])
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise DiagnosticResultError("injection manifest lacks a valid circular specification") from exc
    if not np.isfinite(period_days) or period_days <= 0:
        raise DiagnosticResultError("injection period must be finite and positive")
    if not np.isfinite(rho_star_solar) or rho_star_solar <= 0:
        raise DiagnosticResultError("injected stellar density must be finite and positive")
    return {
        "period_days": period_days,
        "rho_star_solar": rho_star_solar,
        "impact_parameter": float(specification.get("impact_parameter", np.nan)),
        "radius_ratio": float(specification.get("radius_ratio", np.nan)),
        "seed": int(specification.get("seed", -1)),
    }


def positive_weighted_quantiles(values: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    """Summarize finite positive draws while retaining their paired weights."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    keep = np.isfinite(values) & (values > 0.0) & np.isfinite(weights) & (weights >= 0.0)
    if not keep.any() or weights[keep].sum() <= 0.0:
        raise DiagnosticResultError("no finite positive draws for diagnostic summary")
    kept_weights = weights[keep]
    kept_weights = kept_weights / kept_weights.sum()
    quantiles = weighted_quantile(values[keep], kept_weights, [0.025, 0.16, 0.5, 0.84, 0.975])
    return dict(zip(("p2p5", "p16", "p50", "p84", "p97p5"), map(float, quantiles)))


def analyze_result(
    result_fits: Path,
    injection_manifest: Path,
    *,
    planet_index: int = 0,
    warning_threshold_dex: float = 0.10,
) -> dict[str, Any]:
    """Compute weighted circular-density recovery metrics from one FITS result."""
    if planet_index < 0:
        raise DiagnosticResultError("planet index must be non-negative")
    if not np.isfinite(warning_threshold_dex) or warning_threshold_dex <= 0:
        raise DiagnosticResultError("warning threshold must be finite and positive")
    specification = read_specification(injection_manifest)
    suffix = f"_{planet_index}"
    required = ("LN_WT", f"DUR14{suffix}", f"ROR{suffix}", f"IMPACT{suffix}")

    with fits.open(result_fits, memmap=False) as hdul:
        if "SAMPLES" not in hdul:
            raise DiagnosticResultError("result FITS lacks SAMPLES extension")
        samples = hdul["SAMPLES"].data
        names = set(samples.names or ())
        missing = sorted(set(required) - names)
        if missing:
            raise DiagnosticResultError(f"result FITS lacks required columns: {', '.join(missing)}")
        weights = normalize_dynesty_weights(np.asarray(samples["LN_WT"], dtype=float))
        duration_days = np.asarray(samples[f"DUR14{suffix}"], dtype=float)
        radius_ratio = np.asarray(samples[f"ROR{suffix}"], dtype=float)
        impact = np.asarray(samples[f"IMPACT{suffix}"], dtype=float)

    valid = (
        np.isfinite(weights)
        & np.isfinite(duration_days)
        & np.isfinite(radius_ratio)
        & np.isfinite(impact)
        & (weights >= 0)
        & (duration_days > 0)
        & (radius_ratio > 0)
    )
    if not valid.any() or weights[valid].sum() <= 0:
        raise DiagnosticResultError("result FITS has no usable nested samples")
    weights = weights[valid]
    weights = weights / weights.sum()
    duration_days = duration_days[valid]
    radius_ratio = radius_ratio[valid]
    impact = impact[valid]
    rho_circular = macdougall_rho_star_samp(
        specification["period_days"] * 86400.0,
        duration_days * 86400.0,
        radius_ratio,
        impact,
        np.zeros_like(duration_days),
        np.zeros_like(duration_days),
    )
    rho_valid = np.isfinite(rho_circular) & (rho_circular > 0)
    if not rho_valid.any() or weights[rho_valid].sum() <= 0:
        raise DiagnosticResultError("no finite positive circular-density draws")
    rho_weights = weights[rho_valid]
    rho_weights = rho_weights / rho_weights.sum()
    rho_circular = rho_circular[rho_valid]
    rho_summary = positive_weighted_quantiles(rho_circular, rho_weights)
    rho025 = rho_summary["p2p5"]
    rho16 = rho_summary["p16"]
    rho50 = rho_summary["p50"]
    rho84 = rho_summary["p84"]
    rho975 = rho_summary["p97p5"]
    ratio = float(rho50 / specification["rho_star_solar"])
    delta_dex = float(np.log10(ratio))
    interval_tolerance = 16.0 * np.finfo(float).eps * max(
        abs(float(rho16)), abs(float(rho84)), specification["rho_star_solar"], 1.0
    )

    def transit_shape_quantiles(values: np.ndarray) -> dict[str, float]:
        quantiles = weighted_quantile(values, weights, [0.16, 0.5, 0.84])
        return dict(zip(("p16", "p50", "p84"), map(float, quantiles)))

    injected_duration_days = circular_first_to_fourth_duration_days(
        specification["period_days"],
        specification["rho_star_solar"],
        specification["radius_ratio"],
        specification["impact_parameter"],
    )
    counterfactual_density = {
        "paired_posterior": rho_summary,
        "injected_impact_posterior_duration_radius_ratio": positive_weighted_quantiles(
            macdougall_rho_star_samp(
                specification["period_days"] * 86400.0,
                duration_days * 86400.0,
                radius_ratio,
                np.full_like(impact, specification["impact_parameter"]),
                np.zeros_like(duration_days),
                np.zeros_like(duration_days),
            ),
            weights,
        ),
        "posterior_impact_injected_duration_radius_ratio": positive_weighted_quantiles(
            macdougall_rho_star_samp(
                specification["period_days"] * 86400.0,
                np.full_like(duration_days, injected_duration_days * 86400.0),
                radius_ratio,
                impact,
                np.zeros_like(duration_days),
                np.zeros_like(duration_days),
            ),
            weights,
        ),
        "posterior_impact_duration_injected_radius_ratio": positive_weighted_quantiles(
            macdougall_rho_star_samp(
                specification["period_days"] * 86400.0,
                duration_days * 86400.0,
                np.full_like(radius_ratio, specification["radius_ratio"]),
                impact,
                np.zeros_like(duration_days),
                np.zeros_like(duration_days),
            ),
            weights,
        ),
    }

    return {
        "scope": "Diagnostic circular-injection result. It is not a population result or a Sagear replication claim.",
        "result_fits": str(result_fits.resolve()),
        "result_fits_sha256": sha256(result_fits),
        "injection_manifest": str(injection_manifest.resolve()),
        "injection_manifest_sha256": sha256(injection_manifest),
        "planet_index": int(planet_index),
        "injected": specification,
        "nested_rows": int(len(valid)),
        "valid_transit_shape_rows": int(len(duration_days)),
        "valid_density_rows": int(len(rho_circular)),
        "density_effective_sample_size": float(1.0 / np.sum(rho_weights**2)),
        "recovered_rho_circular_solar": {
            "p2p5": float(rho025),
            "p16": float(rho16),
            "p50": float(rho50),
            "p84": float(rho84),
            "p97p5": float(rho975),
        },
        "recovered_to_injected_median_ratio": ratio,
        "log10_recovered_to_injected_median_ratio": delta_dex,
        "injected_within_recovered_68pct_interval": bool(
            rho16 - interval_tolerance
            <= specification["rho_star_solar"]
            <= rho84 + interval_tolerance
        ),
        "injected_within_recovered_95pct_interval": bool(
            rho025 - interval_tolerance
            <= specification["rho_star_solar"]
            <= rho975 + interval_tolerance
        ),
        "warning_threshold_dex": float(warning_threshold_dex),
        "median_bias_exceeds_warning_threshold": bool(abs(delta_dex) > warning_threshold_dex),
        "recovered_transit_shape": {
            "duration_days": transit_shape_quantiles(duration_days),
            "radius_ratio": transit_shape_quantiles(radius_ratio),
            "impact_parameter": transit_shape_quantiles(impact),
        },
        "bias_decomposition": {
            "injected_duration_days": float(injected_duration_days),
            "counterfactual_circular_density_solar": counterfactual_density,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-fits", required=True, type=Path)
    parser.add_argument("--injection-manifest", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--planet-index", type=int, default=0)
    parser.add_argument("--warning-threshold-dex", type=float, default=0.10)
    args = parser.parse_args()
    if args.output_json.exists():
        raise DiagnosticResultError(f"refusing to replace existing output: {args.output_json}")
    output = analyze_result(
        args.result_fits,
        args.injection_manifest,
        planet_index=args.planet_index,
        warning_threshold_dex=args.warning_threshold_dex,
    )
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
