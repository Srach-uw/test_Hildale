from __future__ import annotations

"""Diagnostic circular-density overlap audit for one ALDERAAN planet.

This script is deliberately bounded to a per-planet comparison. It is not a
significance test, a validation of eccentricity extraction, or a population
inference result. It compares a Dynesty-weighted circular transit density with
a declared stellar-density distribution on log10 density.
"""

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from common import normalize_dynesty_weights
from extract_eccentricity_posteriors_direct import (
    DAY_S,
    _sample_frame,
    macdougall_rho_star_samp,
    paired_period_samples,
)


REQUIRED_SAMPLE_COLUMNS = ("P", "T14", "ROR", "IMPACT", "LN_WT")
DIAGNOSTIC_LABEL = "diagnostic_only_not_a_significance_or_population_claim"
DEFAULT_GRID_SIZE = 1024
LOG_QUANTILES = np.array([0.15865525393145707, 0.5, 0.8413447460685429])


class CircularDensityOverlapError(ValueError):
    """Raised when a bounded density-overlap audit cannot be evaluated."""


def sha256_file(path: str | Path) -> str:
    """Return the SHA-256 digest of one immutable diagnostic input file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_finite_positive(values: np.ndarray, label: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise CircularDensityOverlapError(f"{label} must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise CircularDensityOverlapError(f"{label} must contain only finite positive values")
    return values


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, probabilities: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cdf = np.cumsum(sorted_weights)
    cdf /= cdf[-1]
    # Invert the weighted empirical CDF. Interpolating across an almost-zero
    # Dynesty-weight row would manufacture density mass where none exists.
    indices = np.searchsorted(cdf, probabilities, side="left")
    return sorted_values[np.clip(indices, 0, len(sorted_values) - 1)]


def stellar_samples_from_quantiles(
    quantiles: Mapping[str, float],
    *,
    grid_size: int = DEFAULT_GRID_SIZE,
) -> np.ndarray:
    """Build a deterministic log-density distribution from declared 16/50/84 quantiles.

    The approximation linearly interpolates the quantile function in log10
    density between the supplied 16th, 50th, and 84th percentiles and extends
    the two outer segments linearly to the bounded evaluation grid. It is a
    transparent diagnostic approximation, not a replacement for posterior
    samples when those samples are available.
    """
    required = ("p16", "p50", "p84")
    missing = [key for key in required if key not in quantiles]
    if missing:
        raise CircularDensityOverlapError("stellar quantiles missing: " + ", ".join(missing))
    values = np.array([quantiles[key] for key in required], dtype=float)
    _require_finite_positive(values, "stellar quantiles")
    if not values[0] <= values[1] <= values[2]:
        raise CircularDensityOverlapError("stellar quantiles must satisfy p16 <= p50 <= p84")
    if not isinstance(grid_size, int) or grid_size < 32:
        raise CircularDensityOverlapError("grid_size must be an integer of at least 32")

    probabilities = (np.arange(grid_size, dtype=float) + 0.5) / grid_size
    log_values = np.log10(values)
    # Bounded tails use the adjacent measured slope rather than inventing a family.
    slopes = np.diff(log_values) / np.diff(LOG_QUANTILES)
    log_samples = np.interp(probabilities, LOG_QUANTILES, log_values)
    low = probabilities < LOG_QUANTILES[0]
    high = probabilities > LOG_QUANTILES[-1]
    log_samples[low] = log_values[0] + slopes[0] * (probabilities[low] - LOG_QUANTILES[0])
    log_samples[high] = log_values[-1] + slopes[-1] * (probabilities[high] - LOG_QUANTILES[-1])
    return np.power(10.0, log_samples)


def _validated_transit_samples(samples: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    missing = [column for column in REQUIRED_SAMPLE_COLUMNS if column not in samples.columns]
    if missing:
        raise CircularDensityOverlapError("paired ALDERAAN samples missing columns: " + ", ".join(missing))
    if samples.empty:
        raise CircularDensityOverlapError("paired ALDERAAN samples must contain at least one row")
    arrays = tuple(samples[column].to_numpy(dtype=float) for column in REQUIRED_SAMPLE_COLUMNS)
    period, duration, ror, impact, log_weight = arrays
    if not all(len(array) == len(period) for array in arrays):
        raise CircularDensityOverlapError("paired ALDERAAN columns must have equal length")
    if not np.all(np.isfinite(period)) or np.any(period <= 0.0):
        raise CircularDensityOverlapError("P must contain finite positive days")
    if not np.all(np.isfinite(duration)) or np.any(duration <= 0.0):
        raise CircularDensityOverlapError("T14 must contain finite positive days")
    if not np.all(np.isfinite(ror)) or np.any(ror <= 0.0):
        raise CircularDensityOverlapError("ROR must contain finite positive values")
    if not np.all(np.isfinite(impact)) or np.any(impact < 0.0):
        raise CircularDensityOverlapError("IMPACT must contain finite nonnegative values")
    if not np.any(np.isfinite(log_weight)):
        raise CircularDensityOverlapError("LN_WT must contain at least one finite Dynesty log weight")
    return period, duration, ror, impact, log_weight


def read_paired_samples_from_results(
    results_fits_path: str | Path,
    planet_index: int,
) -> pd.DataFrame:
    """Read one row-paired ALDERAAN planet posterior from a result FITS file.

    The returned P, T14, ROR, IMPACT, and LN_WT columns are in the units used
    by ``summarize_overlap``. This performs no fitting or resampling.
    """
    try:
        from astropy.io import fits
    except ModuleNotFoundError as error:
        raise CircularDensityOverlapError(
            "astropy is required to read an ALDERAAN result FITS file"
        ) from error
    results_fits_path = Path(results_fits_path)
    if not results_fits_path.is_file():
        raise FileNotFoundError(f"ALDERAAN result FITS file does not exist: {results_fits_path}")
    if not isinstance(planet_index, int) or planet_index < 0:
        raise CircularDensityOverlapError("planet_index must be a nonnegative integer")
    try:
        with fits.open(results_fits_path, memmap=False) as hdul:
            if "SAMPLES" not in hdul:
                raise CircularDensityOverlapError("ALDERAAN result FITS is missing the SAMPLES HDU")
            npl = int(hdul[0].header.get("NPL", 0))
            if planet_index >= npl:
                raise CircularDensityOverlapError(
                    f"planet_index {planet_index} is outside the result FITS NPL={npl}"
                )
            samples = _sample_frame(hdul["SAMPLES"].data)
            required = (
                f"DUR14_{planet_index}",
                f"ROR_{planet_index}",
                f"IMPACT_{planet_index}",
                "LN_WT",
            )
            missing = [column for column in required if column not in samples]
            if missing:
                raise CircularDensityOverlapError(
                    "ALDERAAN result FITS is missing paired columns: "
                    + ", ".join(missing)
                )
            period = paired_period_samples(samples, hdul, planet_index)
    except OSError as error:
        raise CircularDensityOverlapError(f"could not read ALDERAAN result FITS: {error}") from error
    return pd.DataFrame(
        {
            "P": period,
            "T14": samples[f"DUR14_{planet_index}"].to_numpy(float),
            "ROR": samples[f"ROR_{planet_index}"].to_numpy(float),
            "IMPACT": samples[f"IMPACT_{planet_index}"].to_numpy(float),
            "LN_WT": samples["LN_WT"].to_numpy(float),
        }
    )


def summarize_overlap(
    samples: pd.DataFrame,
    *,
    stellar_density_samples: Sequence[float] | np.ndarray | None = None,
    stellar_density_quantiles: Mapping[str, float] | None = None,
    tolerance_dex: float = 0.1,
    grid_size: int = DEFAULT_GRID_SIZE,
) -> dict[str, Any]:
    """Summarize a bounded, independent circular-versus-stellar density comparison.

    ``P`` and ``T14`` are interpreted as days. The paired transit rows are
    weighted by ``LN_WT``. The stellar distribution is either supplied as
    positive density samples in solar units or deterministically reconstructed
    from explicitly declared p16/p50/p84 quantiles.
    """
    if not np.isfinite(tolerance_dex) or tolerance_dex < 0.0:
        raise CircularDensityOverlapError("tolerance_dex must be finite and nonnegative")
    if not isinstance(grid_size, int) or grid_size < 32 or grid_size > 4096:
        raise CircularDensityOverlapError("grid_size must be an integer from 32 through 4096")
    if (stellar_density_samples is None) == (stellar_density_quantiles is None):
        raise CircularDensityOverlapError(
            "provide exactly one of stellar_density_samples or stellar_density_quantiles"
        )

    period, duration, ror, impact, log_weight = _validated_transit_samples(samples)
    weights = normalize_dynesty_weights(log_weight)
    rho_circular = macdougall_rho_star_samp(
        period * DAY_S,
        duration * DAY_S,
        ror,
        impact,
        np.zeros_like(period),
        np.zeros_like(period),
    )
    valid = np.isfinite(rho_circular) & (rho_circular > 0.0) & np.isfinite(weights) & (weights > 0.0)
    if not valid.any():
        raise CircularDensityOverlapError("no finite positive circular-density samples remain")
    rho_circular = rho_circular[valid]
    weights = weights[valid]
    weights /= weights.sum()

    if stellar_density_samples is not None:
        stellar = _require_finite_positive(np.asarray(stellar_density_samples, dtype=float), "stellar density samples")
        stellar_source = "declared_samples"
    else:
        assert stellar_density_quantiles is not None
        stellar = stellar_samples_from_quantiles(stellar_density_quantiles, grid_size=grid_size)
        stellar_source = "declared_p16_p50_p84_quantiles"

    probabilities = (np.arange(grid_size, dtype=float) + 0.5) / grid_size
    circular_grid = _weighted_quantile(rho_circular, weights, probabilities)
    stellar_grid = np.quantile(stellar, probabilities, method="linear")
    # The Cartesian grid is the deterministic independent product distribution.
    delta_grid = np.log10(circular_grid[:, None]) - np.log10(stellar_grid[None, :])
    delta = delta_grid.ravel()
    delta_p16, delta_p50, delta_p84 = np.quantile(delta, LOG_QUANTILES, method="linear")
    circ_p16, circ_p50, circ_p84 = _weighted_quantile(rho_circular, weights, LOG_QUANTILES)
    star_p16, star_p50, star_p84 = np.quantile(stellar, LOG_QUANTILES, method="linear")
    ess = float(1.0 / np.sum(weights**2))
    return {
        "diagnostic_label": DIAGNOSTIC_LABEL,
        "units": {"P": "days", "T14": "days", "density": "rho_sun", "delta": "dex"},
        "configuration": {"tolerance_dex": float(tolerance_dex), "grid_size": grid_size},
        "input": {
            "alderaan_rows": int(len(samples)),
            "finite_circular_rows": int(valid.sum()),
            "stellar_density_source": stellar_source,
            "stellar_density_input_rows": int(len(stellar)),
        },
        "dynesty": {"weight_column": "LN_WT", "effective_sample_size": ess},
        "rho_circular_solar": {"p16": float(circ_p16), "p50": float(circ_p50), "p84": float(circ_p84)},
        "rho_stellar_solar": {"p16": float(star_p16), "p50": float(star_p50), "p84": float(star_p84)},
        "delta_log10_rho_circular_minus_stellar_dex": {
            "p16": float(delta_p16),
            "p50": float(delta_p50),
            "p84": float(delta_p84),
            "probability_abs_delta_within_tolerance": float(np.mean(np.abs(delta) <= tolerance_dex)),
            "zero_in_68pct_credible_interval": bool(delta_p16 <= 0.0 <= delta_p84),
        },
    }


def write_audit(
    samples_path: str | Path | None,
    output_prefix: str | Path,
    *,
    results_fits_path: str | Path | None = None,
    planet_index: int | None = None,
    stellar_density_samples_path: str | Path | None = None,
    stellar_density_column: str = "rho_star_solar",
    stellar_density_quantiles: Mapping[str, float] | None = None,
    tolerance_dex: float = 0.1,
    grid_size: int = DEFAULT_GRID_SIZE,
) -> tuple[Path, Path, dict[str, Any]]:
    """Read declared inputs and write one-row CSV plus a provenance JSON manifest."""
    if (samples_path is None) == (results_fits_path is None):
        raise CircularDensityOverlapError(
            "provide exactly one of samples_path or results_fits_path"
        )
    samples_hash: str
    results_metadata: dict[str, object] | None = None
    if samples_path is not None:
        samples_path = Path(samples_path)
        if not samples_path.exists():
            raise FileNotFoundError(f"ALDERAAN samples file does not exist: {samples_path}")
        samples = pd.read_csv(samples_path)
        samples_hash = sha256_file(samples_path)
    else:
        if planet_index is None:
            raise CircularDensityOverlapError("planet_index is required with results_fits_path")
        results_fits_path = Path(results_fits_path)
        samples = read_paired_samples_from_results(results_fits_path, planet_index)
        samples_hash = sha256_file(results_fits_path)
        results_metadata = {
            "results_fits_path": str(results_fits_path.resolve()),
            "results_fits_sha256": samples_hash,
            "planet_index": int(planet_index),
            "extraction": "row_paired_period_duration_ror_impact_dynesty_weight",
        }
    stellar_hash: str | None = None
    stellar_samples: np.ndarray | None = None
    if stellar_density_samples_path is not None:
        stellar_density_samples_path = Path(stellar_density_samples_path)
        if not stellar_density_samples_path.exists():
            raise FileNotFoundError(f"stellar density file does not exist: {stellar_density_samples_path}")
        stellar_frame = pd.read_csv(stellar_density_samples_path)
        if stellar_density_column not in stellar_frame.columns:
            raise CircularDensityOverlapError(
                f"stellar density file missing column: {stellar_density_column}"
            )
        stellar_samples = stellar_frame[stellar_density_column].to_numpy(dtype=float)
        stellar_hash = sha256_file(stellar_density_samples_path)
    summary = summarize_overlap(
        samples,
        stellar_density_samples=stellar_samples,
        stellar_density_quantiles=stellar_density_quantiles,
        tolerance_dex=tolerance_dex,
        grid_size=grid_size,
    )
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = prefix.with_name(prefix.name + "_overlap.csv")
    json_path = prefix.with_name(prefix.name + "_manifest.json")
    flat = {
        "diagnostic_label": summary["diagnostic_label"],
        "circular_rho_p16_solar": summary["rho_circular_solar"]["p16"],
        "circular_rho_p50_solar": summary["rho_circular_solar"]["p50"],
        "circular_rho_p84_solar": summary["rho_circular_solar"]["p84"],
        "stellar_rho_p16_solar": summary["rho_stellar_solar"]["p16"],
        "stellar_rho_p50_solar": summary["rho_stellar_solar"]["p50"],
        "stellar_rho_p84_solar": summary["rho_stellar_solar"]["p84"],
        "delta_log10_rho_p16_dex": summary["delta_log10_rho_circular_minus_stellar_dex"]["p16"],
        "delta_log10_rho_p50_dex": summary["delta_log10_rho_circular_minus_stellar_dex"]["p50"],
        "delta_log10_rho_p84_dex": summary["delta_log10_rho_circular_minus_stellar_dex"]["p84"],
        "probability_abs_delta_within_tolerance": summary["delta_log10_rho_circular_minus_stellar_dex"]["probability_abs_delta_within_tolerance"],
        "zero_in_68pct_credible_interval": summary["delta_log10_rho_circular_minus_stellar_dex"]["zero_in_68pct_credible_interval"],
        "dynesty_effective_sample_size": summary["dynesty"]["effective_sample_size"],
    }
    pd.DataFrame([flat]).to_csv(csv_path, index=False)
    manifest = {
        **summary,
        "provenance": {
            "samples_path": str(samples_path.resolve()) if samples_path else None,
            "samples_sha256": samples_hash if samples_path else None,
            "result_fits": results_metadata,
            "stellar_samples_path": str(stellar_density_samples_path.resolve()) if stellar_density_samples_path else None,
            "stellar_samples_sha256": stellar_hash,
            "stellar_quantiles": dict(stellar_density_quantiles) if stellar_density_quantiles else None,
            "script": Path(__file__).name,
        },
        "outputs": {"csv": str(csv_path.resolve()), "manifest": str(json_path.resolve())},
    }
    json_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return csv_path, json_path, manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Diagnostic per-planet circular-density overlap audit.")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--samples", help="CSV with paired P, T14, ROR, IMPACT, and LN_WT columns; P and T14 are days.")
    source.add_argument("--results-fits", help="ALDERAAN result FITS to read with a required --planet-index.")
    parser.add_argument("--planet-index", type=int, help="Zero-based ALDERAAN planet index, required with --results-fits.")
    parser.add_argument("--output-prefix", required=True, help="Path prefix for diagnostic CSV and JSON outputs.")
    stellar = parser.add_mutually_exclusive_group(required=True)
    stellar.add_argument("--stellar-samples", help="CSV containing declared stellar-density samples in rho_sun.")
    stellar.add_argument("--stellar-p16", type=float, help="Declared stellar-density 16th percentile in rho_sun.")
    parser.add_argument("--stellar-column", default="rho_star_solar")
    parser.add_argument("--stellar-p50", type=float)
    parser.add_argument("--stellar-p84", type=float)
    parser.add_argument("--tolerance-dex", type=float, default=0.1)
    parser.add_argument("--grid-size", type=int, default=DEFAULT_GRID_SIZE)
    args = parser.parse_args()
    quantiles = None
    if args.stellar_p16 is not None:
        if args.stellar_p50 is None or args.stellar_p84 is None:
            parser.error("--stellar-p16 requires --stellar-p50 and --stellar-p84")
        quantiles = {"p16": args.stellar_p16, "p50": args.stellar_p50, "p84": args.stellar_p84}
    elif args.stellar_p50 is not None or args.stellar_p84 is not None:
        parser.error("--stellar-p50 and --stellar-p84 require --stellar-p16")
    csv_path, json_path, _ = write_audit(
        args.samples,
        args.output_prefix,
        results_fits_path=args.results_fits,
        planet_index=args.planet_index,
        stellar_density_samples_path=args.stellar_samples,
        stellar_density_column=args.stellar_column,
        stellar_density_quantiles=quantiles,
        tolerance_dex=args.tolerance_dex,
        grid_size=args.grid_size,
    )
    print(f"DIAGNOSTIC_CSV={csv_path}")
    print(f"DIAGNOSTIC_MANIFEST={json_path}")


if __name__ == "__main__":
    main()
