from __future__ import annotations

"""Direct MacDougall-style eccentricity extraction from ALDERAAN samples.

The default ``e_max=0.95`` reproduces the upper limit documented in the
commented Sagear analysis source.  Use ``--e-max 0.92`` for the sampling-limit
sensitivity adopted by MacDougall, Gilbert & Petigura (2023).
"""

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import linear_sum_assignment

from common import load_config, normalize_dynesty_weights, output_dir, read_berger_table2, root_path
from extract_eccentricity_posteriors import absolute_density_error, coverage_summary, result_file_for_target, stable_seed


DAY_S = 86400.0
G_SI = 6.67430e-11
RHO_SUN_KG_M3 = 1408.0
FORMALISM = "MacDougall_Gilbert_Petigura_2023_equation_3_direct_importance"
EQUATION = (
    "rho_star_samp=(3*pi/(G*P^2))*"
    "((((1+r)^2-b^2)/sin^2((pi*T14/P)*(1+e*sin(omega))/sqrt(1-e^2))+b^2)^(3/2))"
)
EXCLUSION_COLUMNS = [
    "koi_target",
    "kepid",
    "kepoi_name",
    "koi_period",
    "results_file",
    "stage",
    "reason",
]


def macdougall_rho_star_samp(
    period_s: float | np.ndarray,
    dur14_s: np.ndarray,
    ror: np.ndarray,
    impact: np.ndarray,
    eccentricity: np.ndarray,
    omega: np.ndarray,
) -> np.ndarray:
    """Evaluate the exact MacDougall et al. (2023) Equation 3 in rho_sun."""
    period_s, dur14_s, ror, impact, eccentricity, omega = np.broadcast_arrays(
        np.asarray(period_s, dtype=float),
        np.asarray(dur14_s, dtype=float),
        np.asarray(ror, dtype=float),
        np.asarray(impact, dtype=float),
        np.asarray(eccentricity, dtype=float),
        np.asarray(omega, dtype=float),
    )
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        chord_sq = (1.0 + ror) ** 2 - impact**2
        velocity_factor = (1.0 + eccentricity * np.sin(omega)) / np.sqrt(1.0 - eccentricity**2)
        sine_sq = np.sin((np.pi * dur14_s / period_s) * velocity_factor) ** 2
        a_over_r_sq = chord_sq / sine_sq + impact**2
        rho_solar = (
            (3.0 * np.pi / (G_SI * period_s**2))
            * np.power(a_over_r_sq, 1.5)
            / RHO_SUN_KG_M3
        )
    invalid = (
        ~np.isfinite(rho_solar)
        | (period_s <= 0.0)
        | (dur14_s <= 0.0)
        | (chord_sq <= 0.0)
        | (eccentricity < 0.0)
        | (eccentricity >= 1.0)
        | (sine_sq <= 0.0)
        | (a_over_r_sq <= 0.0)
    )
    return np.where(invalid, np.nan, rho_solar)


def stellar_density_parameters(
    planet: pd.Series,
    *,
    allow_missing_error_fallback: bool = False,
) -> tuple[float, float, float]:
    rho_true = 10.0 ** float(planet["rho_log"])
    err_hi = absolute_density_error(planet.get("rho_log_upper", np.nan))
    err_lo = absolute_density_error(planet.get("rho_log_lower", np.nan))
    invalid_hi = not np.isfinite(err_hi) or err_hi <= 0.0
    invalid_lo = not np.isfinite(err_lo) or err_lo <= 0.0
    if invalid_hi or invalid_lo:
        if not allow_missing_error_fallback:
            raise ValueError("missing or invalid published stellar-density uncertainty")
        if invalid_hi:
            err_hi = 0.13 * rho_true
        if invalid_lo:
            err_lo = 0.13 * rho_true
    return rho_true, err_hi, err_lo


def density_log_likelihood(
    rho_model: np.ndarray,
    rho_true: float,
    err_hi: float,
    err_lo: float,
    mode: str,
) -> np.ndarray:
    """Gaussian density log likelihood with an explicit error convention."""
    if mode == "symmetric-average":
        sigma = np.full_like(np.asarray(rho_model, dtype=float), 0.5 * (err_hi + err_lo))
    elif mode == "split":
        sigma = np.where(np.asarray(rho_model) >= rho_true, err_hi, err_lo)
    else:
        raise ValueError(f"Unknown density error mode: {mode}")
    with np.errstate(invalid="ignore", divide="ignore"):
        # The normalization cancels for symmetric errors but not when the
        # adopted sigma changes across the two sides of the split likelihood.
        return -0.5 * ((np.asarray(rho_model, dtype=float) - rho_true) / sigma) ** 2 - np.log(sigma)


def weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    quantiles: Iterable[float],
) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    quantiles = np.asarray(list(quantiles), dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights >= 0.0)
    if not valid.any() or weights[valid].sum() <= 0.0:
        return np.full(quantiles.shape, np.nan)
    order = np.argsort(values[valid], kind="mergesort")
    values = values[valid][order]
    weights = weights[valid][order]
    cdf = np.cumsum(weights) - 0.5 * weights
    cdf /= weights.sum()
    return np.interp(quantiles, cdf, values, left=values[0], right=values[-1])


def weighted_posterior_grid(
    eccentricity: np.ndarray,
    omega: np.ndarray,
    weights: np.ndarray,
    e_grid: np.ndarray,
    omega_grid: np.ndarray,
) -> np.ndarray:
    """Deposit weighted proposals into nearest e and periodic omega cells."""
    eccentricity = np.asarray(eccentricity, dtype=float)
    omega = np.mod(np.asarray(omega, dtype=float), 2.0 * np.pi)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(eccentricity) & np.isfinite(omega) & np.isfinite(weights) & (weights > 0.0)
    posterior = np.zeros((len(e_grid), len(omega_grid)), dtype=float)
    if not valid.any():
        return posterior
    e_idx = np.abs(eccentricity[valid, None] - e_grid[None, :]).argmin(axis=1)
    omega_idx = np.floor(omega[valid] * len(omega_grid) / (2.0 * np.pi)).astype(int)
    omega_idx = np.clip(omega_idx, 0, len(omega_grid) - 1)
    np.add.at(posterior, (e_idx, omega_idx), weights[valid])
    total = posterior.sum()
    if np.isfinite(total) and total > 0.0:
        posterior /= total
    return posterior


def direct_importance_posterior(
    ror: np.ndarray,
    impact: np.ndarray,
    dur14_days: np.ndarray,
    nested_weights: np.ndarray,
    period_days: float,
    rho_true: float,
    err_hi: float,
    err_lo: float,
    *,
    n_proposals: int,
    e_max: float,
    e_grid: np.ndarray,
    omega_grid: np.ndarray,
    density_error_mode: str,
    seed: int,
) -> dict[str, object]:
    """Pair ALDERAAN rows, draw uniform (e, omega), and importance reweight."""
    ror = np.asarray(ror, dtype=float)
    impact = np.asarray(impact, dtype=float)
    dur14_days = np.asarray(dur14_days, dtype=float)
    nested_weights = np.asarray(nested_weights, dtype=float)
    if not (len(ror) == len(impact) == len(dur14_days) == len(nested_weights)):
        raise ValueError("ROR, IMPACT, DUR14, and nested weights must remain row-paired")
    nested_weights = nested_weights / nested_weights.sum()
    rng = np.random.default_rng(seed)
    sample_idx = rng.choice(len(ror), size=n_proposals, replace=True, p=nested_weights)
    eccentricity = rng.uniform(0.0, e_max, size=n_proposals)
    omega = rng.uniform(-0.5 * np.pi, 1.5 * np.pi, size=n_proposals)
    rho_model = macdougall_rho_star_samp(
        period_days * DAY_S,
        dur14_days[sample_idx] * DAY_S,
        ror[sample_idx],
        impact[sample_idx],
        eccentricity,
        omega,
    )
    loglike = density_log_likelihood(rho_model, rho_true, err_hi, err_lo, density_error_mode)
    valid = np.isfinite(loglike) & np.isfinite(rho_model) & (rho_model > 0.0)
    if not valid.any():
        raise ValueError("no finite positive MacDougall density proposals")
    loglike_valid = loglike[valid]
    importance_weights = np.exp(loglike_valid - np.max(loglike_valid))
    weight_sum = importance_weights.sum()
    if not np.isfinite(weight_sum) or weight_sum <= 0.0:
        raise ValueError("importance weights have zero or invalid sum")
    importance_weights /= weight_sum
    e_valid = eccentricity[valid]
    omega_valid = omega[valid]
    posterior = weighted_posterior_grid(e_valid, omega_valid, importance_weights, e_grid, omega_grid)
    if posterior.sum() <= 0.0:
        raise ValueError("weighted posterior grid has zero mass")
    quantiles = weighted_quantile(e_valid, importance_weights, [0.16, 0.5, 0.84])
    return {
        "posterior": posterior,
        "e_pdf": posterior.sum(axis=1),
        "e_quantiles": quantiles,
        "proposal_count": int(n_proposals),
        "valid_proposal_count": int(valid.sum()),
        "valid_proposal_fraction": float(valid.mean()),
        "nested_ess": float(1.0 / np.sum(nested_weights**2)),
        "importance_ess": float(1.0 / np.sum(importance_weights**2)),
    }


def robust_period_regression(
    ttimes: np.ndarray,
    indices: np.ndarray | None = None,
    out_flag: np.ndarray | None = None,
) -> float:
    """Fit transit time against integer epoch with iterative MAD clipping."""
    times = np.asarray(ttimes, dtype=float)
    epochs = np.arange(len(times), dtype=float) if indices is None else np.asarray(indices, dtype=float)
    keep = np.isfinite(times) & np.isfinite(epochs)
    if out_flag is not None:
        keep &= np.asarray(out_flag) == 0
    if keep.sum() < 2:
        return np.nan
    x, y = epochs[keep], times[keep]
    for _ in range(6):
        if len(x) < 2 or np.ptp(x) <= 0.0:
            return np.nan
        slope, intercept = np.polyfit(x, y, 1)
        residual = y - (intercept + slope * x)
        center = np.median(residual)
        mad = 1.4826 * np.median(np.abs(residual - center))
        if not np.isfinite(mad) or mad == 0.0:
            break
        new_keep = np.abs(residual - center) <= 5.0 * mad
        if new_keep.all() or new_keep.sum() < 2:
            break
        x, y = x[new_keep], y[new_keep]
    return float(slope) if np.isfinite(slope) and slope > 0.0 else np.nan


def read_alderaan_planets(hdul: fits.HDUList) -> list[dict[str, float]]:
    npl = int(hdul[0].header.get("NPL", 0))
    planets: list[dict[str, float]] = []
    for index in range(npl):
        name = f"TTIMES_{index:02d}"
        period = np.nan
        if name in hdul:
            data = hdul[name].data
            names = set(data.names or [])
            period = robust_period_regression(
                data["TTIME"],
                data["INDEX"] if "INDEX" in names else None,
                data["OUT_FLAG"] if "OUT_FLAG" in names else None,
            )
        planets.append({"alderaan_index": index, "period_days": period})
    return planets


def match_planets_by_period(
    planets: pd.DataFrame,
    alderaan_planets: list[dict[str, float]],
    period_tol: float,
) -> list[tuple[int, int, float, float]]:
    """Globally minimize relative period mismatch with one-to-one assignment."""
    if planets.empty or not alderaan_planets:
        return []
    koi_periods = pd.to_numeric(planets.reset_index(drop=True)["koi_period"], errors="coerce").to_numpy(float)
    fit_periods = np.asarray([p["period_days"] for p in alderaan_planets], dtype=float)
    cost = np.full((len(koi_periods), len(fit_periods)), 1e6, dtype=float)
    for i, koi_period in enumerate(koi_periods):
        if not np.isfinite(koi_period) or koi_period <= 0.0:
            continue
        valid = np.isfinite(fit_periods) & (fit_periods > 0.0)
        cost[i, valid] = np.abs(koi_period - fit_periods[valid]) / koi_period
    rows, cols = linear_sum_assignment(cost)
    matches = []
    for row, col in zip(rows, cols):
        if cost[row, col] <= period_tol:
            matches.append(
                (int(row), int(alderaan_planets[col]["alderaan_index"]), float(fit_periods[col]), float(cost[row, col]))
            )
    return sorted(matches, key=lambda item: item[0])


def _sample_frame(data: fits.FITS_rec) -> pd.DataFrame:
    array = np.array(data)
    if not array.dtype.isnative:
        array = array.byteswap().view(array.dtype.newbyteorder("="))
    return pd.DataFrame(array)


def _exclusion(planet: pd.Series, results_file: Path, stage: str, reason: str) -> dict[str, object]:
    return {
        "koi_target": planet.get("koi_target"),
        "kepid": planet.get("kepid"),
        "kepoi_name": planet.get("kepoi_name"),
        "koi_period": planet.get("koi_period"),
        "results_file": str(results_file),
        "stage": stage,
        "reason": reason,
    }


def process_target(
    results_file: Path,
    planets: pd.DataFrame,
    out_dir: Path,
    e_grid: np.ndarray,
    omega_grid: np.ndarray,
    *,
    n_proposals: int,
    e_max: float,
    density_error_mode: str,
    period_tol: float,
    min_importance_ess: float,
    allow_density_error_fallback: bool,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    with fits.open(results_file, memmap=False) as hdul:
        if "SAMPLES" not in hdul:
            return [], [_exclusion(row, results_file, "fits", "missing SAMPLES HDU") for _, row in planets.iterrows()]
        samples = _sample_frame(hdul["SAMPLES"].data)
        fit_planets = read_alderaan_planets(hdul)
    if "LN_WT" not in samples:
        return [], [_exclusion(row, results_file, "fits", "missing LN_WT column") for _, row in planets.iterrows()]
    nested_weights = normalize_dynesty_weights(samples["LN_WT"].to_numpy(float))
    matches = match_planets_by_period(planets, fit_planets, period_tol)
    matched = {match[0] for match in matches}
    summaries: list[dict[str, object]] = []
    excluded = [
        _exclusion(row, results_file, "matching", f"no one-to-one period match within relative tolerance {period_tol:g}")
        for index, row in planets.reset_index(drop=True).iterrows()
        if index not in matched
    ]
    for planet_index, alderaan_index, fit_period, relative_difference in matches:
        planet = planets.reset_index(drop=True).iloc[planet_index]
        required = [f"ROR_{alderaan_index}", f"IMPACT_{alderaan_index}", f"DUR14_{alderaan_index}"]
        missing = [column for column in required if column not in samples]
        if missing:
            excluded.append(_exclusion(planet, results_file, "samples", f"missing paired columns: {','.join(missing)}"))
            continue
        try:
            rho_true, err_hi, err_lo = stellar_density_parameters(
                planet,
                allow_missing_error_fallback=allow_density_error_fallback,
            )
            if not np.isfinite(rho_true) or rho_true <= 0.0:
                raise ValueError("invalid stellar density")
            result = direct_importance_posterior(
                samples[required[0]].to_numpy(float),
                samples[required[1]].to_numpy(float),
                samples[required[2]].to_numpy(float),
                nested_weights,
                float(planet["koi_period"]),
                rho_true,
                err_hi,
                err_lo,
                n_proposals=n_proposals,
                e_max=e_max,
                e_grid=e_grid,
                omega_grid=omega_grid,
                density_error_mode=density_error_mode,
                seed=stable_seed("macdougall-direct", planet.get("kepoi_name"), n_proposals, e_max, density_error_mode),
            )
        except (ValueError, FloatingPointError) as exc:
            excluded.append(_exclusion(planet, results_file, "importance", str(exc)))
            continue
        e16, e50, e84 = np.asarray(result["e_quantiles"], dtype=float)
        low_ess = bool(float(result["importance_ess"]) < min_importance_ess)
        qc_reasons = "importance_ess_below_threshold" if low_ess else ""
        out_file = out_dir / f"{planet['koi_target']}_{str(planet['kepoi_name']).replace('.', '_')}_eomega_posterior.npz"
        np.savez_compressed(
            out_file,
            e_grid=e_grid,
            omega_grid=omega_grid,
            posterior=result["posterior"],
            e_pdf=result["e_pdf"],
            include_transit_prior=False,
            formalism=FORMALISM,
            formalism_equation=EQUATION,
            proposal_e_prior=f"Uniform(0,{e_max:g})",
            proposal_omega_prior="Uniform(-pi/2,3pi/2)",
            density_error_mode=density_error_mode,
            no_zeta_kde=True,
            posterior_source="alderaan_direct_importance",
            impact_mode="alderaan",
            e_max=e_max,
            proposal_count=result["proposal_count"],
            valid_proposal_count=result["valid_proposal_count"],
            nested_ess=result["nested_ess"],
            importance_ess=result["importance_ess"],
            kepid=int(planet["kepid"]),
            kepoi_name=str(planet["kepoi_name"]),
            koi_target=str(planet["koi_target"]),
            alderaan_planet_index=alderaan_index,
            alderaan_period_days=fit_period,
            period_relative_difference=relative_difference,
        )
        summaries.append(
            {
                "koi_target": planet["koi_target"],
                "kepid": int(planet["kepid"]),
                "kepoi_name": planet["kepoi_name"],
                "koi_period": float(planet["koi_period"]),
                "alderaan_planet_index": alderaan_index,
                "alderaan_period_days": fit_period,
                "period_relative_difference": relative_difference,
                "disk": planet.get("disk"),
                "system": planet.get("system"),
                "e16": e16,
                "e50": e50,
                "e84": e84,
                "posterior_file": str(out_file),
                "posterior_source": "alderaan_direct_importance",
                "impact_mode": "alderaan",
                "formalism": FORMALISM,
                "include_transit_prior": False,
                "density_error_mode": density_error_mode,
                "e_max": e_max,
                "nested_sample_count": len(samples),
                "nested_ess": result["nested_ess"],
                "proposal_count": result["proposal_count"],
                "valid_proposal_count": result["valid_proposal_count"],
                "valid_proposal_fraction": result["valid_proposal_fraction"],
                "importance_ess": result["importance_ess"],
                "minimum_importance_ess": min_importance_ess,
                "qc_importance_ess_low": low_ess,
                "qc_primary_exclude": low_ess,
                "qc_reasons": qc_reasons,
                "qc_manifest_version": "direct_importance_v1",
                "zeta_median": np.nan,
                "zeta_p16": np.nan,
                "zeta_p84": np.nan,
                "rho_true_solar": rho_true,
                "rho_err_hi_solar": err_hi,
                "rho_err_lo_solar": err_lo,
            }
        )
    return summaries, excluded


def _add_density(sample: pd.DataFrame, berger: pd.DataFrame) -> pd.DataFrame:
    density_columns = ["rho_log", "rho_log_upper", "rho_log_lower"]
    lookup = berger[["kepid", *density_columns]].drop_duplicates("kepid")
    merged = sample.merge(lookup, on="kepid", how="left", suffixes=("", "_berger"))
    for column in density_columns:
        alternate = f"{column}_berger"
        if alternate in merged:
            if column in sample:
                merged[column] = merged[column].combine_first(merged[alternate])
            else:
                merged[column] = merged[alternate]
            merged = merged.drop(columns=alternate)
    return merged


def main() -> None:
    parser = argparse.ArgumentParser(description="Exact MacDougall-style direct importance extractor for ALDERAAN samples.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--target", default=None)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--results-dir", default=None)
    parser.add_argument("--posterior-subdir", default="eccentricity_posteriors_direct")
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--coverage-out", default=None)
    parser.add_argument("--excluded-out", default=None, help="Explicit CSV manifest; written even when no planets are excluded.")
    parser.add_argument("--n-proposals", type=int, default=150_000)
    parser.add_argument(
        "--e-max",
        type=float,
        default=0.95,
        help="Uniform-e proposal upper bound: 0.95 in Sagear commented source (default); use 0.92 for MacDougall sensitivity.",
    )
    parser.add_argument("--density-error-mode", choices=["symmetric-average", "split"], default="symmetric-average")
    parser.add_argument(
        "--allow-density-error-fallback",
        action="store_true",
        help="Diagnostic only: replace missing density errors with 13%% instead of excluding the planet.",
    )
    parser.add_argument("--period-tol", type=float, default=0.01)
    parser.add_argument(
        "--min-importance-ess",
        type=float,
        default=100.0,
        help="Flag (do not silently delete) posterior rows with smaller direct-importance ESS.",
    )
    parser.add_argument("--max-targets", type=int, default=None)
    args = parser.parse_args()
    if not 0.0 < args.e_max < 1.0:
        parser.error("--e-max must be strictly between 0 and 1")
    if args.n_proposals <= 0:
        parser.error("--n-proposals must be positive")
    if args.min_importance_ess <= 0:
        parser.error("--min-importance-ess must be positive")

    cfg = load_config(args.config)
    sample_path = Path(args.sample) if args.sample else output_dir() / "canonical_sample_replication.csv"
    if not sample_path.exists():
        raise FileNotFoundError(
            "No canonical replication sample exists yet. Pass --sample explicitly, or complete the "
            "classifier reconstruction and write outputs/canonical_sample_replication.csv. The old "
            "canonical_sample_diagnostic.csv is intentionally not consumed implicitly."
        )
    sample = _add_density(pd.read_csv(sample_path), read_berger_table2(cfg))
    out_dir = output_dir() / args.posterior_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    flat_results_dir = Path(args.results_dir) if args.results_dir else None
    project = root_path(cfg, "alderaan_project")
    if flat_results_dir is None and project is None:
        raise ValueError("Either --results-dir or alderaan_project must be configured")
    run_id = args.run_id or cfg["alderaan"]["run_id"]
    e_grid = np.linspace(0.0, args.e_max, int(cfg["alderaan"]["eccentricity_grid_size"]))
    omega_grid = np.linspace(0.0, 2.0 * np.pi, int(cfg["alderaan"]["omega_grid_size"]), endpoint=False)
    targets = sorted(sample["koi_target"].dropna().astype(str).unique())
    if args.target:
        targets = [args.target]
    if args.max_targets is not None:
        targets = targets[: args.max_targets]

    summaries: list[dict[str, object]] = []
    exclusions: list[dict[str, object]] = []
    for target in targets:
        target_planets = sample[sample["koi_target"].astype(str) == target].sort_values("koi_period").reset_index(drop=True)
        results_file = result_file_for_target(target, flat_results_dir, project, run_id)
        if not results_file.exists():
            exclusions.extend(
                _exclusion(row, results_file, "discovery", "ALDERAAN results file not found")
                for _, row in target_planets.iterrows()
            )
            continue
        target_summary, target_excluded = process_target(
            results_file,
            target_planets,
            out_dir,
            e_grid,
            omega_grid,
            n_proposals=args.n_proposals,
            e_max=args.e_max,
            density_error_mode=args.density_error_mode,
            period_tol=args.period_tol,
            min_importance_ess=args.min_importance_ess,
            allow_density_error_fallback=args.allow_density_error_fallback,
        )
        summaries.extend(target_summary)
        exclusions.extend(target_excluded)

    summary_path = Path(args.summary_out) if args.summary_out else output_dir() / "eccentricity_posterior_summary_direct.csv"
    coverage_path = Path(args.coverage_out) if args.coverage_out else output_dir() / "eccentricity_posterior_coverage_direct.csv"
    excluded_path = Path(args.excluded_out) if args.excluded_out else output_dir() / "eccentricity_posterior_excluded_direct.csv"
    for path in (summary_path, coverage_path, excluded_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(summaries)
    summary.to_csv(summary_path, index=False)
    coverage_summary(sample, summary if not summary.empty else pd.DataFrame({"kepoi_name": []})).to_csv(coverage_path, index=False)
    pd.DataFrame(exclusions, columns=EXCLUSION_COLUMNS).to_csv(excluded_path, index=False)
    print(f"Wrote {len(summary)} direct posterior summaries: {summary_path}")
    print(f"Wrote coverage: {coverage_path}")
    print(f"Wrote {len(exclusions)} explicit exclusions: {excluded_path}")


if __name__ == "__main__":
    main()
