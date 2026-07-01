from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from common import (
    circular_duration_days,
    load_config,
    normalize_dynesty_weights,
    output_dir,
    read_berger_table2,
    root_path,
    stellar_a_over_r,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract e,omega posterior grids from ALDERAAN results.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--target", default=None, help="Optional single KOI target such as K00752.")
    parser.add_argument("--run-id", default=None, help="Override config alderaan.run_id for project-style results.")
    parser.add_argument("--summary-out", default=None, help="Output CSV path for posterior summary.")
    parser.add_argument("--coverage-out", default=None, help="Output CSV path for coverage summary.")
    parser.add_argument("--posterior-subdir", default="eccentricity_posteriors", help="Subdirectory under outputs for e,omega posterior grids.")
    parser.add_argument(
        "--results-dir",
        default=None,
        help=(
            "Directory containing flat ALDERAAN *-results.fits files. If omitted, use "
            "alderaan_project/Results/<run_id>/<target>/<target>-results.fits."
        ),
    )
    parser.add_argument("--period-tol", type=float, default=0.01, help="Maximum relative period mismatch.")
    parser.add_argument("--max-targets", type=int, default=None, help="Diagnostic limit on number of systems to process.")
    parser.add_argument(
        "--impact-mode",
        choices=["geometric", "alderaan", "catalog"],
        default="geometric",
        help=(
            "How to handle impact parameter in the circular-duration calculation. "
            "Default geometric marginalizes b ~ U(0, 1+Rp/R*) because ALDERAAN b can be prior-dominated."
        ),
    )
    parser.add_argument(
        "--include-transit-prior",
        action="store_true",
        default=False,
        help=(
            "Legacy/diagnostic mode: multiply individual e,omega posteriors by the geometric transit probability. "
            "Leave off for Sagear-style hierarchical inference, which applies the reciprocal transit probability "
            "inside the population likelihood."
        ),
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = output_dir() / args.posterior_subdir
    out.mkdir(parents=True, exist_ok=True)

    flat_results_dir = Path(args.results_dir) if args.results_dir else None
    project = root_path(cfg, "alderaan_project")
    if flat_results_dir is None and project is None:
        raise ValueError("Either --results-dir or alderaan_project must be configured")
    run_id = args.run_id or cfg["alderaan"]["run_id"]
    sample_path = Path(args.sample) if args.sample else output_dir() / "canonical_sample_diagnostic.csv"
    if not sample_path.exists():
        raise FileNotFoundError(
            "Missing sample table for mapping ALDERAAN results back to planets. "
            f"Run diagnose_sample.py --allow-fallback-gmm first, or pass --sample explicitly: {sample_path}"
        )
    sample = pd.read_csv(sample_path)
    berger = read_berger_table2(cfg)
    sample = sample.merge(berger[["kepid", "rho_log", "rho_log_upper", "rho_log_lower"]], on="kepid", how="left", suffixes=("", "_berger2"))

    e_grid = np.linspace(0.001, 0.95, int(cfg["alderaan"]["eccentricity_grid_size"]))
    omega_grid = np.linspace(0.0, 2.0 * np.pi, int(cfg["alderaan"]["omega_grid_size"]))

    summaries = []
    targets = sorted(sample["koi_target"].dropna().unique())
    if args.target:
        targets = [args.target]
    if args.max_targets is not None:
        targets = targets[: args.max_targets]

    checked_targets = 0
    missing_results = 0

    for target in targets:
        results_file = result_file_for_target(target, flat_results_dir, project, run_id)
        if not results_file.exists():
            missing_results += 1
            continue
        checked_targets += 1
        target_rows = sample[sample["koi_target"] == target].sort_values("koi_period").reset_index(drop=True)
        target_summaries = process_target(
            results_file,
            target_rows,
            e_grid,
            omega_grid,
            out,
            cfg,
            args.include_transit_prior,
            args.period_tol,
            args.impact_mode,
        )
        summaries.extend(target_summaries)
        if checked_targets % 100 == 0:
            print(f"Processed {checked_targets} ALDERAAN systems; matched {len(summaries)} planets")

    if summaries:
        summary = pd.DataFrame(summaries)
        summary_path = Path(args.summary_out) if args.summary_out else output_dir() / "eccentricity_posterior_summary.csv"
        coverage_path = Path(args.coverage_out) if args.coverage_out else output_dir() / "eccentricity_posterior_coverage.csv"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        coverage_path.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(summary_path, index=False)
        coverage = coverage_summary(sample, summary)
        coverage.to_csv(coverage_path, index=False)
        print(f"Wrote {len(summaries)} posterior summaries")
        print(f"Wrote: {summary_path}")
        print(f"Wrote: {coverage_path}")
        print("\nCoverage by disk/system:")
        print(coverage.to_string(index=False))
        print(f"\nSystems checked: {checked_targets}; systems with no results file: {missing_results}")
    else:
        print("No ALDERAAN results found. Run the validation ALDERAAN batch first.")


def result_file_for_target(target: str, flat_results_dir: Path | None, project: Path | None, run_id: str) -> Path:
    if flat_results_dir is not None:
        return flat_results_dir / f"{target}-results.fits"
    assert project is not None
    return project / "Results" / run_id / target / f"{target}-results.fits"


def process_target(
    results_file: Path,
    planets: pd.DataFrame,
    e_grid: np.ndarray,
    omega_grid: np.ndarray,
    out_dir: Path,
    cfg: dict,
    include_transit_prior: bool,
    period_tol: float,
    impact_mode: str,
) -> list[dict]:
    with fits.open(results_file, memmap=False) as hdul:
        samples = pd.DataFrame(np.array(hdul["SAMPLES"].data).byteswap().newbyteorder())
        alderaan_planets = read_alderaan_planets(hdul)

    logwt = samples["LN_WT"].to_numpy(float)
    weights = normalize_dynesty_weights(logwt)
    n_draw = min(int(cfg["alderaan"]["posterior_resample_size"]), len(samples))
    rng = np.random.default_rng(42)
    idx = rng.choice(np.arange(len(samples)), size=n_draw, replace=True, p=weights)
    s = samples.iloc[idx].reset_index(drop=True)

    summaries = []
    matches = match_planets_by_period(planets, alderaan_planets, period_tol)
    for planet_idx, alderaan_idx, period_fit, rel_diff in matches:
        planet = planets.iloc[planet_idx]
        suffix = f"_{alderaan_idx}"
        ror = s[f"ROR{suffix}"].to_numpy(float)
        impact_alderaan = s[f"IMPACT{suffix}"].to_numpy(float)
        dur14 = s[f"DUR14{suffix}"].to_numpy(float)
        # ALDERAAN stores duration in days in the model parameterization.
        impact = impact_draws(planet, ror, impact_alderaan, impact_mode, n_draw)
        rho_draw = draw_rho_log(planet, n_draw)
        a_rs = stellar_a_over_r(rho_draw, np.full(n_draw, planet["koi_period"], dtype=float))
        tcirc = circular_duration_days(
            np.full(n_draw, planet["koi_period"], dtype=float), impact, ror, a_rs
        )
        zeta = dur14 / tcirc
        zeta = zeta[np.isfinite(zeta) & (zeta > 0)]
        if len(zeta) < 30:
            continue

        post = posterior_grid_from_zeta(zeta, e_grid, omega_grid, include_transit_prior)
        e_pdf = post.sum(axis=1)
        e_pdf = e_pdf / np.trapz(e_pdf, e_grid)
        e_cdf = np.cumsum(e_pdf)
        e_cdf = e_cdf / e_cdf[-1]
        e50 = np.interp(0.5, e_cdf, e_grid)
        e16 = np.interp(0.16, e_cdf, e_grid)
        e84 = np.interp(0.84, e_cdf, e_grid)

        out_file = out_dir / f"{planet['koi_target']}_{planet['kepoi_name'].replace('.', '_')}_eomega_posterior.npz"
        np.savez_compressed(
            out_file,
            e_grid=e_grid,
            omega_grid=omega_grid,
            posterior=post,
            e_pdf=e_pdf,
            include_transit_prior=bool(include_transit_prior),
            kepid=int(planet["kepid"]),
            kepoi_name=str(planet["kepoi_name"]),
            koi_target=str(planet["koi_target"]),
            alderaan_planet_index=int(alderaan_idx),
            alderaan_period_days=float(period_fit),
            period_relative_difference=float(rel_diff),
            zeta_median=float(np.nanmedian(zeta)),
            zeta_p16=float(np.nanpercentile(zeta, 16)),
            zeta_p84=float(np.nanpercentile(zeta, 84)),
            impact_mode=str(impact_mode),
        )
        summaries.append(
            {
                "koi_target": planet["koi_target"],
                "kepid": int(planet["kepid"]),
                "kepoi_name": planet["kepoi_name"],
                "koi_period": float(planet["koi_period"]),
                "alderaan_planet_index": int(alderaan_idx),
                "alderaan_period_days": float(period_fit),
                "period_relative_difference": float(rel_diff),
                "disk": planet.get("disk"),
                "system": planet.get("system"),
                "e16": e16,
                "e50": e50,
                "e84": e84,
                "zeta_median": float(np.nanmedian(zeta)),
                "zeta_p16": float(np.nanpercentile(zeta, 16)),
                "zeta_p84": float(np.nanpercentile(zeta, 84)),
                "n_zeta": int(len(zeta)),
                "impact_mode": str(impact_mode),
                "posterior_file": str(out_file),
            }
        )
    return summaries


def impact_draws(
    planet: pd.Series,
    ror: np.ndarray,
    impact_alderaan: np.ndarray,
    impact_mode: str,
    n_draw: int,
) -> np.ndarray:
    if impact_mode == "alderaan":
        return impact_alderaan
    if impact_mode == "catalog":
        b = float(planet.get("koi_impact", np.nan))
        if np.isfinite(b) and b >= 0:
            return np.full(n_draw, b, dtype=float)
    r_med = float(np.nanmedian(ror))
    if not np.isfinite(r_med) or r_med <= 0:
        r_med = float(planet.get("koi_ror", 0.02))
    rng = np.random.default_rng(abs(hash(("impact", str(planet.get("kepoi_name"))))) % (2**32))
    return rng.uniform(0.0, max(1e-3, 1.0 + r_med - 1e-3), size=n_draw)


def read_alderaan_planets(hdul: fits.HDUList) -> list[dict[str, float]]:
    npl = int(hdul[0].header.get("NPL", 0))
    planets = []
    for n in range(npl):
        hdu_name = f"TTIMES_{n:02d}"
        period = np.nan
        if hdu_name in hdul:
            period = period_from_ttimes(hdul[hdu_name].data["TTIME"])
        planets.append({"alderaan_index": n, "period_days": period})
    return planets


def period_from_ttimes(ttimes: np.ndarray) -> float:
    t = np.sort(np.asarray(ttimes, dtype=float))
    t = t[np.isfinite(t)]
    if len(t) < 2:
        return np.nan
    dt = np.diff(t)
    dt = dt[np.isfinite(dt) & (dt > 0.05)]
    if len(dt) == 0:
        return np.nan
    p0 = np.nanmin(dt)
    multiples = np.maximum(1, np.round(dt / p0))
    return float(np.nanmedian(dt / multiples))


def match_planets_by_period(
    planets: pd.DataFrame,
    alderaan_planets: list[dict[str, float]],
    period_tol: float,
) -> list[tuple[int, int, float, float]]:
    candidates = []
    for pi, planet in planets.reset_index(drop=True).iterrows():
        p_koi = float(planet["koi_period"])
        for ap in alderaan_planets:
            p_fit = float(ap["period_days"])
            if not np.isfinite(p_fit) or p_fit <= 0:
                continue
            rel = abs(p_koi - p_fit) / p_koi
            if rel <= period_tol:
                candidates.append((rel, int(pi), int(ap["alderaan_index"]), p_fit))
    candidates.sort(key=lambda x: x[0])
    used_planets: set[int] = set()
    used_alderaan: set[int] = set()
    matches = []
    for rel, pi, ai, p_fit in candidates:
        if pi in used_planets or ai in used_alderaan:
            continue
        used_planets.add(pi)
        used_alderaan.add(ai)
        matches.append((pi, ai, p_fit, rel))
    matches.sort(key=lambda x: x[0])
    return matches


def draw_rho_log(planet: pd.Series, n_draw: int) -> np.ndarray:
    mean_log = float(planet["rho_log"])
    rho = 10.0**mean_log
    err_hi = absolute_density_error(planet.get("rho_log_upper", np.nan))
    err_lo = absolute_density_error(planet.get("rho_log_lower", np.nan))
    if not np.isfinite(err_hi) or err_hi <= 0:
        err_hi = 0.13 * rho
    if not np.isfinite(err_lo) or err_lo <= 0:
        err_lo = 0.13 * rho
    rng = np.random.default_rng(abs(hash(str(planet.get("kepoi_name")))) % (2**32))
    draws = np.empty(n_draw, dtype=float)
    filled = 0
    while filled < n_draw:
        batch = max(1024, n_draw - filled)
        z = np.abs(rng.standard_normal(batch))
        high = rng.random(batch) < err_hi / (err_hi + err_lo)
        vals = np.where(high, rho + z * err_hi, rho - z * err_lo)
        vals = vals[np.isfinite(vals) & (vals > 1e-8)]
        take = min(len(vals), n_draw - filled)
        if take > 0:
            draws[filled : filled + take] = vals[:take]
            filled += take
    return np.log10(draws)


def absolute_density_error(value: object) -> float:
    value = float(value)
    if not np.isfinite(value):
        return np.nan
    # Berger+2020 table2 stores rho, E_rho, and e_rho as log10 values in
    # solar-density units. E_rho/e_rho are log10 absolute errors, not bounds.
    return 10.0**value


def posterior_grid_from_zeta(
    zeta_samples: np.ndarray,
    e_grid: np.ndarray,
    omega_grid: np.ndarray,
    include_transit_prior: bool,
) -> np.ndarray:
    e, omega = np.meshgrid(e_grid, omega_grid, indexing="ij")
    zeta_model = np.sqrt(1.0 - e**2) / (1.0 + e * np.sin(omega))
    likelihood = density_from_samples(zeta_samples, zeta_model)
    posterior = likelihood
    if include_transit_prior:
        posterior = posterior * (1.0 + e * np.sin(omega)) / (1.0 - e**2)
    posterior = np.clip(posterior, 0.0, np.inf)
    norm = posterior.sum()
    if norm <= 0:
        posterior[:] = 1.0 / posterior.size
    else:
        posterior /= norm
    return posterior


def density_from_samples(samples: np.ndarray, query: np.ndarray) -> np.ndarray:
    samples = np.asarray(samples, dtype=float)
    samples = samples[np.isfinite(samples) & (samples > 0)]
    if len(samples) < 30:
        return np.ones_like(query, dtype=float)
    lo = max(0.05, float(np.nanpercentile(samples, 0.2)) * 0.7)
    hi = min(6.0, float(np.nanpercentile(samples, 99.8)) * 1.3)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = 0.05, 6.0
    bins = np.geomspace(lo, hi, 360)
    hist, edges = np.histogram(samples, bins=bins, density=True)
    kernel_x = np.arange(-5, 6)
    kernel = np.exp(-0.5 * (kernel_x / 1.5) ** 2)
    kernel /= kernel.sum()
    hist = np.convolve(hist, kernel, mode="same")
    centers = np.sqrt(edges[:-1] * edges[1:])
    return np.interp(query.ravel(), centers, hist, left=0.0, right=0.0).reshape(query.shape)


def coverage_summary(sample: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    matched = set(summary["kepoi_name"].astype(str))
    rows = []
    for (disk, system), sub in sample.groupby(["disk", "system"]):
        n_total = len(sub)
        n_matched = int(sub["kepoi_name"].astype(str).isin(matched).sum())
        rows.append(
            {
                "disk": disk,
                "system": system,
                "sample_planets": n_total,
                "posterior_planets": n_matched,
                "missing_planets": n_total - n_matched,
                "coverage_fraction": n_matched / n_total if n_total else np.nan,
            }
        )
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
