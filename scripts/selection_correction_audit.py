from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import output_dir
from hierarchical_rayleigh import POPULATIONS, fit_from_mass_matrix


SELECTION_MODES = [
    "none",
    "legacy_forward_norm",
    "manuscript_reciprocal",
    "manuscript_reciprocal_with_norm",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compare transit-selection correction conventions. The Sagear manuscript's "
            "commented HBM equation uses reciprocal transit probability, while the older "
            "pipeline used forward transit probability plus population normalization."
        )
    )
    parser.add_argument(
        "--summary",
        default=None,
        help="Merged posterior summary. Defaults to paired-exact QC-primary output.",
    )
    parser.add_argument("--tag", default="QC_PRIMARY")
    parser.add_argument("--n-sigma", type=int, default=2000)
    args = parser.parse_args()

    out = output_dir()
    summary_path = Path(args.summary) if args.summary else out / "eccentricity_posterior_summary_merged_paired_exact_qcprimary.csv"
    summary = pd.read_csv(summary_path)
    sigmas = np.linspace(1e-4, 1.0, args.n_sigma)

    rows = []
    for mode in SELECTION_MODES:
        for disk, system, population in POPULATIONS:
            sub = summary[(summary["disk"] == disk) & (summary["system"] == system)].reset_index(drop=True)
            if len(sub) < 5:
                rows.append({"selection_mode": mode, "population": population, "n": len(sub), "status": "skipped_n_lt_5"})
                continue
            masses, e_grid = load_mass_matrix(sub["posterior_file"].tolist(), mode)
            if mode == "legacy_forward_norm":
                fit = fit_from_mass_matrix(masses, e_grid, sigmas, apply_transit_selection=True)
            elif mode == "manuscript_reciprocal_with_norm":
                fit = fit_from_mass_matrix_with_custom_norm(
                    masses,
                    e_grid,
                    sigmas,
                    normalizer_mode="reciprocal_omega_average",
                )
            else:
                fit = fit_from_mass_matrix(masses, e_grid, sigmas, apply_transit_selection=False)
            rows.append({"selection_mode": mode, "population": population, "n": len(sub), "status": "ok", **fit})

    result = pd.DataFrame(rows)
    result_path = out / f"selection_correction_rayleigh_comparison_{args.tag}.csv"
    pivot_path = out / f"selection_correction_rayleigh_pivot_{args.tag}.csv"
    result.to_csv(result_path, index=False)
    result.pivot(index="selection_mode", columns="population", values="expected_e").reset_index().to_csv(pivot_path, index=False)
    print(f"Wrote {result_path}")
    print(f"Wrote {pivot_path}")
    print(result.pivot(index="selection_mode", columns="population", values="expected_e").reset_index().to_string(index=False))


def load_mass_matrix(files: list[str], selection_mode: str) -> tuple[np.ndarray, np.ndarray]:
    masses = []
    e_grid = None
    for file in files:
        data = np.load(file)
        include_transit_prior = bool(np.asarray(data["include_transit_prior"]).item()) if "include_transit_prior" in data.files else None
        if include_transit_prior is not False:
            raise ValueError(
                f"Selection audit requires include_transit_prior=False posteriors; got {include_transit_prior}: {file}"
            )
        this_e = np.asarray(data["e_grid"], dtype=float)
        this_omega = np.asarray(data["omega_grid"], dtype=float)
        posterior = np.asarray(data["posterior"], dtype=float)
        posterior = np.nan_to_num(np.clip(posterior, 0.0, np.inf), nan=0.0, posinf=0.0, neginf=0.0)
        total = posterior.sum()
        posterior = np.full_like(posterior, 1.0 / posterior.size) if total <= 0 else posterior / total
        weight = transit_probability_weight(this_e, this_omega)
        if selection_mode == "none":
            e_mass = posterior.sum(axis=1)
        elif selection_mode == "legacy_forward_norm":
            e_mass = (posterior * weight).sum(axis=1)
        elif selection_mode in {"manuscript_reciprocal", "manuscript_reciprocal_with_norm"}:
            e_mass = (posterior / weight).sum(axis=1)
        else:
            raise ValueError(f"Unknown selection mode: {selection_mode}")
        if e_grid is None:
            e_grid = this_e
        elif not np.array_equal(e_grid, this_e):
            e_mass = np.interp(e_grid, this_e, e_mass, left=0.0, right=0.0)
        masses.append(e_mass)
    if e_grid is None:
        raise ValueError("No posterior files supplied")
    return np.vstack(masses), e_grid


def transit_probability_weight(e_grid: np.ndarray, omega_grid: np.ndarray) -> np.ndarray:
    e, omega = np.meshgrid(e_grid, omega_grid, indexing="ij")
    weight = (1.0 + e * np.sin(omega)) / (1.0 - e**2)
    return np.clip(weight, 1e-12, np.inf)


def fit_from_mass_matrix_with_custom_norm(
    mass_matrix: np.ndarray,
    e_grid: np.ndarray,
    sigmas: np.ndarray,
    normalizer_mode: str,
) -> dict:
    rays = []
    normalizers = []
    for sigma in sigmas:
        ray = (e_grid / sigma**2) * np.exp(-(e_grid**2) / (2.0 * sigma**2))
        ray = ray / trapezoid(ray, e_grid)
        rays.append(ray)
        if normalizer_mode == "reciprocal_omega_average":
            # Mean over omega of 1 / p_transit(e, omega) is sqrt(1-e^2)
            # for the dimensionless factor used in the manuscript comments.
            normalizers.append(trapezoid(ray * np.sqrt(np.clip(1.0 - e_grid**2, 0.0, np.inf)), e_grid))
        else:
            raise ValueError(f"Unknown normalizer mode: {normalizer_mode}")
    rays = np.vstack(rays).T
    normalizers = np.asarray(normalizers, dtype=float)
    terms = (mass_matrix @ rays) / normalizers
    log_terms = np.log(np.clip(terms, 1e-300, None))
    lls = np.sum(log_terms, axis=0)
    return fit_from_ll(sigmas, lls)


def fit_from_ll(sigmas: np.ndarray, lls: np.ndarray) -> dict:
    best = int(np.argmax(lls))
    weights = posterior_weights_from_ll(sigmas, lls)
    sigma_q16, sigma_q50, sigma_q84 = weighted_quantile(sigmas, weights, [0.16, 0.5, 0.84])
    expected_grid = sigmas * np.sqrt(np.pi / 2.0)
    e_q16, e_q50, e_q84 = weighted_quantile(expected_grid, weights, [0.16, 0.5, 0.84])
    boundary = bool(best == 0 or best == len(sigmas) - 1 or weights[0] > 0.01 or weights[-1] > 0.01)
    return {
        "sigma_rayleigh": sigma_q50,
        "sigma_rayleigh_lo": sigma_q16,
        "sigma_rayleigh_hi": sigma_q84,
        "sigma_map": sigmas[best],
        "expected_e": e_q50,
        "expected_e_lo": e_q16,
        "expected_e_hi": e_q84,
        "expected_e_map": expected_grid[best],
        "ll_max": lls[best],
        "transit_selection_applied": True,
        "sigma_at_grid_lower_edge": bool(best == 0),
        "sigma_at_grid_upper_edge": bool(best == len(sigmas) - 1),
        "boundary_flag": boundary,
        "posterior_mass_at_lower_edge": float(weights[0]),
        "posterior_mass_at_upper_edge": float(weights[-1]),
    }


def trapezoid(y: np.ndarray, x: np.ndarray | None = None, axis: int = -1) -> np.ndarray:
    fn = getattr(np, "trapezoid", None)
    if fn is None:
        fn = np.trapz
    return fn(y, x=x, axis=axis)


def posterior_weights_from_ll(sigmas: np.ndarray, lls: np.ndarray) -> np.ndarray:
    widths = grid_cell_widths(sigmas)
    raw = np.exp(lls - np.nanmax(lls)) * widths
    total = raw.sum()
    if not np.isfinite(total) or total <= 0:
        return np.full(len(sigmas), 1.0 / len(sigmas))
    return raw / total


def grid_cell_widths(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    edges = np.empty(len(x) + 1, dtype=float)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = max(0.0, x[0] - (edges[1] - x[0]))
    edges[-1] = x[-1] + (x[-1] - edges[-2])
    return np.diff(edges)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: list[float]) -> tuple[float, ...]:
    order = np.argsort(values)
    v = np.asarray(values)[order]
    w = np.asarray(weights)[order]
    cdf = np.cumsum(w)
    cdf = cdf / cdf[-1]
    return tuple(float(np.interp(q, cdf, v)) for q in quantiles)


if __name__ == "__main__":
    main()
