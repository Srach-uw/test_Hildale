from __future__ import annotations

import numpy as np

from hierarchical_rayleigh import (
    eccentricity_mass_from_posterior,
    fit_from_mass_matrix,
)


TRUE_SIGMA = 0.30
E_MAX = 0.95


def selected_rayleigh_point_posteriors(
    *,
    count: int,
    seed: int,
    e_grid: np.ndarray,
    omega_grid: np.ndarray,
) -> np.ndarray:
    """Return exact-point interim posteriors after a declared transit selection."""
    rng = np.random.default_rng(seed)
    accepted_e: list[float] = []
    accepted_omega: list[float] = []
    probability_ceiling = (1.0 + E_MAX) / (1.0 - E_MAX**2)

    while len(accepted_e) < count:
        eccentricity = rng.rayleigh(TRUE_SIGMA, size=4096)
        omega = rng.uniform(0.0, 2.0 * np.pi, size=len(eccentricity))
        selection = (1.0 + eccentricity * np.sin(omega)) / (1.0 - eccentricity**2)
        keep = (eccentricity < E_MAX) & (
            rng.uniform(0.0, 1.0, size=len(eccentricity)) < selection / probability_ceiling
        )
        accepted_e.extend(eccentricity[keep])
        accepted_omega.extend(omega[keep])

    eccentricity = np.asarray(accepted_e[:count])
    omega = np.asarray(accepted_omega[:count])
    e_index = np.abs(eccentricity[:, None] - e_grid).argmin(axis=1)
    omega_index = np.abs(omega[:, None] - omega_grid).argmin(axis=1)
    posterior = np.zeros((count, len(e_grid), len(omega_grid)), dtype=float)
    posterior[np.arange(count), e_index, omega_index] = 1.0
    return posterior


def fitted_sigma_for_mode(
    posterior: np.ndarray,
    e_grid: np.ndarray,
    omega_grid: np.ndarray,
    selection_mode: str,
) -> float:
    mass = np.vstack(
        [
            eccentricity_mass_from_posterior(row, e_grid, omega_grid, selection_mode)
            for row in posterior
        ]
    )
    fit = fit_from_mass_matrix(
        mass,
        e_grid,
        np.linspace(0.05, 0.55, 1001),
        selection_mode=selection_mode,
    )
    return float(fit["sigma_rayleigh"])


def test_forward_normalized_selection_recovers_declared_transiting_population() -> None:
    """Calibrate conventions against simulated data, never a published population value.

    The underlying population is e ~ Rayleigh(0.30), omega ~ Uniform(0, 2pi).
    A planet enters the observed sample with probability proportional to
    (1 + e sin omega) / (1 - e**2), conditional on e < 0.95. The point
    posteriors represent exact observations under a uniform interim e/omega
    measure. This is a generative calibration control, not a reconstruction of
    the manuscript-reciprocal convention.
    """
    e_grid = np.linspace(0.0, E_MAX, 241)
    omega_grid = np.linspace(0.0, 2.0 * np.pi, 72, endpoint=False)
    posterior = selected_rayleigh_point_posteriors(
        count=600,
        seed=20260924,
        e_grid=e_grid,
        omega_grid=omega_grid,
    )

    forward = fitted_sigma_for_mode(
        posterior,
        e_grid,
        omega_grid,
        "legacy_forward_norm",
    )
    reciprocal = fitted_sigma_for_mode(
        posterior,
        e_grid,
        omega_grid,
        "manuscript_reciprocal",
    )

    assert abs(forward - TRUE_SIGMA) < 0.015
    assert abs(reciprocal - TRUE_SIGMA) > 0.04
    assert abs(forward - TRUE_SIGMA) < abs(reciprocal - TRUE_SIGMA)
