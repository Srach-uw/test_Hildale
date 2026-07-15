from __future__ import annotations

import numpy as np
import pandas as pd

from population_system_robustness import analyze_population, grouped_log_terms


def test_grouped_log_terms_preserve_planet_likelihood_sum() -> None:
    e_grid = np.linspace(0.001, 0.95, 80)
    sigmas = np.linspace(0.001, 0.5, 100)
    mass = np.zeros((3, len(e_grid)))
    mass[0, 3] = 1.0
    mass[1, 5] = 1.0
    mass[2, 7] = 1.0
    hosts, grouped = grouped_log_terms(
        mass, e_grid, sigmas, pd.Series(["K1", "K1", "K2"]), "manuscript_reciprocal"
    )
    assert hosts.tolist() == ["K1", "K2"]
    assert grouped.shape == (2, len(sigmas))

    summary = pd.DataFrame({"koi_target": ["K1", "K1", "K2"]})
    row, leverage, trials = analyze_population(
        summary,
        mass,
        e_grid,
        sigmas,
        population="thin_multis",
        selection_mode="manuscript_reciprocal",
        trials=20,
        seed=2,
    )
    assert row["n_planets"] == 3
    assert row["n_hosts"] == 2
    assert len(leverage) == 2
    assert len(trials) == 20
