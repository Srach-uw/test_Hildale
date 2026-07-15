from __future__ import annotations

import pandas as pd

from factorial_validation_diagnostics import (
    clean_pairs,
    cluster_bootstrap,
    leverage_table,
    population_summaries,
)


def synthetic_pairs() -> pd.DataFrame:
    rows = []
    for index, delta in enumerate([0.01, -0.02, 0.03, 0.00, 0.02, -0.01]):
        rows.append(
            {
                "comparison_id": "test",
                "effect": "limb_darkening",
                "koi_target": f"K{index // 2:05d}",
                "kepoi_name": f"K{index:05d}.01",
                "disk_baseline": "thin" if index < 4 else "thick",
                "system_baseline": "single",
                "e16_baseline": 0.1,
                "e50_baseline": 0.2,
                "e84_baseline": 0.4,
                "e16_comparison": 0.1,
                "e50_comparison": 0.2 + delta,
                "e84_comparison": 0.4,
                "impact_p50_baseline": 0.4,
                "impact_p50_comparison": 0.42,
                "delta_e": delta,
                "delta_impact": 0.02,
                "delta_t14_hr": 0.01,
                "fractional_delta_rp_over_rs": 0.01,
                "qc_primary_exclude_baseline": False,
                "qc_primary_exclude_comparison": False,
            }
        )
    return pd.DataFrame(rows)


def test_population_summary_and_leverage_are_system_clustered() -> None:
    pairs = clean_pairs(synthetic_pairs())
    summary = population_summaries(pairs, n_bootstrap=100, seed=7)
    thin = summary[summary["population"].eq("thin_single")].iloc[0]
    assert thin["n_planets"] == 4
    assert thin["n_systems"] == 2
    leverage = leverage_table(pairs)
    assert set(leverage["koi_target"]) == {"K00000", "K00001", "K00002"}
    assert leverage["leave_system_out_median_delta_e"].notna().all()


def test_cluster_bootstrap_is_seed_deterministic_and_preserves_system_rows() -> None:
    pairs = clean_pairs(synthetic_pairs())
    first = cluster_bootstrap(pairs, "delta_e", n_bootstrap=200, seed=11)
    second = cluster_bootstrap(pairs, "delta_e", n_bootstrap=200, seed=11)
    assert first == second
    assert first[0] <= pairs["delta_e"].median() <= first[1]
