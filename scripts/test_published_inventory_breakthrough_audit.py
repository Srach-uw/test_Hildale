from __future__ import annotations

import pandas as pd

from published_inventory_breakthrough_audit import (
    comparison_rows,
    score_hypotheses,
)


def fixture_fit(values: dict[str, float]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "population": population,
                "n": 10,
                "expected_e": value,
                "expected_e_lo": value - 0.005,
                "expected_e_hi": value + 0.005,
                "boundary_flag": False,
            }
            for population, value in values.items()
        ]
    )


def test_exact_values_score_better_than_offset_values() -> None:
    exact = fixture_fit(
        {
            "thick_singles": 0.066,
            "thin_singles": 0.022,
            "thick_multis": 0.033,
            "thin_multis": 0.030,
        }
    )
    offset = exact.copy()
    offset["expected_e"] += 0.05
    offset["expected_e_lo"] += 0.05
    offset["expected_e_hi"] += 0.05
    comparison = pd.concat(
        [
            comparison_rows(exact, "exact"),
            comparison_rows(offset, "offset"),
        ],
        ignore_index=True,
    )
    scores = score_hypotheses(comparison).set_index("hypothesis")
    assert scores.loc["exact", "rmse_expected_e"] == 0
    assert (
        scores.loc["exact", "rmse_expected_e"]
        < scores.loc["offset", "rmse_expected_e"]
    )
