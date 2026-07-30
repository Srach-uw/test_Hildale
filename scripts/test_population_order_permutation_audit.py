import pandas as pd

from population_order_permutation_audit import PUBLISHED, PUBLISHED_ORDER, audit


def test_identity_is_best_for_published_values() -> None:
    frame = pd.DataFrame(
        {
            "population": PUBLISHED_ORDER,
            "expected_e": [PUBLISHED[name][0] for name in PUBLISHED_ORDER],
            "expected_e_lo": [PUBLISHED[name][1] for name in PUBLISHED_ORDER],
            "expected_e_hi": [PUBLISHED[name][2] for name in PUBLISHED_ORDER],
        }
    )
    result = audit(frame)
    assert result.iloc[0]["source_order"] == "|".join(PUBLISHED_ORDER)
    assert result.iloc[0]["interval_overlaps"] == 4
    assert result.iloc[0]["rmse"] == 0.0


def test_all_24_permutations_are_enumerated() -> None:
    frame = pd.DataFrame(
        {
            "population": PUBLISHED_ORDER,
            "expected_e": [0.01, 0.02, 0.03, 0.04],
            "expected_e_lo": [0.0] * 4,
            "expected_e_hi": [1.0] * 4,
        }
    )
    result = audit(frame)
    assert len(result) == 24
    assert result["source_order"].nunique() == 24
