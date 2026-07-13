from __future__ import annotations

import pandas as pd
import pytest

from common import add_target_and_system


def test_raw_host_multiplicity_survives_planet_level_cuts() -> None:
    raw_nonfp = pd.DataFrame(
        {
            "kepid": [100, 100, 200],
            "kepoi_name": ["K00100.01", "K00100.02", "K00200.01"],
        }
    )
    raw_counts = raw_nonfp.groupby("kepid")["kepoi_name"].size()
    filtered = raw_nonfp.iloc[[0, 2]].copy()

    labeled = add_target_and_system(filtered, raw_counts)

    assert labeled.set_index("kepid").loc[100, "system"] == "multi"
    assert labeled.set_index("kepid").loc[200, "system"] == "single"


def test_filtered_multiplicity_is_explicit_fallback_only() -> None:
    filtered = pd.DataFrame({"kepid": [100], "kepoi_name": ["K00100.01"]})
    labeled = add_target_and_system(filtered)
    assert labeled.loc[0, "system"] == "single"


def test_missing_raw_host_count_fails_closed() -> None:
    filtered = pd.DataFrame({"kepid": [100], "kepoi_name": ["K00100.01"]})
    with pytest.raises(ValueError, match="missing KIC identifiers"):
        add_target_and_system(filtered, pd.Series({200: 1}))
