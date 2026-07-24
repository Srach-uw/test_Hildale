from __future__ import annotations

import pandas as pd
import pytest

from targeted_posterior_sensitivity import summarize_deltas


def test_summarize_deltas_groups_population_and_variant() -> None:
    paired = pd.DataFrame(
        {
            "variant": ["a", "a", "a", "b"],
            "population": ["thin_single", "thin_single", "thick_multi", "thin_single"],
            "delta_e50": [0.1, -0.2, 0.3, 0.01],
        }
    )
    out = summarize_deltas(paired).set_index(["variant", "population"])
    assert out.loc[("a", "thin_single"), "n_planets"] == 2
    assert out.loc[("a", "thin_single"), "median_abs_delta_e50"] == pytest.approx(0.15)
    assert out.loc[("b", "thin_single"), "max_abs_delta_e50"] == 0.01
