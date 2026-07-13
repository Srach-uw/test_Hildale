from __future__ import annotations

import pandas as pd
import pytest

from common import add_target_and_system
from merge_posterior_summaries import annotate_qc
from system_definition_diagnostics import assign_system


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


def test_system_definition_diagnostic_does_not_invent_single_label() -> None:
    sample = pd.DataFrame({"kepid": [100], "system": ["multi"]})
    with pytest.raises(ValueError, match="lacks multiplicity"):
        assign_system(sample, "raw_nonfp", pd.Series({200: "single"}))


def test_unknown_system_completeness_is_qc_excluded() -> None:
    summary = pd.DataFrame(
        {
            "kepoi_name": ["K00100.01"],
            "zeta_median": [1.0],
            "zeta_p16": [0.9],
            "zeta_p84": [1.1],
            "e84": [0.2],
            "period_relative_difference": [0.0],
        }
    )
    audited = annotate_qc(summary)
    assert bool(audited.loc[0, "incomplete_system_unknown"])
    assert bool(audited.loc[0, "incomplete_system"])
    assert bool(audited.loc[0, "qc_primary_exclude"])
    assert "incomplete_system_unknown" in audited.loc[0, "qc_reasons"]
