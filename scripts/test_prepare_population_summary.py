from __future__ import annotations

import pandas as pd
import pytest

from prepare_population_summary import attach_context


def test_attach_context_uses_sample_labels_and_intersection() -> None:
    posterior = pd.DataFrame(
        {
            "kepoi_name": ["K1.01", "K2.01"],
            "posterior_file": ["a.npz", "b.npz"],
            "posterior_source": ["direct", "direct"],
            "impact_mode": ["alderaan", "alderaan"],
            "disk": ["wrong", "wrong"],
        }
    )
    sample = pd.DataFrame(
        {
            "kepoi_name": ["K1.01"],
            "koi_target": ["K1"],
            "kepid": [1],
            "disk": ["thin"],
            "system": ["single"],
        }
    )
    got = attach_context(posterior, sample)
    assert got["kepoi_name"].tolist() == ["K1.01"]
    assert got.loc[0, "disk"] == "thin"
    assert got.loc[0, "sample_context_rows_dropped"] == 1


def test_duplicate_sample_context_fails() -> None:
    posterior = pd.DataFrame(
        {"kepoi_name": ["K1.01"], "posterior_file": ["a"], "posterior_source": ["direct"], "impact_mode": ["alderaan"]}
    )
    sample = pd.DataFrame(
        {
            "kepoi_name": ["K1.01", "K1.01"],
            "koi_target": ["K1", "K1"],
            "kepid": [1, 1],
            "disk": ["thin", "thin"],
            "system": ["single", "single"],
        }
    )
    with pytest.raises(ValueError, match="duplicate"):
        attach_context(posterior, sample)
