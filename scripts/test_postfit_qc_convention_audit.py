import numpy as np
import pandas as pd

from postfit_qc_convention_audit import summarize, weighted_fraction


def test_weighted_fraction_normalizes_weights():
    assert weighted_fraction(np.array([True, False]), np.array([2.0, 1.0])) == 2.0 / 3.0


def test_summary_counts_convention_disagreements():
    frame = pd.DataFrame(
        {
            "disk": ["thin", "thin"],
            "system": ["single", "single"],
            "weighted_grazing_fraction": [0.06, 0.01],
            "raw_grazing_fraction": [0.01, 0.06],
            "weighted_grazing_exclude": [True, False],
            "raw_grazing_exclude": [False, True],
            "convention_disagrees": [True, True],
        }
    )
    result = summarize(frame)
    row = result.iloc[0]
    assert row["weighted_excluded"] == 1
    assert row["raw_excluded"] == 1
    assert row["convention_disagreements"] == 2
