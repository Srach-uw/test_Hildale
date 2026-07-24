from __future__ import annotations

import pandas as pd

from missing_coverage_selection_bias_audit import build_target_level_frame


def test_target_level_coverage_is_system_based() -> None:
    sample = pd.DataFrame(
        {
            "koi_target": ["K00001", "K00001", "K00002"],
            "kepid": [1, 1, 2],
            "disk": ["thin", "thin", "thick"],
            "system": ["multi", "multi", "single"],
            "koi_model_snr": [20.0, 30.0, 10.0],
            "koi_period": [2.0, 5.0, 3.0],
            "koi_prad": [1.0, 2.0, 1.5],
            "koi_impact": [0.2, 0.8, 0.4],
            "P_thick": [0.1, 0.1, 0.9],
            "berger_logg": [4.5, 4.5, 4.4],
            "berger_rad": [1.0, 1.0, 0.9],
            "rho_log": [0.0, 0.0, 0.1],
        }
    )
    missing = pd.DataFrame({"koi_target": ["K00001"]})
    result = build_target_level_frame(sample, missing)
    assert result.set_index("koi_target").loc["K00001", "coverage_state"] == "missing"
    assert result.set_index("koi_target").loc["K00002", "coverage_state"] == "available"
