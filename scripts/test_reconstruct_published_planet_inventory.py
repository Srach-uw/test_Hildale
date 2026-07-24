from __future__ import annotations

import pandas as pd

from reconstruct_published_planet_inventory import (
    add_result_coverage,
    build_pre_visual_inventory,
)


def test_multiplicity_comes_from_pre_cut_koi_count() -> None:
    koi = pd.DataFrame(
        {
            "kepid": [1, 1, 2],
            "kepoi_name": ["K00001.01", "K00001.02", "K00002.01"],
            "koi_disposition": ["CONFIRMED", "CONFIRMED", "CANDIDATE"],
            "koi_period": [10.0, 150.0, 20.0],
            "koi_count": [2, 2, 1],
        }
    )
    hosts = pd.DataFrame(
        {
            "kepid": [1, 2],
            "disk_published": ["thin", "thick"],
            "p_thick_published": [0.1, 0.9],
            "has_measured_velocity": [True, False],
        }
    )
    result = build_pre_visual_inventory(koi, hosts)
    assert result.set_index("kepoi_name").loc["K00001.01", "system"] == "multi"
    assert result.set_index("kepoi_name").loc["K00002.01", "system"] == "single"


def test_result_coverage_prefers_archive() -> None:
    inventory = pd.DataFrame({"koi_target": ["K00001", "K00002", "K00003"]})
    result = add_result_coverage(
        inventory,
        archive_targets={"K00001"},
        cloud_targets={"K00001", "K00002"},
    )
    assert result["preferred_result_source"].tolist() == [
        "original_archive",
        "cloud_missing_run",
        "missing",
    ]
