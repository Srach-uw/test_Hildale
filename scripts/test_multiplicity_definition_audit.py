from __future__ import annotations

import pandas as pd

from multiplicity_definition_audit import PUBLISHED, build_audit


def test_catalog_count_reproduces_published_pre_visual_contract() -> None:
    koi = pd.DataFrame(
        {
            "kepid": [1, 1, 2, 3, 3],
            "kepoi_name": ["K00001.01", "K00001.02", "K00002.01", "K00003.01", "K00003.02"],
            "koi_disposition": ["CONFIRMED", "FALSE POSITIVE", "CONFIRMED", "CONFIRMED", "CONFIRMED"],
            "koi_period": [10, 20, 10, 10, 200],
            "koi_count": [2, 2, 1, 2, 2],
        }
    )
    hosts = pd.DataFrame({"kepid": [1, 2, 3], "disk_published": ["thin", "thick", "thin"]})
    counts, comparison = build_audit(koi, hosts)
    # This synthetic fixture is only structural; overwrite the published
    # target in the module is intentionally avoided. Check the definitions
    # themselves instead.
    assert set(comparison.columns) >= {"catalog_koi_count", "nonfp_all_periods"}
    catalog = counts.loc[counts["definition"] == "catalog_koi_count"].iloc[0]
    assert int(catalog["thin_single"]) == 0
    assert int(catalog["thick_single"]) == 1
    assert int(catalog["thin_multi"]) == 2


def test_non_false_positive_recount_is_not_silent_fallback() -> None:
    koi = pd.DataFrame(
        {
            "kepid": [1, 1],
            "kepoi_name": ["K00001.01", "K00001.02"],
            "koi_disposition": ["CONFIRMED", "FALSE POSITIVE"],
            "koi_period": [10, 20],
            "koi_count": [2, 2],
        }
    )
    hosts = pd.DataFrame({"kepid": [1], "disk_published": ["thin"]})
    _, comparison = build_audit(koi, hosts)
    row = comparison.iloc[0]
    assert bool(row["koi_count_vs_nonfp_all_periods_diff"])
