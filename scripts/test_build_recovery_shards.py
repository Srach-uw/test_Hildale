from __future__ import annotations

import pandas as pd
import pytest

from build_recovery_shards import (
    add_runtime_proxy,
    assign_lpt_shards,
    validate_partition,
)


def fixtures() -> tuple[pd.DataFrame, pd.DataFrame]:
    targets = pd.DataFrame(
        {
            "koi_target": ["K00001", "K00002", "K00003", "K00004"],
            "prior_launch_status": [
                "not_in_previous_592_target_launch",
                "previously_launched_no_usable_result",
                "not_in_previous_592_target_launch",
                "not_in_previous_592_target_launch",
            ],
        }
    )
    catalog = pd.DataFrame(
        {
            "koi_id": ["K00001", "K00002", "K00002", "K00003", "K00004"],
            "npl": [1, 2, 2, 1, 1],
            "depth": [1000.0, 100.0, 200.0, 800.0, 900.0],
        }
    )
    return targets, catalog


def test_partition_is_complete_unique_and_deterministic() -> None:
    targets, catalog = fixtures()
    enriched = add_runtime_proxy(targets, catalog)
    first = assign_lpt_shards(enriched, 2)
    second = assign_lpt_shards(enriched, 2)
    validate_partition(targets, first)
    assert first[["koi_target", "shard_id"]].equals(
        second[["koi_target", "shard_id"]]
    )
    assert first.groupby("shard_id")["runtime_proxy"].sum().max() - first.groupby(
        "shard_id"
    )["runtime_proxy"].sum().min() <= enriched["runtime_proxy"].max()


def test_duplicate_target_is_rejected() -> None:
    targets, catalog = fixtures()
    enriched = add_runtime_proxy(targets, catalog)
    duplicate = pd.concat([enriched, enriched.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="unique"):
        assign_lpt_shards(duplicate, 2)


def test_missing_catalog_target_is_rejected() -> None:
    targets, catalog = fixtures()
    with pytest.raises(ValueError, match="missing"):
        add_runtime_proxy(targets, catalog.loc[catalog["koi_id"] != "K00004"])
