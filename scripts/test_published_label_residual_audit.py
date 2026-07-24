from __future__ import annotations

import pandas as pd
import pytest

from published_label_residual_audit import (
    build_migration,
    fraction_true,
    population_key,
    require_columns,
)


def test_population_key_pluralizes_system() -> None:
    assert population_key("thin", "single") == "thin_singles"
    assert population_key("thick", "multi") == "thick_multis"


def test_fraction_true_handles_csv_boolean_strings() -> None:
    values = pd.Series(["True", "false", "1", "0", None])
    assert fraction_true(values) == pytest.approx(0.5)


def test_migration_preserves_system_and_both_disk_labels() -> None:
    summary = pd.DataFrame(
        {
            "kepoi_name": ["K1.01", "K2.01", "K2.02"],
            "kepid": [1, 2, 2],
            "system": ["single", "multi", "multi"],
            "disk_reconstructed": ["thin", "thin", "thin"],
            "disk_published": ["thick", "thin", "thin"],
        }
    )
    result = build_migration(summary)
    migrated = result.loc[
        (result["system"] == "single")
        & (result["disk_reconstructed"] == "thin")
        & (result["disk_published"] == "thick")
    ].iloc[0]
    assert migrated["planets"] == 1
    assert migrated["hosts"] == 1


def test_required_columns_fail_loudly() -> None:
    with pytest.raises(ValueError, match="missing required columns"):
        require_columns(pd.DataFrame({"a": [1]}), {"a", "b"}, "fixture")
