import pandas as pd
import pytest

from summary_inventory_contract_audit import audit_inventory


def frame(names):
    return pd.DataFrame(
        {
            "kepoi_name": names,
            "koi_target": [name.split(".")[0] for name in names],
            "disk": ["thin"] * len(names),
            "system": ["single"] * len(names),
        }
    )


def test_exact_inventory_passes():
    result = audit_inventory(frame(["K00001.01", "K00002.01"]), frame(["K00001.01"]))
    assert result["summary_planets_outside_inventory"] == 0
    assert result["inventory_planets_without_summary"] == 1


def test_mixed_archive_is_visible():
    result = audit_inventory(frame(["K00001.01"]), frame(["K00001.01", "K99999.01"]))
    assert result["summary_planets_outside_inventory"] == 1
    assert result["summary_hosts_outside_inventory"] == 1


def test_duplicate_rows_are_counted():
    inventory = frame(["K00001.01", "K00002.01"])
    summary = frame(["K00001.01", "K00001.01"])
    result = audit_inventory(inventory, summary)
    assert result["summary_duplicate_planet_rows"] == 1
