from __future__ import annotations

import pandas as pd
import pytest

from make_disk_swap_diagnostic import swap_disk_labels


def test_swap_is_explicit_and_reversible() -> None:
    source = pd.DataFrame(
        {
            "disk": ["thin", "thick"],
            "system": ["single", "multi"],
            "kepoi_name": ["K00001.01", "K00002.01"],
        }
    )
    swapped = swap_disk_labels(source)
    assert swapped["disk"].tolist() == ["thick", "thin"]
    assert swapped["disk_published_unswapped"].tolist() == ["thin", "thick"]
    assert set(swapped["disk_assignment_mode"]) == {
        "inverted_for_diagnostic_only"
    }
    assert source["disk"].tolist() == ["thin", "thick"]


def test_swap_rejects_unknown_labels() -> None:
    source = pd.DataFrame({"disk": ["halo"], "system": ["single"]})
    with pytest.raises(ValueError, match="thin/thick"):
        swap_disk_labels(source)
