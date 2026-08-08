import pandas as pd
import pytest

from sagear_residual_attribution import group_count_rows, require_unique


def test_require_unique_rejects_duplicate_planets():
    frame = pd.DataFrame({"kepoi_name": ["K1.01", "K1.01"]})
    with pytest.raises(ValueError, match="duplicate"):
        require_unique(frame, "kepoi_name", "sample")


def test_group_count_rows_retains_zero_count_categories():
    frame = pd.DataFrame(
        {"disk": ["thin", "thin"], "system": ["single", "single"]}
    )
    rows = group_count_rows(frame, "test", "test radius")
    keyed = {(row["disk"], row["system"]): row for row in rows}
    assert keyed[("thin", "single")]["count"] == 2
    assert keyed[("thick", "multi")]["count"] == 0
