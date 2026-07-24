import pandas as pd

from true_rv_subset_audit import build


def test_true_rv_subset_is_host_level(tmp_path):
    summary = pd.DataFrame({"kepid": [1, 1, 2], "disk": ["thin", "thin", "thick"], "system": ["single"] * 3})
    inventory = pd.DataFrame({"kepid": [1, 2], "has_measured_velocity": [True, False]})
    summary_path = tmp_path / "summary.csv"
    inventory_path = tmp_path / "inventory.csv"
    output_path = tmp_path / "out.csv"
    summary.to_csv(summary_path, index=False)
    inventory.to_csv(inventory_path, index=False)
    result = build(summary_path, inventory_path, output_path)
    assert result["kepid"].tolist() == [1, 1]
