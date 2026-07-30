import importlib.util
from pathlib import Path

import pytest


def load_module():
    path = Path(__file__).parents[1] / "cloud" / "recovery_142" / "target_pairs.py"
    spec = importlib.util.spec_from_file_location("recovery_target_pairs", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_quoted_standard_csv_header_and_values_are_accepted(tmp_path: Path) -> None:
    targets = tmp_path / "targets.csv"
    targets.write_text(
        '"target_index","koi_target","kepid"\n'
        '"1","K00001","12345"\n',
        encoding="utf-8",
    )
    assert load_module().load_pairs(targets) == [("K00001", "12345")]


def test_duplicate_target_is_rejected(tmp_path: Path) -> None:
    targets = tmp_path / "targets.csv"
    targets.write_text(
        "target_index,koi_target,kepid\n1,K00001,12345\n2,K00001,99999\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate"):
        load_module().load_pairs(targets)
