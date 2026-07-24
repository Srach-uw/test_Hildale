from __future__ import annotations

import pandas as pd

from prepare_published_inventory_missing_bundle import classify_retry_targets, sha256


def test_retry_classification_is_explicit() -> None:
    targets = pd.DataFrame({"koi_target": ["K00001", "K00002"]})
    classified = classify_retry_targets(targets, {"K00002"})
    assert classified["prior_launch_status"].tolist() == [
        "not_in_previous_592_target_launch",
        "previously_launched_no_usable_result",
    ]


def test_manifest_hash_is_deterministic(tmp_path) -> None:
    path = tmp_path / "payload.txt"
    path.write_text("same bytes\n", encoding="ascii")
    first = sha256(path)
    second = sha256(path)
    assert first == second
