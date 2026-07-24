from pathlib import Path

import pandas as pd

from prepare_published_inventory_missing_bundle import sha256


def test_manifest_hash_and_size_round_trip(tmp_path: Path) -> None:
    payload = tmp_path / "payload.txt"
    payload.write_text("replication\n", encoding="ascii")
    manifest = pd.DataFrame(
        [
            {
                "file": payload.name,
                "bytes": payload.stat().st_size,
                "sha256": sha256(payload),
            }
        ]
    )
    row = manifest.iloc[0]
    assert int(row["bytes"]) == payload.stat().st_size
    assert str(row["sha256"]) == sha256(payload)
