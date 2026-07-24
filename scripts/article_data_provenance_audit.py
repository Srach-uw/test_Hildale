"""Verify that the local Sagear host table is the journal-hosted article data."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd

from common import load_config, read_sagear2026_kinematic_hosts, root_path
from published_sagear_audit import EXPECTED, validate_published_hosts


EXPECTED_SHA256 = "92280ede0c828413abb7c8314ac6f35b0ccc3e68aabbcf6a94075e84b30e76ed"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit(official_path: Path, canonical_path: Path) -> dict[str, object]:
    official_hash = sha256(official_path)
    canonical_hash = sha256(canonical_path)
    official_bytes = official_path.read_bytes()
    canonical_bytes = canonical_path.read_bytes()
    result: dict[str, object] = {
        "official_path": str(official_path),
        "canonical_path": str(canonical_path),
        "official_sha256": official_hash,
        "canonical_sha256": canonical_hash,
        "official_hash_matches_recorded": official_hash == EXPECTED_SHA256,
        "byte_identical": official_bytes == canonical_bytes,
        "official_bytes": len(official_bytes),
        "canonical_bytes": len(canonical_bytes),
    }
    if not result["official_hash_matches_recorded"]:
        raise AssertionError(f"official article-data hash changed: {official_hash}")
    if not result["byte_identical"]:
        raise AssertionError("canonical host table is not byte-identical to the journal data")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--official", type=Path, required=True)
    parser.add_argument("--canonical", type=Path, default=None)
    args = parser.parse_args()
    cfg = load_config()
    canonical = args.canonical or root_path(cfg, "sagear2026_kinematic_hosts")
    result = audit(args.official.resolve(), canonical.resolve())
    hosts = read_sagear2026_kinematic_hosts(cfg)
    counts = validate_published_hosts(hosts)
    result["host_counts"] = counts
    result["expected_counts"] = EXPECTED
    pd.DataFrame([result]).to_json(
        Path("outputs/article_data_provenance_audit.json"),
        orient="records",
        indent=2,
    )
    print(result)


if __name__ == "__main__":
    main()
