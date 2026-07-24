"""Compare the archived ALDERAAN v0.1.0 tree with the local run tree.

The comparison is semantic for text files: CRLF/LF differences are ignored.
This prevents a version or packaging claim from resting on file timestamps or
platform-specific line endings.
"""
from __future__ import annotations

import argparse
from pathlib import Path


TEXT_SUFFIXES = {".py", ".yml", ".md", ".txt", ".csv", ".gitignore"}


def normalized(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", required=True, type=Path)
    ap.add_argument("--local", required=True, type=Path)
    args = ap.parse_args()
    missing: list[str] = []
    extra: list[str] = []
    different: list[str] = []
    release_files = {p.relative_to(args.release) for p in args.release.rglob("*") if p.is_file()}
    local_files = {p.relative_to(args.local) for p in args.local.rglob("*") if p.is_file()}
    for rel in sorted(release_files | local_files):
        suffix = Path(rel).suffix
        if suffix not in TEXT_SUFFIXES and Path(rel).name != ".gitignore":
            continue
        left = args.release / rel
        right = args.local / rel
        if not left.exists():
            extra.append(str(rel))
        elif not right.exists():
            missing.append(str(rel))
        elif normalized(left) != normalized(right):
            different.append(str(rel))
    print(f"release={args.release}")
    print(f"local={args.local}")
    print(f"semantic_text_files_compared={len(release_files & local_files)}")
    print(f"missing_from_local={len(missing)}")
    print(f"extra_in_local={len(extra)}")
    print(f"semantic_differences={len(different)}")
    for label, values in (("MISSING", missing), ("EXTRA", extra), ("DIFFERENT", different)):
        for value in values:
            print(f"{label}\t{value}")
    if missing or different:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
