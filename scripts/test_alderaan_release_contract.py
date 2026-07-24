from pathlib import Path

from audit_alderaan_zenodo_release import normalized


ROOT = Path(__file__).resolve().parent
RELEASE = ROOT / "inputs" / "sagear_2026_published" / "alderaan-0.1.0_extract" / "alderaan-0.1.0"
LOCAL = ROOT / "external" / "ssagear_alderaan"
TEXT_SUFFIXES = {".py", ".yml", ".md", ".txt", ".csv", ".gitignore"}


def test_published_alderaan_release_matches_local_tree():
    release_files = {
        path.relative_to(RELEASE)
        for path in RELEASE.rglob("*")
        if path.is_file() and (path.suffix in TEXT_SUFFIXES or path.name == ".gitignore")
    }
    local_files = {
        path.relative_to(LOCAL)
        for path in LOCAL.rglob("*")
        if path.is_file() and (path.suffix in TEXT_SUFFIXES or path.name == ".gitignore")
    }
    assert release_files == local_files
    assert all(normalized(RELEASE / rel) == normalized(LOCAL / rel) for rel in release_files)
