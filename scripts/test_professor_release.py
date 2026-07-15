from __future__ import annotations

import check_professor_release as release


def test_release_surface_has_no_forbidden_markers() -> None:
    assert release.scan() == []
