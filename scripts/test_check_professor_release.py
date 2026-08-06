from __future__ import annotations

import check_professor_release


def test_release_candidate_files_include_unique_sorted_git_paths(monkeypatch) -> None:
    captured: list[list[str]] = []

    def fake_check_output(command, **kwargs):
        captured.append(command)
        return b"scripts/z.py\0README.md\0scripts/z.py\0"

    monkeypatch.setattr(
        check_professor_release.subprocess,
        "check_output",
        fake_check_output,
    )
    paths = check_professor_release.release_candidate_files()

    assert [path.relative_to(check_professor_release.ROOT).as_posix() for path in paths] == [
        "README.md",
        "scripts/z.py",
    ]
    assert captured == [
        [
            "git",
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
            "-z",
        ]
    ]
