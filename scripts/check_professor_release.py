from __future__ import annotations

import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {
    ".bib", ".cff", ".csv", ".json", ".md", ".ps1", ".py", ".sh",
    ".tex", ".toml", ".tsv", ".txt", ".yml", ".yaml",
}
FORBIDDEN = {
    "personal Windows home path": re.compile(r"C:[\\/]Users[\\/][^\\/\s`]+", re.IGNORECASE),
    "personal cloud project": re.compile(r"\bproject-[0-9a-f-]{20,}\b", re.IGNORECASE),
    "personal email": re.compile(r"\b[A-Z0-9._%+-]+@gmail\.com\b", re.IGNORECASE),
    "personal VM username": re.compile(r"shreshth" r"_rach1", re.IGNORECASE),
    "Unicode em dash": re.compile("\N{EM DASH}"),
    "merge conflict marker": re.compile(r"^(?:<{7}|>{7})", re.MULTILINE),
    "AI assistant marker": re.compile(
        r"\b(?:as an AI|AI[- ]generated|large language model|ChatGPT|Claude|Codex)\b",
        re.IGNORECASE,
    ),
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "likely cloud credential": re.compile(
        r"(?:AIza[0-9A-Za-z_-]{30,}|gh[pousr]_[0-9A-Za-z]{30,}|sk-[0-9A-Za-z]{20,})"
    ),
}


def release_candidate_files() -> list[Path]:
    output = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        cwd=ROOT,
        text=False,
    )
    names = sorted(
        {item.decode("utf-8") for item in output.split(b"\0") if item}
    )
    return [ROOT / name for name in names]


def scan() -> list[str]:
    findings: list[str] = []
    for path in release_candidate_files():
        if not path.is_file():
            continue
        if path.name.lower() in {"report_hildale.md", "recovery_live_checkpoint.md"}:
            findings.append(f"{path.relative_to(ROOT)}: private project record")
            continue
        if path == Path(__file__).resolve():
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            findings.append(f"{path.relative_to(ROOT)}: not valid UTF-8")
            continue
        for label, pattern in FORBIDDEN.items():
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                findings.append(f"{path.relative_to(ROOT)}:{line}: {label}")
    return findings


def main() -> int:
    findings = scan()
    if findings:
        print("Professor release check failed:")
        for finding in findings:
            print(f"  {finding}")
        return 1
    print("Professor release check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
