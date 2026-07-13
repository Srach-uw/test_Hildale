from __future__ import annotations

import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {
    ".bib", ".cff", ".csv", ".json", ".md", ".ps1", ".py", ".sh",
    ".tex", ".toml", ".tsv", ".txt", ".yml", ".yaml",
}
TOOL_NAME_PATTERN = (
    r"\b(?:cl" + r"aude|co" + r"dex|chat" + r"gpt|fa" + r"ble|son" + r"net)\b"
)
TEMP_TOOL_PATTERN = (
    r"(?:AppData[\\/]Local[\\/]Temp[\\/](?:cl" + r"aude|co" + r"dex)"
    r"|Documents[\\/]Co" + r"dex)"
)
FORBIDDEN = {
    "internal tool name": re.compile(TOOL_NAME_PATTERN, re.IGNORECASE),
    "temporary processing path": re.compile(TEMP_TOOL_PATTERN, re.IGNORECASE),
    "personal Windows home path": re.compile(r"C:[\\/]Users[\\/][^\\/\s`]+", re.IGNORECASE),
    "Unicode em dash": re.compile("\N{EM DASH}"),
    "merge conflict marker": re.compile(r"^(?:<{7}|>{7})", re.MULTILINE),
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "likely cloud credential": re.compile(
        r"(?:AIza[0-9A-Za-z_-]{30,}|gh[pousr]_[0-9A-Za-z]{30,}|sk-[0-9A-Za-z]{20,})"
    ),
}


def tracked_files() -> list[Path]:
    output = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT, text=False)
    return [ROOT / item.decode("utf-8") for item in output.split(b"\0") if item]


def scan() -> list[str]:
    findings: list[str] = []
    for path in tracked_files():
        if path == Path(__file__).resolve():
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES or not path.is_file():
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
