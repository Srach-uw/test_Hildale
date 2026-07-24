"""Classify ALDERAAN target logs without changing scientific outputs.

The classifier is deliberately conservative. It preserves the raw log path and
stage so a retry can be reviewed against the original method rather than
silently replacing a failed result.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


RULES = (
    ("celerite_linalg", re.compile(r"LinAlgError|failed to factorize|celerite2", re.I)),
    ("low_quality_transits", re.compile(r"over 50% of transits|low quality", re.I)),
    ("limb_darkening", re.compile(r"LD_U1|LD_U2|limb.?dark", re.I)),
    ("missing_lightcurve", re.compile(r"missing_lightcurve|No long-cadence FITS", re.I)),
    ("missing_or_invalid_input", re.compile(r"FileNotFoundError|No such file|not found", re.I)),
    ("posterior_fit", re.compile(r"fit_transit|dynesty|results\.fits", re.I)),
    ("unknown", re.compile(r"Traceback|Error|Exception", re.I)),
)


def classify(text: str, status: str = "") -> str:
    for name, pattern in RULES:
        if pattern.search(text) or pattern.search(status):
            return name
    return "no_error_signature"


def infer_stage(path: Path, text: str) -> str:
    match = re.search(r"failed_([a-z_]+)_exit_\d+", text)
    if match:
        return match.group(1)
    for stage in ("detrend", "noise", "transit_fit", "download"):
        if stage in text.lower() or stage in path.name.lower():
            return stage
    return "unknown"


def target_from_log(path: Path) -> str | None:
    match = re.match(r"(K\d{5})\.(?:stdout|stderr)\.log$", path.name)
    return match.group(1) if match else None


def terminal_error_text(text: str) -> str:
    """Return the final traceback/error region, excluding earlier warnings."""
    traceback_starts = [match.start() for match in re.finditer(r"(?m)^Traceback \(most recent call last\):", text)]
    if traceback_starts:
        return text[traceback_starts[-1] :]
    lines = text.splitlines()
    return "\n".join(lines[-200:])


def find_status(target: str, status_root: Path | None) -> str:
    if status_root is None or not status_root.exists():
        return ""
    matches = sorted(status_root.rglob(f"{target}.status"))
    if not matches:
        return ""
    values = [path.read_text(encoding="utf-8", errors="replace").strip() for path in matches]
    return values[-1]


def has_result(target: str, results_root: Path | None) -> bool:
    if results_root is None or not results_root.exists():
        return False
    return any(path.stat().st_size > 0 for path in results_root.rglob(f"{target}-results.fits"))


def build_taxonomy(
    log_paths: Iterable[Path],
    status_root: Path | None = None,
    results_root: Path | None = None,
) -> pd.DataFrame:
    grouped: dict[str, list[Path]] = {}
    for path in sorted(log_paths):
        if path.is_file() and path.suffix.lower() in {".log", ".txt"}:
            target = target_from_log(path)
            if target is not None:
                grouped.setdefault(target, []).append(path)

    rows = []
    for target, paths in grouped.items():
        status = find_status(target, status_root)
        result_present = has_result(target, results_root)
        texts = [path.read_text(encoding="utf-8", errors="replace") for path in paths]
        combined = "\n".join(texts)
        terminal = terminal_error_text(combined)
        complete = result_present or status == "complete"

        if complete:
            outcome = "success"
            family = "none"
            stage = "complete"
        else:
            # Stage-specific status is authoritative. The final traceback is
            # used for the family so recovered fallback warnings do not become
            # false failures.
            outcome = "failed" if status.startswith("failed_") else "incomplete"
            family = classify(terminal, status)
            stage = infer_stage(paths[-1], status + "\n" + terminal)

        rows.append(
            {
                "target": target,
                "outcome": outcome,
                "status": status,
                "result_present": result_present,
                "stage": stage,
                "failure_family": family,
                "has_traceback": bool(re.search(r"Traceback", terminal, re.I)),
                "logs": "|".join(str(path) for path in paths),
            }
        )

    return pd.DataFrame(
        rows,
        columns=[
            "target",
            "outcome",
            "status",
            "result_present",
            "stage",
            "failure_family",
            "has_traceback",
            "logs",
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--logs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--status-root", type=Path)
    parser.add_argument("--results-root", type=Path)
    args = parser.parse_args()

    frame = build_taxonomy(
        args.logs.rglob("*"),
        status_root=args.status_root,
        results_root=args.results_root,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)
    print(f"classified_targets,{len(frame)}")
    if not frame.empty:
        print(frame.groupby(["outcome", "stage", "failure_family"]).size().to_string())


if __name__ == "__main__":
    main()
