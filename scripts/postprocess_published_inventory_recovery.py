"""Postprocess the exact 142-system published-inventory recovery run.

This is the canonical bridge from a cloud results archive to the completed
Sagear population comparison. It preserves the dynesty-weighted branch and the
equal-raw-row historical diagnostic as separate products and will not run the
final hierarchy while exact-inventory posterior coverage remains incomplete.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
import tarfile
from pathlib import Path

import pandas as pd


BRANCHES = {
    "dynesty": {
        "base": "outputs/eccentricity_posterior_summary_dynesty_published_inventory_pre_visual_qc.csv",
        "allow_non_dynesty": False,
    },
    "equal_raw_diagnostic": {
        "base": "outputs/eccentricity_posterior_summary_equal_nested_published_inventory_pre_visual_qc.csv",
        "allow_non_dynesty": True,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_extract(archive: Path, destination: Path) -> None:
    destination = destination.resolve()
    with tarfile.open(archive, "r:*") as handle:
        for member in handle.getmembers():
            target = (destination / member.name).resolve()
            if destination != target and destination not in target.parents:
                raise ValueError(f"Archive member escapes destination: {member.name}")
        handle.extractall(destination)


def find_results_dir(extracted: Path) -> tuple[Path, int]:
    result_files = sorted(extracted.rglob("*-results.fits"))
    if not result_files:
        raise FileNotFoundError("No ALDERAAN *-results.fits files found in recovery archive")
    result_roots = {
        parent
        for path in result_files
        for parent in path.parents
        if parent.name == "Results"
    }
    if len(result_roots) != 1:
        raise ValueError(f"Expected one Results root, found: {sorted(map(str, result_roots))}")
    return result_roots.pop(), len(result_files)


def run(command: list[str], root: Path) -> None:
    print("+ " + subprocess.list2cmdline(command), flush=True)
    subprocess.run(command, cwd=root, check=True)


def coverage_status(
    sample_path: Path,
    summary_path: Path,
    queue_path: Path,
) -> tuple[pd.DataFrame, int]:
    sample = pd.read_csv(sample_path)
    summary = pd.read_csv(summary_path)
    queue = pd.read_csv(queue_path)
    available = set(summary["kepoi_name"].astype(str))
    rows = []
    for (disk, system), group in sample.groupby(["disk", "system"], dropna=False):
        covered = int(group["kepoi_name"].astype(str).isin(available).sum())
        rows.append(
            {
                "disk": disk,
                "system": system,
                "inventory_planets": len(group),
                "posterior_planets": covered,
                "missing_planets": len(group) - covered,
            }
        )
    table = pd.DataFrame(rows).sort_values(["disk", "system"]).reset_index(drop=True)
    missing = int(table["missing_planets"].sum())
    table["recovery_queue_systems"] = int(queue["koi_target"].nunique())
    return table, missing


def write_status_report(
    path: Path,
    *,
    archive: Path,
    archive_hash: str,
    recovered_fits: int,
    branch_tables: dict[str, pd.DataFrame],
    hierarchy_ran: bool,
) -> None:
    lines = [
        "# Published-inventory recovery status",
        "",
        f"- archive: `{archive}`",
        f"- archive SHA256: `{archive_hash}`",
        f"- recovered ALDERAAN result FITS: {recovered_fits}",
        f"- final hierarchy executed: {'yes' if hierarchy_ran else 'no'}",
        "",
    ]
    for branch, table in branch_tables.items():
        lines.extend(
            [
                f"## {branch}",
                "",
                "| disk | system | inventory | posteriors | missing |",
                "|---|---|---:|---:|---:|",
            ]
        )
        for _, row in table.iterrows():
            lines.append(
                f"| {row['disk']} | {row['system']} | "
                f"{int(row['inventory_planets'])} | "
                f"{int(row['posterior_planets'])} | "
                f"{int(row['missing_planets'])} |"
            )
        lines.append("")
    if not hierarchy_ran:
        lines.extend(
            [
                "The hierarchy was deliberately withheld because exact-inventory",
                "posterior coverage is still incomplete. Inspect the extraction",
                "exclusion manifests and target-level cloud logs before any retry.",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("tmp/published_inventory_recovery_postprocess"),
    )
    parser.add_argument(
        "--sample",
        type=Path,
        default=Path("outputs/sagear2026_planet_inventory_pre_visual_qc.csv"),
    )
    parser.add_argument(
        "--recovery-sample",
        type=Path,
        default=Path(
            "cloud_published_inventory_missing_batch/"
            "published_inventory_missing_population_rows.csv"
        ),
    )
    parser.add_argument(
        "--queue",
        type=Path,
        default=Path(
            "cloud_published_inventory_missing_batch/targets_missing_launchable.csv"
        ),
    )
    parser.add_argument(
        "--run-id",
        default="sagear_published_inventory_missing",
    )
    parser.add_argument("--n-proposals", type=int, default=150_000)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    archive = args.archive.resolve()
    if not archive.is_file():
        raise FileNotFoundError(archive)
    work = (root / args.work_dir).resolve() if not args.work_dir.is_absolute() else args.work_dir.resolve()
    extracted = work / "extracted"
    products = work / "products"
    extracted.mkdir(parents=True, exist_ok=True)
    products.mkdir(parents=True, exist_ok=True)

    archive_hash = sha256(archive)
    marker = extracted / ".archive_sha256"
    if marker.exists():
        prior = marker.read_text(encoding="ascii").strip()
        if prior != archive_hash:
            raise ValueError(
                f"Work directory already contains a different archive: {prior}"
            )
    else:
        safe_extract(archive, extracted)
        marker.write_text(archive_hash + "\n", encoding="ascii")

    results_dir, recovered_fits = find_results_dir(extracted)
    sample = (root / args.sample).resolve()
    recovery_sample = (root / args.recovery_sample).resolve()
    queue = (root / args.queue).resolve()
    for required in (sample, recovery_sample, queue):
        if not required.is_file():
            raise FileNotFoundError(required)

    log_roots = [path for path in extracted.rglob("logs") if path.is_dir()]
    status_roots = [path for path in extracted.rglob("status") if path.is_dir()]
    if len(log_roots) == 1:
        command = [
            sys.executable,
            str(root / "alderaan_failure_taxonomy.py"),
            "--logs",
            str(log_roots[0]),
            "--results-root",
            str(results_root),
            "--output",
            str(products / "alderaan_failure_taxonomy.csv"),
        ]
        if len(status_roots) == 1:
            command.extend(["--status-root", str(status_roots[0])])
        run(command, root)

    branch_tables: dict[str, pd.DataFrame] = {}
    branch_merged: dict[str, Path] = {}
    for branch, config in BRANCHES.items():
        branch_dir = products / branch
        branch_dir.mkdir(parents=True, exist_ok=True)
        summary = branch_dir / "recovered_summary.csv"
        coverage = branch_dir / "recovered_coverage.csv"
        excluded = branch_dir / "recovered_exclusions.csv"
        posterior_dir = branch_dir / "posteriors"
        nested_mode = "equal" if config["allow_non_dynesty"] else "dynesty"
        run(
            [
                sys.executable,
                str(root / "extract_eccentricity_posteriors_direct.py"),
                "--sample",
                str(recovery_sample),
                "--results-dir",
                str(results_dir),
                "--run-id",
                args.run_id,
                "--posterior-subdir",
                str(posterior_dir),
                "--summary-out",
                str(summary),
                "--coverage-out",
                str(coverage),
                "--excluded-out",
                str(excluded),
                "--n-proposals",
                str(args.n_proposals),
                "--nested-weight-mode",
                nested_mode,
                "--posterior-sampling-mode",
                "weighted_grid",
                "--density-source",
                "berger2020_table2",
            ],
            root,
        )

        merged = branch_dir / "exact_inventory_summary.csv"
        qc = branch_dir / "exact_inventory_summary_qc.csv"
        merged_coverage = branch_dir / "exact_inventory_coverage.csv"
        manifest = branch_dir / "exact_inventory_manifest.csv"
        overlap = branch_dir / "overlap_audit.csv"
        base = (root / str(config["base"])).resolve()
        run(
            [
                sys.executable,
                str(root / "assemble_uniform_paired_posteriors.py"),
                "--archive",
                str(base),
                "--new",
                str(summary),
                "--sample",
                str(sample),
                "--new-excluded",
                str(excluded),
                "--out",
                str(merged),
                "--qc-out",
                str(qc),
                "--coverage-out",
                str(merged_coverage),
                "--manifest-out",
                str(manifest),
                "--overlap-out",
                str(overlap),
            ],
            root,
        )
        table, _ = coverage_status(sample, merged, queue)
        table.to_csv(branch_dir / "completion_gate.csv", index=False)
        branch_tables[branch] = table
        branch_merged[branch] = merged

    missing_by_branch = {
        branch: int(table["missing_planets"].sum())
        for branch, table in branch_tables.items()
    }
    hierarchy_ran = all(value == 0 for value in missing_by_branch.values())
    if hierarchy_ran:
        for branch, merged in branch_merged.items():
            command = [
                sys.executable,
                str(root / "hierarchical_rayleigh.py"),
                "--summary",
                str(merged),
                "--selection-mode",
                "manuscript_reciprocal",
                "--diagnostics",
                "--out-tag",
                f"PUBLISHED_RECOVERY_{branch.upper()}",
            ]
            if BRANCHES[branch]["allow_non_dynesty"]:
                command.append("--allow-non-dynesty-weights")
            run(command, root)

    report = products / "published_inventory_recovery_status.md"
    write_status_report(
        report,
        archive=archive,
        archive_hash=archive_hash,
        recovered_fits=recovered_fits,
        branch_tables=branch_tables,
        hierarchy_ran=hierarchy_ran,
    )
    print(f"Recovery report: {report}")
    print(f"Missing planets by branch: {missing_by_branch}")
    if not hierarchy_ran:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
