from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

from alderaan_batch import build_alderaan_catalog
from common import load_config, output_dir, root_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare a cloud-ready ALDERAAN bundle for launchable missing-posterior systems.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--queue", default=None, help="Planet-level missing launchable queue CSV.")
    parser.add_argument("--sample", default=None, help="Canonical sample CSV.")
    parser.add_argument("--out", default=None, help="Output bundle directory.")
    parser.add_argument("--run-id", default="sagear_missing")
    parser.add_argument("--jobs", type=int, default=30)
    parser.add_argument("--shards", type=int, default=4)
    parser.add_argument("--copy-to-external-output", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = Path(args.out) if args.out else Path(__file__).resolve().parent / "cloud_missing_batch"
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = output_dir()
    queue_path = Path(args.queue) if args.queue else outputs / "posterior_queue_missing_launchable_alderaan_best.csv"
    sample_path = Path(args.sample) if args.sample else outputs / "canonical_sample_old_astropy_rawcc.csv"
    queue = pd.read_csv(queue_path)
    sample = pd.read_csv(sample_path)

    targets = build_targets(queue, sample)
    expanded = sample[sample["koi_target"].isin(targets["koi_target"])].copy()
    catalog = build_alderaan_catalog(expanded, cfg)

    targets_path = out_dir / "targets_missing_launchable.csv"
    catalog_path = out_dir / "sagear_missing_catalog.csv"
    expanded_path = out_dir / "catalog_rows_missing_launchable.csv"
    manifest_path = out_dir / "cloud_missing_manifest.md"
    validate_path = out_dir / "validate_bundle.py"

    targets.to_csv(targets_path, index=False)
    catalog.to_csv(catalog_path)
    expanded.to_csv(expanded_path, index=False)
    write_shards(out_dir, targets, args.shards)
    write_scripts(out_dir, args.run_id, args.jobs)
    write_validator(validate_path)
    write_readme(out_dir, args.run_id, args.jobs, targets, catalog)
    write_manifest(manifest_path, args.run_id, args.jobs, targets, catalog, queue)

    copy_root = Path(args.copy_to_external_output) if args.copy_to_external_output else None
    if copy_root:
        copy_root.mkdir(parents=True, exist_ok=True)
        for path in [targets_path, catalog_path, expanded_path, manifest_path]:
            (copy_root / path.name).write_bytes(path.read_bytes())

    print(f"Wrote cloud missing-posterior bundle: {out_dir}")
    print(f"Missing launchable planet rows: {len(queue)}")
    print(f"Target systems: {len(targets)}")
    print(f"ALDERAAN catalog rows after system expansion: {len(catalog)}")
    print(targets.groupby(["disk", "system"]).agg(targets=("koi_target", "count"), missing_planets=("missing_planets", "sum")).to_string())


def build_targets(queue: pd.DataFrame, sample: pd.DataFrame) -> pd.DataFrame:
    q = (
        queue.groupby("koi_target")
        .agg(
            kepid=("kepid", "first"),
            disk=("disk", lambda x: ",".join(sorted(set(map(str, x))))),
            system=("system", lambda x: ",".join(sorted(set(map(str, x))))),
            P_thick=("P_thick", "median"),
            missing_planets=("kepoi_name", "count"),
            missing_kepoi_names=("kepoi_name", lambda x: ",".join(map(str, x))),
            min_missing_period=("koi_period", "min"),
            max_missing_period=("koi_period", "max"),
            max_missing_snr=("koi_model_snr", "max"),
        )
        .reset_index()
    )
    meta = (
        sample.groupby("koi_target")
        .agg(
            target_sample_planets=("kepoi_name", "count"),
            all_kepoi_names=("kepoi_name", lambda x: ",".join(map(str, x))),
            min_period=("koi_period", "min"),
            max_period=("koi_period", "max"),
            max_snr=("koi_model_snr", "max"),
        )
        .reset_index()
    )
    out = q.merge(meta, on="koi_target", how="left")
    out = out.sort_values(["disk", "system", "koi_target"]).reset_index(drop=True)
    out.insert(0, "target_index", range(len(out)))
    return out


def write_shards(out_dir: Path, targets: pd.DataFrame, shards: int) -> None:
    shard_dir = out_dir / "shards"
    shard_dir.mkdir(exist_ok=True)
    shards = max(1, int(shards))
    n = len(targets)
    size = max(1, math.ceil(n / shards))
    summary = []
    for i in range(shards):
        sub = targets.iloc[i * size : (i + 1) * size].copy()
        if sub.empty:
            continue
        path = shard_dir / f"targets_shard_{i:03d}.csv"
        sub.to_csv(path, index=False)
        summary.append(
            {
                "shard": i,
                "file": str(path.name),
                "targets": len(sub),
                "missing_planets": int(sub["missing_planets"].sum()),
            }
        )
    pd.DataFrame(summary).to_csv(out_dir / "shard_summary.csv", index=False)


def write_scripts(out_dir: Path, run_id: str, jobs: int) -> None:
    (out_dir / "setup_vm.sh").write_text(
        """#!/usr/bin/env bash
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
sudo apt-get update
sudo apt-get install -y git curl wget bzip2 build-essential parallel rsync

if [ ! -d "$HOME/miniforge3" ]; then
  curl -L -o /tmp/miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
fi

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda config --set channel_priority flexible

if [ ! -d "$HOME/alderaan" ]; then
  git clone https://github.com/gjgilbert/alderaan "$HOME/alderaan"
fi

if ! conda env list | awk '{print $1}' | grep -qx alderaan; then
  if command -v mamba >/dev/null 2>&1; then
    mamba env create -n alderaan -f "$HOME/alderaan/environment.yml"
  else
    conda env create -n alderaan -f "$HOME/alderaan/environment.yml"
  fi
fi

echo "Ready. Next:"
echo "source ~/miniforge3/etc/profile.d/conda.sh"
echo "conda activate alderaan"
echo "bash validate_bundle.py"
echo "JOBS=30 bash run_batch.sh"
""",
        encoding="utf-8",
    )

    (out_dir / "run_one_target.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail

export OMP_NUM_THREADS="${{OMP_NUM_THREADS:-1}}"
export MKL_NUM_THREADS="${{MKL_NUM_THREADS:-1}}"
export OPENBLAS_NUM_THREADS="${{OPENBLAS_NUM_THREADS:-1}}"
export NUMEXPR_NUM_THREADS="${{NUMEXPR_NUM_THREADS:-1}}"

TARGET="${{1:?target required}}"
KEPID="${{2:?kepid required}}"
PROJECT_DIR="${{PROJECT_DIR:-$PWD/alderaan_project}}"
RUN_ID="${{RUN_ID:-{run_id}}}"
MISSION="${{MISSION:-Kepler}}"
ALDERAAN_REPO="${{ALDERAAN_REPO:-$HOME/alderaan}}"
DATA_DIR="$PROJECT_DIR/Data"
CATALOG_NAME="${{CATALOG_NAME:-sagear_missing_catalog.csv}}"
RESULTS_FITS="$PROJECT_DIR/Results/$RUN_ID/$TARGET/$TARGET-results.fits"
STATUS_DIR="$PROJECT_DIR/status"

mkdir -p "$DATA_DIR" "$STATUS_DIR" logs

if [ -s "$RESULTS_FITS" ]; then
  echo "complete" > "$STATUS_DIR/$TARGET.status"
  echo "$TARGET already complete"
  exit 0
fi

echo "running" > "$STATUS_DIR/$TARGET.status"
echo "[$(date -Is)] Starting $TARGET / KIC $KEPID"

pushd "$DATA_DIR" >/dev/null
python "$ALDERAAN_REPO/bin/get_kepler_data.py" "$KEPID" -c long -t lightcurve -o "get_${{KEPID}}_lc.sh" --cmdtype wget
bash "get_${{KEPID}}_lc.sh" || true
popd >/dev/null

pushd "$ALDERAAN_REPO" >/dev/null
python bin/detrend_and_estimate_ttvs.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR" --data_dir "$DATA_DIR/" --catalog "$CATALOG_NAME"
python bin/analyze_autocorrelated_noise.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR" --data_dir "$DATA_DIR/" --catalog "$CATALOG_NAME"
python bin/fit_transit_shape_simultaneous_nested.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR"
popd >/dev/null

if [ -s "$RESULTS_FITS" ]; then
  echo "complete" > "$STATUS_DIR/$TARGET.status"
  echo "[$(date -Is)] Finished $TARGET"
else
  echo "missing_results" > "$STATUS_DIR/$TARGET.status"
  echo "[$(date -Is)] Finished without results FITS: $TARGET" >&2
  exit 3
fi
""",
        encoding="utf-8",
    )

    (out_dir / "run_batch.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail

export OMP_NUM_THREADS="${{OMP_NUM_THREADS:-1}}"
export MKL_NUM_THREADS="${{MKL_NUM_THREADS:-1}}"
export OPENBLAS_NUM_THREADS="${{OPENBLAS_NUM_THREADS:-1}}"
export NUMEXPR_NUM_THREADS="${{NUMEXPR_NUM_THREADS:-1}}"

JOBS="${{JOBS:-{jobs}}}"
RUN_ID="${{RUN_ID:-{run_id}}}"
TARGET_CSV="${{TARGET_CSV:-targets_missing_launchable.csv}}"
PROJECT_DIR="${{PROJECT_DIR:-$PWD/alderaan_project}}"
ALDERAAN_REPO="${{ALDERAAN_REPO:-$HOME/alderaan}}"

mkdir -p "$PROJECT_DIR/Catalogs" "$PROJECT_DIR/Data" "$PROJECT_DIR/Results" "$PROJECT_DIR/Figures" "$PROJECT_DIR/status" logs
cp sagear_missing_catalog.csv "$PROJECT_DIR/Catalogs/sagear_missing_catalog.csv"
cp "$TARGET_CSV" "$PROJECT_DIR/targets.csv"
chmod +x run_one_target.sh

echo "Running $TARGET_CSV with JOBS=$JOBS RUN_ID=$RUN_ID"
python validate_bundle.py --targets "$TARGET_CSV" --catalog sagear_missing_catalog.csv

tail -n +2 "$TARGET_CSV" | awk -F, '{{print $2","$3}}' \\
  | parallel --colsep ',' -j "$JOBS" --joblog "logs/parallel_joblog_${{RUN_ID}}.tsv" './run_one_target.sh {{1}} {{2}}' \\
  > "logs/batch_${{RUN_ID}}_stdout.log" 2> "logs/batch_${{RUN_ID}}_stderr.log"

bash summarize_progress.sh
echo "Batch complete. Results are under $PROJECT_DIR/Results/$RUN_ID"
""",
        encoding="utf-8",
    )

    (out_dir / "summarize_progress.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
RUN_ID="${{RUN_ID:-{run_id}}}"
PROJECT_DIR="${{PROJECT_DIR:-$PWD/alderaan_project}}"
TARGET_CSV="${{TARGET_CSV:-targets_missing_launchable.csv}}"
TOTAL=$(($(wc -l < "$TARGET_CSV") - 1))
COMPLETE=$(find "$PROJECT_DIR/Results/$RUN_ID" -name '*-results.fits' 2>/dev/null | wc -l || true)
FAILED=$(find "$PROJECT_DIR/status" -name '*.status' -print0 2>/dev/null | xargs -0 grep -l -v '^complete$' 2>/dev/null | wc -l || true)
echo "targets_total,$TOTAL" > progress_summary.csv
echo "results_fits,$COMPLETE" >> progress_summary.csv
echo "failed_or_incomplete_status,$FAILED" >> progress_summary.csv
cat progress_summary.csv
""",
        encoding="utf-8",
    )

    (out_dir / "pack_results.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
RUN_ID="${{RUN_ID:-{run_id}}}"
PROJECT_DIR="${{PROJECT_DIR:-$PWD/alderaan_project}}"
OUT="alderaan_results_${{RUN_ID}}_$(date +%Y%m%d_%H%M%S).tar.gz"
bash summarize_progress.sh
tar -czf "$OUT" -C "$PROJECT_DIR" "Results/$RUN_ID" status -C "$PWD" logs progress_summary.csv targets_missing_launchable.csv sagear_missing_catalog.csv cloud_missing_manifest.md
echo "$OUT"
""",
        encoding="utf-8",
    )

    (out_dir / "create_gcp_spot_vm.sh").write_text(
        """#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_ID:?Set PROJECT_ID first}"
ZONE="${ZONE:-us-central1-a}"
INSTANCE="${INSTANCE:-alderaan-missing-e2-32}"
MACHINE="${MACHINE:-e2-standard-32}"
DISK_SIZE="${DISK_SIZE:-150GB}"
IMAGE_FAMILY="${IMAGE_FAMILY:-ubuntu-2204-lts}"
IMAGE_PROJECT="${IMAGE_PROJECT:-ubuntu-os-cloud}"

gcloud config set project "$PROJECT_ID"
gcloud compute instances create "$INSTANCE" \
  --zone "$ZONE" \
  --machine-type "$MACHINE" \
  --provisioning-model=SPOT \
  --instance-termination-action=STOP \
  --boot-disk-size "$DISK_SIZE" \
  --boot-disk-type pd-ssd \
  --image-family "$IMAGE_FAMILY" \
  --image-project "$IMAGE_PROJECT" \
  --metadata=google-logging-enabled=true

echo "Created $INSTANCE in $ZONE"
echo "Copy bundle with:"
echo "gcloud compute scp --recurse . $INSTANCE:~/sagear_cloud_missing --zone $ZONE"
""",
        encoding="utf-8",
    )

    for script in [
        "setup_vm.sh",
        "run_one_target.sh",
        "run_batch.sh",
        "summarize_progress.sh",
        "pack_results.sh",
        "create_gcp_spot_vm.sh",
    ]:
        (out_dir / script).chmod(0o755)


def write_validator(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", default="targets_missing_launchable.csv")
    parser.add_argument("--catalog", default="sagear_missing_catalog.csv")
    args = parser.parse_args()

    targets = pd.read_csv(args.targets)
    catalog = pd.read_csv(args.catalog, index_col=0)
    failures = []
    if targets["koi_target"].duplicated().any():
        failures.append("duplicate koi_target rows in targets")
    missing_catalog_targets = set(targets["koi_target"]) - set(catalog["koi_id"])
    if missing_catalog_targets:
        failures.append(f"{len(missing_catalog_targets)} target(s) missing from catalog")
    for col in ["koi_id", "kic_id", "npl", "period", "epoch", "depth", "duration", "impact"]:
        if col not in catalog.columns:
            failures.append(f"catalog missing column {col}")
        elif catalog[col].isna().any():
            failures.append(f"catalog column {col} has NaNs")
    for col in ["period", "depth", "duration"]:
        if col in catalog.columns and (pd.to_numeric(catalog[col], errors="coerce") <= 0).any():
            failures.append(f"catalog column {col} has non-positive values")
    npl = catalog.groupby("koi_id").size().rename("actual").reset_index()
    decl = catalog.groupby("koi_id")["npl"].first().rename("declared").reset_index()
    bad = npl.merge(decl, on="koi_id")
    bad = bad[bad["actual"] != bad["declared"]]
    if len(bad):
        failures.append(f"{len(bad)} target(s) have npl mismatch")
    if failures:
        print("VALIDATION FAILED")
        for failure in failures:
            print("-", failure)
        sys.exit(2)
    print(f"VALIDATION OK: {len(targets)} targets, {len(catalog)} catalog rows")


if __name__ == "__main__":
    main()
""",
        encoding="utf-8",
    )
    path.chmod(0o755)


def write_readme(out_dir: Path, run_id: str, jobs: int, targets: pd.DataFrame, catalog: pd.DataFrame) -> None:
    readme = f"""# GCP ALDERAAN Missing-Posterior Bundle

This bundle is for the true missing-posterior queue only.

- Run id: `{run_id}`
- Runnable target systems: `{len(targets)}`
- Missing planet rows covered by those targets: `{int(targets['missing_planets'].sum())}`
- ALDERAAN catalog rows after whole-system expansion: `{len(catalog)}`
- Default parallel jobs: `{jobs}`

## Local Sanity Check

```bash
python validate_bundle.py
```

## Create VM

```bash
export PROJECT_ID=your-gcp-project
export ZONE=us-central1-a
bash create_gcp_spot_vm.sh
gcloud compute scp --recurse . alderaan-missing-e2-32:~/sagear_cloud_missing --zone $ZONE
gcloud compute ssh alderaan-missing-e2-32 --zone $ZONE
```

## Run On VM

```bash
cd ~/sagear_cloud_missing
bash setup_vm.sh
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alderaan
python validate_bundle.py
JOBS={jobs} bash run_batch.sh
bash pack_results.sh
```

For a test slice:

```bash
TARGET_CSV=shards/targets_shard_000.csv JOBS=8 bash run_batch.sh
```

## Retrieve

```bash
gcloud compute scp alderaan-missing-e2-32:~/sagear_cloud_missing/alderaan_results_{run_id}_*.tar.gz . --zone $ZONE
```

Extract locally into `sagear_reproduction/alderaan_project` or a separate results folder, then run:

```powershell
python sagear_reproduction\\extract_eccentricity_posteriors.py `
  --sample sagear_reproduction\\outputs\\canonical_sample_old_astropy_rawcc.csv `
  --run-id {run_id} `
  --posterior-subdir eccentricity_posteriors_{run_id} `
  --summary-out sagear_reproduction\\outputs\\eccentricity_posterior_summary_{run_id}.csv `
  --coverage-out sagear_reproduction\\outputs\\eccentricity_posterior_coverage_{run_id}.csv

python sagear_reproduction\\merge_posterior_summaries.py `
  --new sagear_reproduction\\outputs\\eccentricity_posterior_summary_{run_id}.csv `
  --out sagear_reproduction\\outputs\\eccentricity_posterior_summary_merged_{run_id}.csv `
  --coverage-out sagear_reproduction\\outputs\\eccentricity_posterior_coverage_merged_{run_id}.csv
```
"""
    (out_dir / "README_GCP_MISSING.md").write_text(readme, encoding="utf-8")


def write_manifest(path: Path, run_id: str, jobs: int, targets: pd.DataFrame, catalog: pd.DataFrame, queue: pd.DataFrame) -> None:
    pop = (
        queue.groupby(["disk", "system"])
        .size()
        .rename("missing_planets")
        .reset_index()
        .sort_values(["disk", "system"])
    )
    target_pop = (
        targets.groupby(["disk", "system"])
        .agg(targets=("koi_target", "count"), missing_planets=("missing_planets", "sum"))
        .reset_index()
        .sort_values(["disk", "system"])
    )
    lines = [
        "# Cloud Missing-Posterior ALDERAAN Manifest",
        "",
        f"- Run id: `{run_id}`",
        f"- Default jobs: `{jobs}`",
        f"- Missing launchable planet rows: `{len(queue)}`",
        f"- Target systems: `{len(targets)}`",
        f"- Catalog rows after system expansion: `{len(catalog)}`",
        "",
        "## Missing Planets By Population",
        "",
        pop.to_markdown(index=False),
        "",
        "## Target Systems By Population",
        "",
        target_pop.to_markdown(index=False),
        "",
        "## First 30 Targets",
        "",
        targets.head(30).to_markdown(index=False, floatfmt=".4f"),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
