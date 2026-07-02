from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd

from common import load_config, output_dir, root_path
from alderaan_batch import build_alderaan_catalog


PIPELINE_DIR = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare a self-contained GCP ALDERAAN batch bundle.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None, help="classified canonical sample path")
    parser.add_argument("--out", default=None, help="Output cloud bundle directory")
    parser.add_argument("--max-planets", type=int, default=300)
    parser.add_argument("--jobs", type=int, default=30)
    parser.add_argument(
        "--strategy",
        choices=["stratified", "thick_first", "all"],
        default="stratified",
        help="How to choose targets from canonical_sample.",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    sample_path = Path(args.sample) if args.sample else output_dir() / "canonical_sample_diagnostic.csv"
    if not sample_path.exists():
        raise FileNotFoundError(
            "Missing classified sample for cloud target selection. "
            f"Run diagnose_sample.py --allow-fallback-gmm first, or pass --sample explicitly: {sample_path}"
        )

    out_dir = Path(args.out) if args.out else PIPELINE_DIR / "cloud_batch"
    out_dir.mkdir(parents=True, exist_ok=True)

    sample = pd.read_csv(sample_path)
    targets = choose_targets(sample, args.max_planets, args.strategy)
    selected = sample[sample["koi_target"].isin(targets["koi_target"])]
    catalog = build_alderaan_catalog(selected, cfg)

    targets.to_csv(out_dir / "targets.csv", index=False)
    catalog.to_csv(out_dir / "sagear_cloud_catalog.csv")
    write_scripts(out_dir, jobs=args.jobs, cfg=cfg)
    write_readme(out_dir, args.jobs, args.max_planets)

    print(f"Wrote cloud batch bundle: {out_dir}")
    print(f"Targets: {len(targets)} systems, {int(targets['n_planets'].sum())} planet rows")
    print(targets.groupby(["disk", "system"]).agg(systems=("koi_target", "count"), planets=("n_planets", "sum")))


def choose_targets(sample: pd.DataFrame, max_planets: int, strategy: str) -> pd.DataFrame:
    required = ["koi_period", "koi_time0bk", "koi_depth", "koi_duration", "koi_impact"]
    df = sample.copy()
    for col in required:
        df = df[np.isfinite(pd.to_numeric(df[col], errors="coerce"))]
    df = df[(df["koi_depth"] > 0) & (df["koi_duration"] > 0)]

    systems = (
        df.groupby("koi_target")
        .agg(
            kepid=("kepid", "first"),
            disk=("disk", "first"),
            system=("system", "first"),
            p_thick=("P_thick", "median"),
            n_planets=("kepoi_name", "count"),
            min_period=("koi_period", "min"),
        )
        .reset_index()
    )
    systems = systems[systems["disk"].isin(["thin", "thick"])]

    if strategy == "all":
        ordered = systems.sort_values(["disk", "system", "kepid"])
    elif strategy == "thick_first":
        ordered = systems.assign(rank_disk=np.where(systems["disk"] == "thick", 0, 1)).sort_values(
            ["rank_disk", "system", "kepid"]
        )
    else:
        pieces = []
        bins = [
            ("thick", "single"),
            ("thin", "single"),
            ("thick", "multi"),
            ("thin", "multi"),
        ]
        quota = max(1, max_planets // len(bins))
        for disk, system in bins:
            sub = systems[(systems["disk"] == disk) & (systems["system"] == system)].copy()
            sub = sub.sort_values(["p_thick", "kepid"], ascending=[disk == "thin", True])
            pieces.append(take_until_planets(sub, quota))
        ordered = pd.concat(pieces, ignore_index=True)
        remainder = systems[~systems["koi_target"].isin(ordered["koi_target"])].sort_values(["disk", "system", "kepid"])
        ordered = pd.concat([ordered, remainder], ignore_index=True)

    return take_until_planets(ordered, max_planets).drop(columns=[c for c in ["rank_disk"] if c in ordered])


def take_until_planets(df: pd.DataFrame, max_planets: int) -> pd.DataFrame:
    rows = []
    total = 0
    for _, row in df.iterrows():
        n = int(row["n_planets"])
        if rows and total + n > max_planets:
            continue
        rows.append(row)
        total += n
        if total >= max_planets:
            break
    return pd.DataFrame(rows)


def write_scripts(out_dir: Path, jobs: int, cfg: dict) -> None:
    (out_dir / "setup_vm.sh").write_text(
        """#!/usr/bin/env bash
set -euo pipefail

sudo apt-get update
sudo apt-get install -y git curl wget bzip2 build-essential parallel

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
  mamba env create -n alderaan -f "$HOME/alderaan/environment.yml" || conda env create -n alderaan -f "$HOME/alderaan/environment.yml"
fi

echo "ALDERAAN env ready. Activate with: conda activate alderaan"
""",
        encoding="utf-8",
    )

    (out_dir / "run_one_target.sh").write_text(
        """#!/usr/bin/env bash
set -euo pipefail

TARGET="${1:?target required}"
KEPID="${2:?kepid required}"
PROJECT_DIR="${PROJECT_DIR:-$PWD/alderaan_project}"
RUN_ID="${RUN_ID:-sagear_cloud}"
MISSION="${MISSION:-Kepler}"
ALDERAAN_REPO="${ALDERAAN_REPO:-$HOME/alderaan}"
DATA_DIR="$PROJECT_DIR/Data"
CATALOG_NAME="${CATALOG_NAME:-sagear_cloud_catalog.csv}"
RESULTS_FITS="$PROJECT_DIR/Results/$RUN_ID/$TARGET/$TARGET-results.fits"
STATUS_DIR="$PROJECT_DIR/status"

mkdir -p "$DATA_DIR" "$STATUS_DIR"

if [ -f "$RESULTS_FITS" ]; then
  echo "$TARGET already complete: $RESULTS_FITS"
  echo "complete" > "$STATUS_DIR/$TARGET.status"
  exit 0
fi

echo "[$(date -Is)] Starting $TARGET / KIC $KEPID"
pushd "$DATA_DIR" >/dev/null
python "$ALDERAAN_REPO/bin/get_kepler_data.py" "$KEPID" -c long -t lightcurve -o "get_${KEPID}_lc.sh" --cmdtype curl
bash "get_${KEPID}_lc.sh" || true
popd >/dev/null

pushd "$ALDERAAN_REPO" >/dev/null
python bin/detrend_and_estimate_ttvs.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR" --data_dir "$DATA_DIR/" --catalog "$CATALOG_NAME"
python bin/analyze_autocorrelated_noise.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR"
python bin/fit_transit_shape_simultaneous_nested.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR"
popd >/dev/null

echo "complete" > "$STATUS_DIR/$TARGET.status"
echo "[$(date -Is)] Finished $TARGET"
""",
        encoding="utf-8",
    )

    (out_dir / "run_batch.sh").write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail

JOBS="${{JOBS:-{jobs}}}"
RUN_ID="${{RUN_ID:-sagear_cloud}}"
PROJECT_DIR="${{PROJECT_DIR:-$PWD/alderaan_project}}"
ALDERAAN_REPO="${{ALDERAAN_REPO:-$HOME/alderaan}}"

mkdir -p "$PROJECT_DIR/Catalogs" "$PROJECT_DIR/Data" "$PROJECT_DIR/Results" "$PROJECT_DIR/Figures" "$PROJECT_DIR/status" logs
cp sagear_cloud_catalog.csv "$PROJECT_DIR/Catalogs/sagear_cloud_catalog.csv"
cp targets.csv "$PROJECT_DIR/targets.csv"

chmod +x run_one_target.sh

tail -n +2 targets.csv | awk -F, '{{print $1","$2}}' | parallel --colsep ',' -j "$JOBS" --joblog logs/parallel_joblog.tsv './run_one_target.sh {{1}} {{2}}' > logs/batch_stdout.log 2> logs/batch_stderr.log

echo "Batch complete. Results are under $PROJECT_DIR/Results/$RUN_ID"
""",
        encoding="utf-8",
    )

    (out_dir / "create_gcp_spot_vm.sh").write_text(
        """#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_ID:?Set PROJECT_ID first}"
ZONE="${ZONE:-us-central1-a}"
INSTANCE="${INSTANCE:-alderaan-e2-32}"
MACHINE="${MACHINE:-e2-standard-32}"
DISK_SIZE="${DISK_SIZE:-100GB}"

gcloud config set project "$PROJECT_ID"
gcloud compute instances create "$INSTANCE" \
  --zone "$ZONE" \
  --machine-type "$MACHINE" \
  --provisioning-model=SPOT \
  --instance-termination-action=STOP \
  --boot-disk-size "$DISK_SIZE" \
  --boot-disk-type pd-ssd \
  --image-family ubuntu-2204-lts \
  --image-project ubuntu-os-cloud \
  --maintenance-policy TERMINATE \
  --scopes cloud-platform

echo "VM created. Copy this folder up with:"
echo "gcloud compute scp --recurse . $INSTANCE:~/sagear_cloud_batch --zone $ZONE"
""",
        encoding="utf-8",
    )

    for script in ["setup_vm.sh", "run_one_target.sh", "run_batch.sh", "create_gcp_spot_vm.sh"]:
        (out_dir / script).chmod(0o755)


def write_readme(out_dir: Path, jobs: int, max_planets: int) -> None:
    (out_dir / "README_GCP.md").write_text(
        f"""# GCP ALDERAAN Batch Bundle

This bundle runs ALDERAAN by KOI system target, not by planet row.

Generated target budget: up to {max_planets} planet rows.
Default parallelism: {jobs} target jobs.

## Critical Notes

- ALDERAAN's unit of work is a KOI system target such as `K00085`; multi-planet systems are fit simultaneously.
- Spot VMs are cheaper but can be preempted, so this runner skips targets that already have `TARGET-results.fits`.
- Pricing changes by region and time. Check current Google pages before spending real money:
  - https://cloud.google.com/products/compute/pricing
  - https://cloud.google.com/spot-vms/pricing
  - https://cloud.google.com/free

## VM Setup

```bash
export PROJECT_ID=your-gcp-project
export ZONE=us-central1-a
./create_gcp_spot_vm.sh
gcloud compute scp --recurse . alderaan-e2-32:~/sagear_cloud_batch --zone $ZONE
gcloud compute ssh alderaan-e2-32 --zone $ZONE
```

On the VM:

```bash
cd ~/sagear_cloud_batch
./setup_vm.sh
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alderaan
JOBS=30 ./run_batch.sh
```

## Retrieve Results

```bash
gcloud compute scp --recurse alderaan-e2-32:~/sagear_cloud_batch/alderaan_project/Results ./Results --zone $ZONE
gcloud compute scp --recurse alderaan-e2-32:~/sagear_cloud_batch/logs ./logs --zone $ZONE
```

Then run locally:

```powershell
& "C:\\Users\\shres\\anaconda3\\python.exe" sagear_reproduction\\extract_eccentricity_posteriors.py
& "C:\\Users\\shres\\anaconda3\\python.exe" sagear_reproduction\\hierarchical_rayleigh.py
```
""",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
