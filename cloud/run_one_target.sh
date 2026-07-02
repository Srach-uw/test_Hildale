#!/usr/bin/env bash
set -euo pipefail

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

TARGET="${1:?target required}"
KEPID="${2:?kepid required}"
PROJECT_DIR="${PROJECT_DIR:-$PWD/alderaan_project}"
RUN_ID="${RUN_ID:-sagear_missing}"
MISSION="${MISSION:-Kepler}"
ALDERAAN_REPO="${ALDERAAN_REPO:-$HOME/alderaan}"
# ALDERAAN has no setup.py/installer; bin/*.py do `import alderaan.io`, which
# only resolves if the repo root is on PYTHONPATH.
export PYTHONPATH="$ALDERAAN_REPO${PYTHONPATH:+:$PYTHONPATH}"
DATA_DIR="$PROJECT_DIR/Data"
CATALOG_NAME="${CATALOG_NAME:-sagear_missing_catalog.csv}"
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
# --cmdtype wget, not curl: archive.stsci.edu 301-redirects http->https, and
# ALDERAAN's generated curl command has no -L, so it silently "succeeds" by
# saving the tiny HTML redirect page instead of the real FITS file. wget
# follows redirects by default and needs no extra flag.
python "$ALDERAAN_REPO/bin/get_kepler_data.py" "$KEPID" -c long -t lightcurve -o "get_${KEPID}_lc.sh" --cmdtype wget
bash "get_${KEPID}_lc.sh" || true
popd >/dev/null

pushd "$ALDERAAN_REPO" >/dev/null
python bin/detrend_and_estimate_ttvs.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR" --data_dir "$DATA_DIR/" --catalog "$CATALOG_NAME"
python bin/analyze_autocorrelated_noise.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR"
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
