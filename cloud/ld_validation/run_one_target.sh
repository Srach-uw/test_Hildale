#!/usr/bin/env bash
set -euo pipefail

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

TARGET="${1:?target required}"
KEPID="${2:?kepid required}"
SLOT="${3:-shared}"
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$PWD/alderaan_project}"
RUN_ID="${RUN_ID:-sagear_missing}"
MISSION="${MISSION:-Kepler}"
ALDERAAN_REPO="${ALDERAAN_REPO:-$HOME/alderaan_sagear_pinned}"
CADENCE_MODE="${CADENCE_MODE:-both}"
SEED_OFFSET="${SEED_OFFSET:-0}"
# ALDERAAN has no setup.py/installer; bin/*.py do `import alderaan.io`, which
# only resolves if the repo root is on PYTHONPATH.
export PYTHONPATH="$ALDERAAN_REPO${PYTHONPATH:+:$PYTHONPATH}"
# Give each concurrent GNU-parallel job slot (passed as $3 = {%}) its own
# Theano compiledir. Theano's compile lock is a single global lock per
# compiledir with a polling wait (not a native futex) - at high JOBS with a
# mostly-cold cache this serializes almost all real work. Isolating by slot
# means a slot's lock is never contended by any other slot, while targets
# processed sequentially by the same slot still reuse that slot's cache.
# Falls back to the shared default dir for ad-hoc single-target invocations
# that don't pass a slot number.
if [ "$SLOT" != "shared" ]; then
  export THEANO_FLAGS="base_compiledir=$HOME/.theano_slot_${SLOT}${THEANO_FLAGS:+,$THEANO_FLAGS}"
fi
DATA_DIR="$PROJECT_DIR/Data"
CATALOG_NAME="${CATALOG_NAME:-sagear_missing_catalog.csv}"
RESULTS_FITS="$PROJECT_DIR/Results/$RUN_ID/$TARGET/$TARGET-results.fits"
STATUS_DIR="$PROJECT_DIR/status"

mkdir -p "$DATA_DIR" "$STATUS_DIR" logs "$BUNDLE_DIR/provenance"

case "$CADENCE_MODE" in
  long|both) ;;
  *) echo "CADENCE_MODE must be long or both, got: $CADENCE_MODE" >&2; exit 2 ;;
esac

SEED_HEX="$(printf '%s' "$TARGET" | sha256sum | cut -c1-8)"
export ALDERAAN_SEED="$(( (16#$SEED_HEX + SEED_OFFSET) % 2147483647 ))"

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
if [ "$CADENCE_MODE" = "both" ]; then
  python "$ALDERAAN_REPO/bin/get_kepler_data.py" "$KEPID" -c short -t lightcurve -o "get_${KEPID}_sc.sh" --cmdtype wget
  bash "get_${KEPID}_sc.sh" || true
fi
popd >/dev/null

LC_COUNT="$(find "$DATA_DIR" -maxdepth 1 -type f -name "kplr$(printf '%09d' "$KEPID")-*_llc.fits" | wc -l)"
SC_COUNT="$(find "$DATA_DIR" -maxdepth 1 -type f -name "kplr$(printf '%09d' "$KEPID")-*_slc.fits" | wc -l)"
if [ "$LC_COUNT" -eq 0 ]; then
  echo "No long-cadence FITS downloaded for $TARGET / $KEPID" >&2
  echo "missing_lightcurve" > "$STATUS_DIR/$TARGET.status"
  exit 4
fi

pushd "$ALDERAAN_REPO" >/dev/null
python bin/detrend_and_estimate_ttvs.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR" --data_dir "$DATA_DIR/" --catalog "$CATALOG_NAME"
# --data_dir/--catalog are required=True here too; omitting them makes
# argparse raise SystemExit, which this script's own top-level
# `except SystemExit: warnings.warn(...)` silently swallows instead of
# exiting - execution then continues with MISSION/TARGET/etc never
# assigned, crashing later with a confusing NameError. Found live on the
# smoke-test VM (K00179 failure, 2026-07-02).
python bin/analyze_autocorrelated_noise.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR" --data_dir "$DATA_DIR/" --catalog "$CATALOG_NAME"
python bin/fit_transit_shape_simultaneous_nested.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR"
popd >/dev/null

if [ -s "$RESULTS_FITS" ]; then
  ALDERAAN_COMMIT="$(git -C "$ALDERAAN_REPO" rev-parse HEAD)"
  CATALOG_SHA256="$(sha256sum "$PROJECT_DIR/Catalogs/$CATALOG_NAME" | awk '{print $1}')"
  {
    printf 'target\t%s\n' "$TARGET"
    printf 'kepid\t%s\n' "$KEPID"
    printf 'run_id\t%s\n' "$RUN_ID"
    printf 'alderaan_commit\t%s\n' "$ALDERAAN_COMMIT"
    printf 'catalog_sha256\t%s\n' "$CATALOG_SHA256"
    printf 'rng_seed\t%s\n' "$ALDERAAN_SEED"
    printf 'cadence_mode\t%s\n' "$CADENCE_MODE"
    printf 'long_cadence_files\t%s\n' "$LC_COUNT"
    printf 'short_cadence_files\t%s\n' "$SC_COUNT"
  } > "$BUNDLE_DIR/provenance/${RUN_ID}_${TARGET}.tsv"
  echo "complete" > "$STATUS_DIR/$TARGET.status"
  echo "[$(date -Is)] Finished $TARGET"
else
  echo "missing_results" > "$STATUS_DIR/$TARGET.status"
  echo "[$(date -Is)] Finished without results FITS: $TARGET" >&2
  exit 3
fi
