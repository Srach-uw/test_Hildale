#!/usr/bin/env bash
# Execute only after the fixture hashes are inspected and recorded.
set -euo pipefail

: "${ALDERAAN_REPO:?Set ALDERAAN_REPO to the pinned ALDERAAN checkout}"
: "${PROJECT_DIR:?Set PROJECT_DIR to the copied K00367 fixture project}"

RUN_ID="${RUN_ID:-diagnostic_k00367_exact_seeded}"
TARGET="${TARGET:-K00367}"
RHO_STAR_SOLAR="${RHO_STAR_SOLAR:-1.0198685421655524}"
SEED="${SEED:-20260925}"
MAXCALL="${MAXCALL:-20000}"
OUTDIR="${OUTDIR:-$PWD/direct_chord_output}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUTDIR"
export PYTHONPATH="$ALDERAAN_REPO${PYTHONPATH:+:$PYTHONPATH}"

python "$SCRIPT_DIR/../../scripts/direct_chord_eccentricity_diagnostic.py" \
  --project-dir "$PROJECT_DIR" \
  --run-id "$RUN_ID" \
  --target "$TARGET" \
  --rho-star-solar "$RHO_STAR_SOLAR" \
  --manifest "$OUTDIR/k00367_direct_chord_manifest.json" \
  --execute \
  --seed "$SEED" \
  --maxcall "$MAXCALL" \
  --output "$OUTDIR/k00367_direct_chord_samples.npz"

sha256sum "$OUTDIR/k00367_direct_chord_manifest.json" "$OUTDIR/k00367_direct_chord_samples.npz"
