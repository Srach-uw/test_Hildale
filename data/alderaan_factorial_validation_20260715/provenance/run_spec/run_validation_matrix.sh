#!/usr/bin/env bash
set -euo pipefail

OVERALL_RC=0
for arm in original_lc reference_lc original_lcsc reference_lcsc original_lc_repeat paper_priors_original_lc; do
  echo "[$(date -Is)] starting $arm"
  if bash "run_${arm}.sh"; then
    echo "[$(date -Is)] $arm complete"
  else
    echo "[$(date -Is)] $arm incomplete; continuing to the next validation arm" >&2
    OVERALL_RC=1
  fi
done
if [ "$OVERALL_RC" -ne 0 ]; then
  echo "[$(date -Is)] validation matrix finished with incomplete arms" >&2
  exit 1
fi
echo "[$(date -Is)] validation matrix complete"
