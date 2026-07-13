#!/usr/bin/env bash
set -euo pipefail
for arm in original_lc reference_lc original_lcsc reference_lcsc original_lc_repeat paper_priors_original_lc; do
  echo "[$(date -Is)] starting $arm"
  bash "run_${arm}.sh"
done
echo "[$(date -Is)] validation matrix complete"
