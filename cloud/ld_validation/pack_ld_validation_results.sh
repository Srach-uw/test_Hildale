#!/usr/bin/env bash
set -euo pipefail
OUT="alderaan_factorial_validation_$(date +%Y%m%d_%H%M%S).tar.gz"
tar -czf "$OUT" projects logs provenance *.csv README_LD_VALIDATION.md run_*.sh patch_alderaan_*.py
echo "$OUT"
