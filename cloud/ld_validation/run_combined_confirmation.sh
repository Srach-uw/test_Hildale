#!/usr/bin/env bash
set -euo pipefail

# Run only the post-matrix confirmation arm. Existing completed arms are not
# modified or rerun by this command.
bash run_paper_priors_reference_lcsc.sh
