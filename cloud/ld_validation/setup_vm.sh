#!/usr/bin/env bash
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

ALDERAAN_COMMIT="${ALDERAAN_COMMIT:-7443dff16b7f9092e14a6f0cc1f8948d457c9e0b}"
ALDERAAN_REPO="${ALDERAAN_REPO:-$HOME/alderaan_sagear_pinned}"
if [ ! -d "$ALDERAAN_REPO/.git" ]; then
  git clone https://github.com/gjgilbert/alderaan "$ALDERAAN_REPO"
fi
git -C "$ALDERAAN_REPO" fetch origin "$ALDERAAN_COMMIT"
git -C "$ALDERAAN_REPO" checkout --detach "$ALDERAAN_COMMIT"
ACTUAL_COMMIT="$(git -C "$ALDERAAN_REPO" rev-parse HEAD)"
if [ "$ACTUAL_COMMIT" != "$ALDERAAN_COMMIT" ]; then
  echo "ALDERAAN commit mismatch: expected $ALDERAAN_COMMIT, got $ACTUAL_COMMIT" >&2
  exit 2
fi
python "$PWD/patch_alderaan_repro.py" "$ALDERAAN_REPO"

# Keep the public-code and manuscript-table prior hypotheses in separate clones.
# The latter is sensitivity-only because the paper and public code disagree.
ALDERAAN_PAPER_PRIOR_REPO="${ALDERAAN_PAPER_PRIOR_REPO:-$HOME/alderaan_sagear_paper_priors}"
if [ ! -d "$ALDERAAN_PAPER_PRIOR_REPO/.git" ]; then
  git clone https://github.com/gjgilbert/alderaan "$ALDERAAN_PAPER_PRIOR_REPO"
fi
git -C "$ALDERAAN_PAPER_PRIOR_REPO" fetch origin "$ALDERAAN_COMMIT"
git -C "$ALDERAAN_PAPER_PRIOR_REPO" checkout --detach "$ALDERAAN_COMMIT"
python "$PWD/patch_alderaan_repro.py" "$ALDERAAN_PAPER_PRIOR_REPO"
if [ -f "$PWD/patch_alderaan_paper_priors.py" ]; then
  python "$PWD/patch_alderaan_paper_priors.py" "$ALDERAAN_PAPER_PRIOR_REPO"
fi

if ! conda env list | awk '{print $1}' | grep -qx alderaan; then
  if command -v mamba >/dev/null 2>&1; then
    mamba env create -n alderaan -f "$ALDERAAN_REPO/environment.yml"
  else
    conda env create -n alderaan -f "$ALDERAAN_REPO/environment.yml"
  fi
fi

echo "Ready. Next:"
echo "source ~/miniforge3/etc/profile.d/conda.sh"
echo "conda activate alderaan"
echo "export ALDERAAN_REPO=$ALDERAAN_REPO"
echo "export ALDERAAN_PAPER_PRIOR_REPO=$ALDERAAN_PAPER_PRIOR_REPO"
echo "python validate_bundle.py"
echo "JOBS=30 bash run_batch.sh"
