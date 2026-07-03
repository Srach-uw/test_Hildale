#!/usr/bin/env bash
set -euo pipefail

# Self-heal CRLF in the sibling scripts. A stale/Windows-built bundle zip can
# carry carriage returns that break `set -euo pipefail` ("invalid option name")
# the moment run_batch.sh / run_one_target.sh are executed. Strip them here so
# the run doesn't depend on the zip's line endings being perfect. (This file
# must itself be LF to reach this line; unzip on the VM preserves that, but if
# in doubt run `sed -i 's/\r$//' *.sh` manually right after unzipping.)
sed -i 's/\r$//' "$(dirname "$0")"/*.sh 2>/dev/null || true

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

# Upstream bug (verified 2026-07-02, all commits on main/develop, all tags):
# bin/detrend_and_estimate_ttvs.py unconditionally imports alderaan.validate,
# a module that has never existed in the repo. The two functions it provides
# are only used inside `if MISSION == "Kepler-Validation":`, which our runs
# never enter (we always pass --mission Kepler). Make the import lazy so it
# doesn't crash script startup for the code path we actually use.
DETREND_SCRIPT="$HOME/alderaan/bin/detrend_and_estimate_ttvs.py"
if grep -q '^from alderaan.validate import remove_known_transits, inject_synthetic_transits$' "$DETREND_SCRIPT"; then
  sed -i 's/^from alderaan.validate import remove_known_transits, inject_synthetic_transits$/try:\n    from alderaan.validate import remove_known_transits, inject_synthetic_transits\nexcept ImportError:\n    remove_known_transits = inject_synthetic_transits = None/' "$DETREND_SCRIPT"
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
echo "python validate_bundle.py"
echo "JOBS=30 bash run_batch.sh"
