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

if [ ! -d "$HOME/alderaan" ]; then
  git clone https://github.com/gjgilbert/alderaan "$HOME/alderaan"
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
