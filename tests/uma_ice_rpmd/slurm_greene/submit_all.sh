#!/bin/bash
# Create output dirs and submit the RPMD benchmark array job.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
UMA_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$SCRIPT_DIR"
mkdir -p logs slurm_out

# Ensure LAMMPS data exists for 128 molecules (required by all models)
if [[ ! -f "$UMA_DIR/lammps/data.ice_uma_128" ]]; then
  echo "Building lammps/data.ice_uma_128..."
  python "$UMA_DIR/lammps/build_lammps_ice_data.py" --nx 2 --ny 2 --nz 4 -o "$UMA_DIR/lammps/data.ice_uma_128"
fi

sbatch --job-name=rpmd_bench --array=0-3899 --time=48:00:00 submit_slurm.sh
echo "Submitted. Check logs/ for stdout/stderr."
