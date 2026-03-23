#!/usr/bin/env bash
# Build LAMMPS ice data files for cubic supercells n×n×n, n = 1..5.
# Molecules N = 8 n³ → data.ice_uma_${N} in tests/uma_ice_rpmd/lammps/
#
# Run once from anywhere:
#   bash tests/uma_ice_rpmd/slurm_campaign_100ps/prepare_lammps_data.sh
set -euo pipefail

UMA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$UMA_DIR"

for n in 1 2 3 4 5; do
  molecules=$((8 * n * n * n))
  out="lammps/data.ice_uma_${molecules}"
  if [[ -f "$out" ]]; then
    echo "Exists (skip): $out"
    continue
  fi
  echo "Building $out (ice Ih ${n}×${n}×${n} = ${molecules} molecules)..."
  python lammps/build_lammps_ice_data.py --nx "$n" --ny "$n" --nz "$n" -o "$out"
done

echo "Done. Files under lammps/data.ice_uma_*"
