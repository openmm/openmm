#!/usr/bin/env bash
# Phase 2 minimal dynamics smoke: 5 production RPMD steps (64 H2O, 4 beads).
# Requires OpenMM built with Python ML plugin (PythonForce).
#
#   export OPENMM_PLUGIN_DIR=/path/to/openmm/build/lib/plugins
#   bash run_phase2_minimal.sh
#
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

if [[ -z "${OPENMM_PLUGIN_DIR:-}" ]]; then
  echo "Set OPENMM_PLUGIN_DIR to your OpenMM build plugins directory (contains Python ML plugin)." >&2
  exit 2
fi

export OPENMMML_UMA_RPMD_CHUNK="${OPENMMML_UMA_RPMD_CHUNK:-4}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

echo "=== Phase 2: OpenMM UMA RPMD, 5 steps, order CSV every step ==="
python run_openmm_rpmd_reference.py \
  --nx 2 --ny 2 --nz 2 \
  --beads 4 \
  --dt 0.1 \
  --equil 0 \
  --steps 5 \
  --order-every 1 \
  --report-every-steps 1 \
  --optimize-inference-tf32-only \
  --order-csv "${ROOT}/pipeline_out/phase2_openmm_5step.csv" \
  --platform cuda

echo "=== Optional: static OpenMM vs LAMMPS-style forces (same data file) ==="
python compare_forces_openmm_vs_lammps.py --data "${ROOT}/lammps/data.ice_uma_64" --openmm-platform cuda || true

echo "=== Done ==="
