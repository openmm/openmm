#!/usr/bin/env bash
# OpenMM UMA RPMD only (no TIP4P, no i-PI/LAMMPS).
# Same defaults as step 2 of run_rpmd_32x32.sh.
#
# Usage (from tests/uma_ice_rpmd/):
#   bash run_openmm_uma_rpmd_only.sh
#
# Optional env (same as full pipeline):
#   RPMD_STEPS=10000 RPMD_DT_FS=0.1 OPENMMML_UMA_RPMD_CHUNK=4 bash run_openmm_uma_rpmd_only.sh
#
# After OpenMM rebuild, point plugins at your build tree, e.g.:
#   export OPENMM_PLUGIN_DIR=/media/extradrive/Trajectories/openmm/build/lib/plugins

set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

CHUNK="${OPENMMML_UMA_RPMD_CHUNK:-4}"
STEPS="${RPMD_STEPS:-10000}"
DT_FS="${RPMD_DT_FS:-0.1}"

export OPENMMML_UMA_RPMD_CHUNK="$CHUNK"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

mkdir -p pipeline_out

PROD_PS=$(python3 -c "import sys; print(float(sys.argv[1]) * float(sys.argv[2]) / 1000.0)" "${STEPS}" "${DT_FS}")
echo "=== OpenMM UMA RPMD only: ${STEPS} steps @ ${DT_FS} fs -> ${PROD_PS} ps production ==="
echo "=== Output: pipeline_out/ice_order_openmm_rpmd.csv ==="

python run_openmm_rpmd_reference.py \
  --molecules 32 --beads 32 --dt "$DT_FS" --steps "$STEPS" \
  --platform cuda \
  --rpmd-thermostat pile-g

echo "=== Done. ==="
