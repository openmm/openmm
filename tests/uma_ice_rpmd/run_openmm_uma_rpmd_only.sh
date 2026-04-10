#!/usr/bin/env bash
# OpenMM UMA RPMD only (no TIP4P, no i-PI/LAMMPS).
# Same defaults as step 2 of run_rpmd_32x32.sh.
#
# Usage (from tests/uma_ice_rpmd/):
#   bash run_openmm_uma_rpmd_only.sh
#
# Optional env (same as full pipeline):
#   RPMD_STEPS=10000 RPMD_DT_FS=0.1 OPENMMML_UMA_RPMD_CHUNK=4 bash run_openmm_uma_rpmd_only.sh
#   RPMD_REPORT_EVERY_STEPS=0  # ps-based console interval (faster for long UMA+many-bead runs)
#
# After OpenMM rebuild, point plugins at your build tree, e.g.:
#   export OPENMM_PLUGIN_DIR=/media/extradrive/Trajectories/openmm/build/lib/plugins

set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

CHUNK="${OPENMMML_UMA_RPMD_CHUNK:-4}"
STEPS="${RPMD_STEPS:-10000}"
DT_FS="${RPMD_DT_FS:-0.1}"
EQUIL_PS="${RPMD_EQUIL_PS:-2}"
REPORT_EVERY_STEPS="${RPMD_REPORT_EVERY_STEPS:-5}"

export OPENMMML_UMA_RPMD_CHUNK="$CHUNK"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

if [[ "$EQUIL_PS" == "0" || "$EQUIL_PS" == "0.0" ]]; then
  echo "WARNING: RPMD_EQUIL_PS=0 — no OpenMM NVT equil before production (see docs/IPI_OPENMM_THERMOSTAT_PARITY.md)." >&2
fi

mkdir -p pipeline_out

PROD_PS=$(python3 -c "import sys; print(float(sys.argv[1]) * float(sys.argv[2]) / 1000.0)" "${STEPS}" "${DT_FS}")
echo "=== OpenMM UMA RPMD only: ${STEPS} steps @ ${DT_FS} fs -> ${PROD_PS} ps production | equil=${EQUIL_PS} ps | report_every_steps=${REPORT_EVERY_STEPS} ==="
echo "=== Output: pipeline_out/ice_order_openmm_rpmd.csv ==="

python run_openmm_rpmd_reference.py \
  --nx 2 --ny 2 --nz 2 --beads 32 --dt "$DT_FS" --steps "$STEPS" \
  --equil "$EQUIL_PS" \
  --report-every-steps "$REPORT_EVERY_STEPS" \
  --platform cuda \
  --rpmd-thermostat pile-g

echo "=== Done. ==="
