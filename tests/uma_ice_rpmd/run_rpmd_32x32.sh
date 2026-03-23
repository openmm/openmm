#!/usr/bin/env bash
# Full RPMD pipeline: TIP4P/2005f + OpenMM UMA + i-PI/LAMMPS UMA (all PILE-G / pile_g).
# 32 H2O, 32 beads. Run from tests/uma_ice_rpmd/.
#
# Large batched UMA RPMD can OOM on ~12 GB GPUs; chunk inference:
#   export OPENMMML_UMA_RPMD_CHUNK=4   # or 8
#   export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
#   bash run_rpmd_32x32.sh
#
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

CHUNK="${OPENMMML_UMA_RPMD_CHUNK:-4}"
# 1 ps @ 0.1 fs = 10000 steps (path-integral order CSV; PILE-G / pile_g)
STEPS="${RPMD_STEPS:-10000}"
DT_FS="${RPMD_DT_FS:-0.1}"

export OPENMMML_UMA_RPMD_CHUNK="$CHUNK"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

echo "=== 0) Stop prior RPMD pipeline processes on this machine (best-effort) ==="
for pat in \
  "run_openmm_rpmd_reference.py" \
  "run_openmm_tip4p_rpmd.py" \
  "run_ipi_lammps_uma_rpmd.py" \
  "test_uma_ice_rpmd.py" \
  ; do
  pkill -TERM -f "$pat" 2>/dev/null || true
done
pkill -TERM -f "i-pi.*input.xml" 2>/dev/null || true
pkill -TERM -f "i-pi input.xml" 2>/dev/null || true
sleep 2

mkdir -p pipeline_out

# Avoid nested "..." inside $(...) — a `)` in print(...) can break command substitution.
PROD_PS=$(python3 -c "import sys; print(float(sys.argv[1]) * float(sys.argv[2]) / 1000.0)" "${STEPS}" "${DT_FS}")
echo "=== Settings: steps=${STEPS} dt=${DT_FS} fs -> ${PROD_PS} ps production ==="

echo "=== 1) OpenMM TIP4P/2005f RPMD (PILE-G), order CSV ==="
python run_openmm_tip4p_rpmd.py \
  --molecules 32 --beads 32 --dt "$DT_FS" --steps "$STEPS" \
  --rpmd-thermostat pile-g \
  --order-csv pipeline_out/ice_order_tip4p_rpmd.csv \
  --order-every 10 \
  --platform cuda

echo "=== 2) OpenMM UMA RPMD (PILE-G via test_uma_ice_rpmd defaults) ==="
python run_openmm_rpmd_reference.py \
  --molecules 32 --beads 32 --dt "$DT_FS" --steps "$STEPS" \
  --platform cuda \
  --rpmd-thermostat pile-g

echo "=== 3) i-PI + LAMMPS UMA (thermostat pile_g, match dt-fs) ==="
python run_ipi_lammps_uma_rpmd.py \
  --molecules 32 --beads 32 --dt-fs "$DT_FS" --steps "$STEPS" \
  --device cuda \
  --ipi-thermostat pile_g

echo "=== 4) i-PI trajectory -> order CSV (path-integral order; i-PI 3: traj_00) ==="
python ipi_order_from_traj.py \
  --traj ipi/ice__i-pi.traj_00.xyz --beads 32 --dt-fs "$DT_FS" \
  -o pipeline_out/ice_order_ipi_rpmd.csv

echo "=== 5) Comparison plot (TIP4P + OpenMM UMA + i-PI) ==="
python plot_rpmd_comparison.py \
  --openmm pipeline_out/ice_order_openmm_rpmd.csv \
  --ipi pipeline_out/ice_order_ipi_rpmd.csv \
  --tip4p pipeline_out/ice_order_tip4p_rpmd.csv \
  -o pipeline_out/rpmd_comparison_32x32.png \
  --title "RPMD ice order: TIP4P vs OpenMM UMA vs i-PI+LAMMPS UMA (32 mol, 32 beads, 243 K, PILE-G)"

echo "Done."
