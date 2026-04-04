#!/usr/bin/env bash
# OpenMM UMA RPMD → i-PI + LAMMPS UMA → comparison plot, with timestamped outputs
# so pipeline_out/ice_order_openmm_rpmd.csv and friends are not overwritten.
#
# Prerequisites: conda env with OpenMM, openmmml, fairchem, i-pi, LAMMPS Python API.
#   conda activate fairchem   # or your env
#   export OPENMM_PLUGIN_DIR="${CONDA_PREFIX}/lib/plugins"
#
# Run from anywhere:
#   bash tests/uma_ice_rpmd/run_rpmd_openmm_ipi_stamped.sh
#
# Defaults: minimization + 0.1 ps production (RPMD_STEPS=400 at RPMD_DT_FS=0.25)
# with no extra equilibration. Override RPMD_STEPS / RPMD_DT_FS for longer runs.
#
# i-PI leg: ``run_ipi_lammps_uma_rpmd.py`` now runs LAMMPS+UMA minimize before i-PI by default
# (OpenMM-matched starting geometry). Use ``--no-ipi-minimize`` on that script for legacy raw-data init.
#
# Same knobs as run_rpmd_32x32.sh (env overrides):
#   RPMD_STEPS RPMD_DT_FS RPMD_EQUIL_PS RPMD_REPORT_EVERY_STEPS RPMD_PILE_TAU_FS
#   OPENMMML_UMA_RPMD_CHUNK PYTORCH_CUDA_ALLOC_CONF RPMD_IPI_CLIENT
#
# i-PI still removes stale ipi/ice_*i-pi* before each run; copy ipi/ aside first
# if you need to keep previous bead trajectories.
#
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

CHUNK="${OPENMMML_UMA_RPMD_CHUNK:-4}"
IPI_CLIENT="${RPMD_IPI_CLIENT:-lammps}"
# 400 steps × 0.25 fs = 0.1 ps production (override RPMD_STEPS / RPMD_DT_FS as needed)
STEPS="${RPMD_STEPS:-400}"
DT_FS="${RPMD_DT_FS:-0.25}"
EQUIL_PS="${RPMD_EQUIL_PS:-0}"
REPORT_EVERY_STEPS="${RPMD_REPORT_EVERY_STEPS:-0}"
PILE_TAU_FS="${RPMD_PILE_TAU_FS:-1000}"
OUT_DIR="${RPMD_PIPELINE_OUT:-pipeline_out}"

export OPENMMML_UMA_RPMD_CHUNK="$CHUNK"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

STAMP="${RPMD_STAMP:-$(date +%Y%m%d_%H%M%S)}"
OPENMM_CSV="${OUT_DIR}/ice_order_openmm_rpmd_${STAMP}.csv"
IPI_CSV="${OUT_DIR}/ice_order_ipi_rpmd_${STAMP}.csv"
PLOT_PNG="${OUT_DIR}/rpmd_comparison_${STAMP}.png"
IPI_MD="${OUT_DIR}/ice__i-pi_${STAMP}.md"

if [[ "$EQUIL_PS" == "0" || "$EQUIL_PS" == "0.0" ]]; then
  echo "INFO: OpenMM leg will run minimization + velocity initialization + production (no separate equilibration)." >&2
fi

echo "=== Stamped run: STAMP=${STAMP} | OUT_DIR=${OUT_DIR} ==="
echo "    OpenMM order CSV: ${OPENMM_CSV}"
echo "    i-PI order CSV:   ${IPI_CSV}"
echo "    Plot:             ${PLOT_PNG}"

mkdir -p "$OUT_DIR"

echo "=== Stop stale i-PI (best-effort) ==="
pkill -TERM -f "i-pi.*input.xml" 2>/dev/null || true
pkill -TERM -f "i-pi input.xml" 2>/dev/null || true
sleep 2

PROD_PS=$(python3 -c "import sys; print(float(sys.argv[1]) * float(sys.argv[2]) / 1000.0)" "${STEPS}" "${DT_FS}")
OPENMM_CENTROID_FRICTION=$(python3 -c "import sys; print(1000.0/float(sys.argv[1]))" "${PILE_TAU_FS}")
echo "=== Settings: steps=${STEPS} dt=${DT_FS} fs -> ${PROD_PS} ps production | equil=${EQUIL_PS} ps | chunk=${CHUNK} | client=${IPI_CLIENT} | pile_tau=${PILE_TAU_FS} fs ==="

echo "=== 1) OpenMM UMA RPMD ==="
python run_openmm_rpmd_reference.py \
  --nx 2 --ny 2 --nz 2 --beads 32 --dt "$DT_FS" --steps "$STEPS" \
  --equil "$EQUIL_PS" \
  --report-every-steps "$REPORT_EVERY_STEPS" \
  --platform cuda \
  --rpmd-thermostat pile-g \
  --rpmd-centroid-friction "$OPENMM_CENTROID_FRICTION" \
  --order-csv "$OPENMM_CSV"
python rpmd_output_validation.py \
  --csv "$OPENMM_CSV" \
  --expected-final-ps "$PROD_PS" \
  --label "OpenMM stamped order CSV"

echo "=== 2) i-PI + UMA (${IPI_CLIENT}) ==="
python run_ipi_lammps_uma_rpmd.py \
  --client "$IPI_CLIENT" \
  --molecules 64 --beads 32 --dt-fs "$DT_FS" --steps "$STEPS" \
  --device cuda \
  --ipi-thermostat pile_g \
  --ipi-tau-fs "$PILE_TAU_FS" \
  --order-csv "$IPI_CSV"
python rpmd_output_validation.py \
  --csv "$IPI_CSV" \
  --expected-final-ps "$PROD_PS" \
  --label "i-PI stamped order CSV"
if [[ -f "ipi/ice__i-pi.md" ]]; then
  cp "ipi/ice__i-pi.md" "$IPI_MD"
else
  echo "ERROR: missing i-PI thermo file ipi/ice__i-pi.md; stamped CSV lacks kinetic temperature merge." >&2
  exit 1
fi
python ipi_order_from_traj.py \
  --traj "ipi/ice__i-pi.traj_00.xyz" \
  --beads 32 \
  --dt-fs "$DT_FS" \
  --thermo "$IPI_MD" \
  -o "$IPI_CSV"
python rpmd_output_validation.py \
  --csv "$IPI_CSV" \
  --expected-final-ps "$PROD_PS" \
  --label "i-PI stamped order CSV"

echo "=== 3) Comparison plot ==="
python plot_rpmd_comparison.py \
  --openmm "$OPENMM_CSV" \
  --ipi "$IPI_CSV" \
  -o "$PLOT_PNG" \
  --title "RPMD ice order: OpenMM UMA vs i-PI+LAMMPS UMA (64 mol, 32 beads, 243 K, PILE-G)"

echo "Done. Outputs:"
echo "  ${OPENMM_CSV}"
echo "  ${IPI_CSV}"
echo "  ${PLOT_PNG}"
