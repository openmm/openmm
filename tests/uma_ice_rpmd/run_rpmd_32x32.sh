#!/usr/bin/env bash
# RPMD pipeline: OpenMM UMA + i-PI/LAMMPS UMA only (PILE-G / pile_g). TIP4P is skipped.
# 64 H2O (ice Ih supercell 2×2×2), 32 beads. Run from tests/uma_ice_rpmd/.
#
# Force client (step 2): LAMMPS by default (fix ipi + fix external UMA).
#   export RPMD_IPI_CLIENT=python    # ASE Python driver instead
#
# Large batched UMA RPMD can OOM on ~12 GB GPUs; chunk inference:
#   export OPENMMML_UMA_RPMD_CHUNK=4   # or 8
#   export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
#   bash run_rpmd_32x32.sh
# Console progress every N OpenMM steps (default ps-based auto interval). For long runs use:
#   RPMD_REPORT_EVERY_STEPS=0
# Exact i-PI pile_g centroid tau (fs) used to derive the OpenMM centroid thermostat coupling:
#   RPMD_PILE_TAU_FS=1000
#
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

CHUNK="${OPENMMML_UMA_RPMD_CHUNK:-4}"
IPI_CLIENT="${RPMD_IPI_CLIENT:-lammps}"
# 1 ps @ 0.25 fs = 4000 steps (path-integral order CSV; PILE-G / pile_g)
STEPS="${RPMD_STEPS:-4000}"
DT_FS="${RPMD_DT_FS:-0.25}"
# Default to no separate equilibration so the two legs run the same short parity protocol.
EQUIL_PS="${RPMD_EQUIL_PS:-0}"
# Console progress: every N integrator steps (many beads → expensive). 0 = ps-based auto interval.
REPORT_EVERY_STEPS="${RPMD_REPORT_EVERY_STEPS:-0}"
PILE_TAU_FS="${RPMD_PILE_TAU_FS:-1000}"

export OPENMMML_UMA_RPMD_CHUNK="$CHUNK"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

if [[ "$EQUIL_PS" == "0" || "$EQUIL_PS" == "0.0" ]]; then
  echo "WARNING: RPMD_EQUIL_PS=0 — OpenMM has no NVT equilibration before production order rows." >&2
  echo "         Early centroid T_K and structural CVs often disagree with i-PI+LAMMPS." >&2
  echo "         Use default (2 ps) or docs/RPMD_DYNAMICS_PARITY_RUN_METADATA.md for logging." >&2
fi

echo "=== 0) Stop prior RPMD pipeline processes on this machine (best-effort) ==="
for pat in \
  "run_openmm_rpmd_reference.py" \
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
OPENMM_CENTROID_FRICTION=$(python3 -c "import sys; print(1000.0/float(sys.argv[1]))" "${PILE_TAU_FS}")
echo "=== Settings: steps=${STEPS} dt=${DT_FS} fs -> ${PROD_PS} ps production | equil=${EQUIL_PS} ps | report_every_steps=${REPORT_EVERY_STEPS} | pile_tau=${PILE_TAU_FS} fs ==="

echo "=== 1) OpenMM UMA RPMD (PILE-G via test_uma_ice_rpmd defaults) ==="
python run_openmm_rpmd_reference.py \
  --nx 2 --ny 2 --nz 2 --beads 32 --dt "$DT_FS" --steps "$STEPS" \
  --equil "$EQUIL_PS" \
  --report-every-steps "$REPORT_EVERY_STEPS" \
  --platform cuda \
  --rpmd-thermostat pile-g \
  --rpmd-centroid-friction "$OPENMM_CENTROID_FRICTION"
python rpmd_output_validation.py \
  --csv "pipeline_out/ice_order_openmm_rpmd.csv" \
  --expected-final-ps "$PROD_PS" \
  --label "OpenMM order CSV"

echo "=== 2) i-PI + UMA RPMD + path-integral order CSV (same estimator as OpenMM UMA) ==="
python run_ipi_lammps_uma_rpmd.py \
  --client "$IPI_CLIENT" \
  --molecules 64 --beads 32 --dt-fs "$DT_FS" --steps "$STEPS" \
  --device cuda \
  --ipi-thermostat pile_g \
  --ipi-tau-fs "$PILE_TAU_FS" \
  --order-csv pipeline_out/ice_order_ipi_rpmd.csv
python rpmd_output_validation.py \
  --csv "pipeline_out/ice_order_ipi_rpmd.csv" \
  --expected-final-ps "$PROD_PS" \
  --label "i-PI order CSV"

echo "=== 3) Comparison plot (OpenMM UMA vs i-PI+LAMMPS UMA; no TIP4P) ==="
python plot_rpmd_comparison.py \
  --openmm pipeline_out/ice_order_openmm_rpmd.csv \
  --ipi pipeline_out/ice_order_ipi_rpmd.csv \
  -o pipeline_out/rpmd_comparison_32x32.png \
  --title "RPMD ice order: OpenMM UMA vs i-PI+LAMMPS UMA (64 mol, 32 beads, 243 K, PILE-G)"

echo "Done."
