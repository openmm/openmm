#!/usr/bin/env bash
# Coarse OpenMM-only thermostat sweep for RPMD dynamics parity debugging.
# Varies --rpmd-centroid-friction and --rpmd-friction; writes separate order CSVs.
#
# Usage (from tests/uma_ice_rpmd):
#   bash sweep_openmm_rpmd_thermostat.sh
#
# Environment overrides:
#   SWEEP_STEPS=3000          # production steps (each run)
#   SWEEP_EQUIL_PS=0.5        # OpenMM NVT equilibration before production (ps)
#   SWEEP_BEADS=8
#   SWEEP_DT_FS=0.1
#   OPENMMML_UMA_RPMD_CHUNK=4
#
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

STEPS="${SWEEP_STEPS:-3000}"
EQUIL="${SWEEP_EQUIL_PS:-0.5}"
BEADS="${SWEEP_BEADS:-8}"
DT_FS="${SWEEP_DT_FS:-0.1}"
CHUNK="${OPENMMML_UMA_RPMD_CHUNK:-4}"
OUTDIR="${SWEEP_OUTDIR:-$ROOT/pipeline_out/thermostat_sweep}"

export OPENMMML_UMA_RPMD_CHUNK="$CHUNK"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

mkdir -p "$OUTDIR"
SUMMARY="$OUTDIR/summary.txt"
: > "$SUMMARY"
echo "sweep_openmm_rpmd_thermostat: STEPS=$STEPS EQUIL_PS=$EQUIL BEADS=$BEADS DT_FS=$DT_FS CHUNK=$CHUNK" | tee -a "$SUMMARY"

# Centroid friction (1/ps) × internal PILE friction (1/ps)
for cf in 0.25 0.5 1.0; do
  for inf in 0.5 1.0 2.0; do
    tag="cf${cf}_if${inf}"
    csv="$OUTDIR/ice_order_openmm_${tag}.csv"
    log="$OUTDIR/log_${tag}.txt"
    echo "--- Running OpenMM cf=$cf if=$inf -> $csv (log $log) ---" | tee -a "$SUMMARY"
    if python run_openmm_rpmd_reference.py \
      --nx 2 --ny 2 --nz 2 \
      --beads "$BEADS" \
      --dt "$DT_FS" \
      --steps "$STEPS" \
      --equil "$EQUIL" \
      --platform cuda \
      --rpmd-thermostat pile-g \
      --rpmd-centroid-friction "$cf" \
      --rpmd-friction "$inf" \
      --order-csv "$csv" \
      --order-every 10 \
      --report-every-steps 50 \
      >"$log" 2>&1; then
      echo "  OK $tag" | tee -a "$SUMMARY"
    else
      echo "  FAIL $tag (log $log)" | tee -a "$SUMMARY"
    fi
  done
done

echo "Done. Inspect early T_K / Q6 in CSVs; compare to i-PI tau in ipi/input.xml (default 1000 fs)." | tee -a "$SUMMARY"
echo "Plot example: python plot_rpmd_comparison.py --openmm $OUTDIR/ice_order_openmm_cf0.5_if1.0.csv --ipi pipeline_out/ice_order_ipi_rpmd.csv -o $OUTDIR/compare_cf0.5_if1.0.png"
