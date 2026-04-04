#!/bin/bash
# Dispatch a single RPMD job by array index.
# Mapping: id % 3 = model (tip4p, openmm_uma, lammps_uma)
#          (id // 3) % 13 = bead index, (id // 3) // 13 = replica
set -euo pipefail

ID="${SLURM_ARRAY_TASK_ID:-0}"
BEADS=(8 10 12 14 16 18 20 22 24 26 28 30 32)

model_idx=$((ID % 3))
id2=$((ID / 3))
bead_idx=$((id2 % 13))
replica=$((id2 / 13))
bead=${BEADS[$bead_idx]}

case $model_idx in
  0) model="tip4p" ;;
  1) model="openmm_uma" ;;
  2) model="lammps_uma" ;;
  *) echo "Invalid model_idx $model_idx"; exit 1 ;;
esac

SEED=$((284759 + replica))
OUTDIR="${SLURM_SUBMIT_DIR:-.}/slurm_out/${model}/replica_$(printf '%03d' "$replica")/bead_$(printf '%02d' "$bead")"
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$OUTDIR"

# For LAMMPS UMA: per-task port and work dir to avoid i-PI collisions
PORT=$((64511 + (ID % 10000)))
WORKDIR="${TMPDIR:-/tmp}/rpmd_${SLURM_JOB_ID:-$$}_${ID}"

echo "Task $ID: model=$model bead=$bead replica=$replica seed=$SEED outdir=$OUTDIR"

cd "$SCRIPT_DIR"

case $model in
  tip4p)
    python run_openmm_tip4p_rpmd.py \
      --data lammps/data.ice_uma_128 \
      --molecules 128 \
      --beads "$bead" \
      --dt 0.1 \
      --prod 100 \
      --order-csv "$OUTDIR/order.csv" \
      --output "$OUTDIR" \
      --seed "$SEED"
    ;;
  openmm_uma)
    python run_openmm_rpmd_reference.py \
      --molecules 128 \
      --beads "$bead" \
      --dt 0.1 \
      --prod 100 \
      --order-csv "$OUTDIR/order.csv" \
      --seed "$SEED"
    ;;
  lammps_uma)
    # Copy project to isolated work dir (i-PI writes to ipi/; avoid cross-task collisions)
    mkdir -p "$WORKDIR"
    cp -r "$SCRIPT_DIR"/* "$WORKDIR/"
    cd "$WORKDIR"
    python run_ipi_lammps_uma_rpmd.py \
      --molecules 128 \
      --beads "$bead" \
      --prod 100 \
      --dt-fs 0.1 \
      --seed "$SEED" \
      --port "$PORT" \
      --data lammps/data.ice_uma_128 \
      --order-csv "$OUTDIR/order.csv"
    rm -rf "$WORKDIR"
    ;;
  *)
    echo "Unknown model $model"; exit 1 ;;
esac

echo "Done: $OUTDIR/order.csv"
