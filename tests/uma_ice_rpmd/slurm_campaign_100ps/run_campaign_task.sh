#!/usr/bin/env bash
# One SLURM array task: RPMD 100 ps (default), cubic supercell 1³–5³, beads 4–64, 100 replicas.
#
# Decode SLURM_ARRAY_TASK_ID (0–2499):
#   replica   = ID % 100
#   rest      = ID / 100   (0..24)
#   bead_idx  = rest % 5    → 4,8,16,32,64
#   size_idx  = rest / 5    → edge n = 1..5
#
# Model: CAMPAIGN_MODEL=tip4p | openmm_uma | lammps_uma (required)
#
# Optional overrides for dry runs:
#   RPMD_STEPS=1000     → pass --steps to drivers (overrides production length)
#   PROD_PS=1.0         → production time in ps if RPMD_STEPS unset (default 100)
#   RPMD_DT_FS=0.1      → timestep fs (default 0.1)
#   ORDER_EVERY=100     → order CSV row stride (TIP4P / OpenMM UMA order output)
#   CAMPAIGN_SEED_BASE  → default 284759; seed = BASE + ID
set -euo pipefail

UMA_DIR="$(cd "$(dirname "$0")/.." && pwd)"
CAMPAIGN_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$UMA_DIR"

MODEL="${CAMPAIGN_MODEL:-}"
if [[ -z "$MODEL" ]]; then
  echo "ERROR: Set CAMPAIGN_MODEL to tip4p, openmm_uma, or lammps_uma" >&2
  exit 1
fi

ID="${SLURM_ARRAY_TASK_ID:-0}"
BEADS_LIST=(4 8 16 32 64)

replica=$((ID % 100))
rest=$((ID / 100))
bead_idx=$((rest % 5))
size_idx=$((rest / 5))
n=$((size_idx + 1))
beads=${BEADS_LIST[$bead_idx]}
nmol=$((8 * n * n * n))

BASE_SEED="${CAMPAIGN_SEED_BASE:-284759}"
SEED=$((BASE_SEED + ID))

DT_FS="${RPMD_DT_FS:-0.1}"
ORDER_EVERY="${ORDER_EVERY:-100}"
PROD_PS="${PROD_PS:-100}"

DATA_FILE="lammps/data.ice_uma_${nmol}"
OUTDIR="${CAMPAIGN_DIR}/slurm_out/${MODEL}/n${n}_b${beads}/replica_$(printf '%03d' "$replica")"
mkdir -p "$OUTDIR"

# Record metadata for downstream analysis
{
  echo "SLURM_ARRAY_TASK_ID=${ID}"
  echo "CAMPAIGN_MODEL=${MODEL}"
  echo "supercell_edge=${n}"
  echo "molecules=${nmol}"
  echo "beads=${beads}"
  echo "replica=${replica}"
  echo "seed=${SEED}"
  echo "prod_ps_default=${PROD_PS}"
  echo "dt_fs=${DT_FS}"
} > "$OUTDIR/task_meta.txt"

if [[ ! -f "$DATA_FILE" ]]; then
  echo "ERROR: Missing ${UMA_DIR}/${DATA_FILE}. Run: bash slurm_campaign_100ps/prepare_lammps_data.sh" >&2
  exit 1
fi

echo "Task ${ID}: model=${MODEL} n=${n}³ nmol=${nmol} beads=${beads} replica=${replica} seed=${SEED}"
echo "OUTDIR=${OUTDIR}"

# --- TIP4P OpenMM RPMD ---
if [[ "$MODEL" == "tip4p" ]]; then
  extra=()
  if [[ -n "${RPMD_STEPS:-}" ]]; then
    extra+=(--steps "${RPMD_STEPS}")
  else
    extra+=(--prod "${PROD_PS}")
  fi
  python run_openmm_tip4p_rpmd.py \
    --data "${DATA_FILE}" \
    --molecules "${nmol}" \
    --beads "${beads}" \
    --dt "${DT_FS}" \
    "${extra[@]}" \
    --rpmd-thermostat pile-g \
    --order-csv "${OUTDIR}/order.csv" \
    --order-every "${ORDER_EVERY}" \
    --seed "${SEED}" \
    --platform cuda
  echo "Done: ${OUTDIR}/order.csv"
  exit 0
fi

# --- OpenMM UMA RPMD (delegates to test_uma_ice_rpmd) ---
if [[ "$MODEL" == "openmm_uma" ]]; then
  extra=()
  if [[ -n "${RPMD_STEPS:-}" ]]; then
    extra+=(--steps "${RPMD_STEPS}")
  else
    extra+=(--prod "${PROD_PS}")
  fi
  python run_openmm_rpmd_reference.py \
    --nx "$n" --ny "$n" --nz "$n" \
    --beads "${beads}" \
    --dt "${DT_FS}" \
    "${extra[@]}" \
    --order-csv "${OUTDIR}/order.csv" \
    --seed "${SEED}" \
    --platform cuda
  echo "Done: ${OUTDIR}/order.csv"
  exit 0
fi

# --- i-PI + LAMMPS UMA (isolated workdir; unique port) ---
if [[ "$MODEL" == "lammps_uma" ]]; then
  PORT=$((64511 + (ID % 10000)))
  WORKDIR="${TMPDIR:-/tmp}/rpmd_campaign_${SLURM_JOB_ID:-local}_${ID}"
  rm -rf "$WORKDIR"
  mkdir -p "$WORKDIR"
  # Avoid copying accumulated slurm_out (large); prefer rsync when available.
  if command -v rsync >/dev/null 2>&1; then
    rsync -a \
      --exclude 'slurm_campaign_100ps/slurm_out/' \
      --exclude '.git/' \
      --exclude '__pycache__/' \
      "${UMA_DIR}/" "${WORKDIR}/"
  else
    cp -a "${UMA_DIR}/." "${WORKDIR}/"
    rm -rf "${WORKDIR}/slurm_campaign_100ps/slurm_out"
  fi
  cd "$WORKDIR"

  extra_s=()
  if [[ -n "${RPMD_STEPS:-}" ]]; then
    extra_s+=(--steps "${RPMD_STEPS}")
  else
    extra_s+=(--prod "${PROD_PS}")
  fi

  python run_ipi_lammps_uma_rpmd.py \
    --molecules "${nmol}" \
    --beads "${beads}" \
    --dt-fs "${DT_FS}" \
    "${extra_s[@]}" \
    --seed "${SEED}" \
    --port "${PORT}" \
    --data "${DATA_FILE}" \
    --device cuda \
    --ipi-thermostat pile_g \
    --order-csv "${OUTDIR}/order.csv"

  rm -rf "$WORKDIR"
  cd "$UMA_DIR"
  echo "Done: ${OUTDIR}/order.csv"
  exit 0
fi

echo "ERROR: Unknown CAMPAIGN_MODEL=${MODEL}" >&2
exit 1
