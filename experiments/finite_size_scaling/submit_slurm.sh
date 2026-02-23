#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32GB
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# Submit as:  sbatch --job-name=fss_N250 --array=0-499 --time=12:00:00 submit_slurm.sh 250
#             sbatch --job-name=fss_N1000 --array=0-124 --time=24:00:00 submit_slurm.sh 1000
#
# For baseline (no cavity):
#             sbatch --job-name=fss_base_N250 --array=0-499 --time=12:00:00 submit_slurm.sh 250 --no-cavity
#
# The first positional argument is N (number of molecules).
# Optional --no-cavity flag runs the lambda=0 baseline.

set -euo pipefail

N=${1:?Usage: sbatch submit_slurm.sh <N_molecules> [--no-cavity]}
shift
EXTRA_FLAGS="$*"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

REPLICA_ID=${SLURM_ARRAY_TASK_ID}

mkdir -p logs

echo "=== SLURM Array Job ==="
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Array ID:   ${SLURM_ARRAY_TASK_ID}"
echo "Node:       $(hostname)"
echo "GPU:        $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "N:          ${N}"
echo "Replica:    ${REPLICA_ID}"
echo "Extra:      ${EXTRA_FLAGS}"
echo ""

# Load modules (adjust for your environment)
module purge 2>/dev/null || true
module load cuda 2>/dev/null || true

# Activate conda/venv if needed (uncomment and adjust)
# source /path/to/conda/etc/profile.d/conda.sh && conda activate openmm
# source /path/to/venv/bin/activate

python3 run_single.py --N "${N}" --replica "${REPLICA_ID}" ${EXTRA_FLAGS}

echo ""
echo "=== Done ==="
