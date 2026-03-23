#!/bin/bash
#SBATCH --job-name=rpmd100_tip4p
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=48GB
#SBATCH --time=48:00:00
#SBATCH --array=0-2499
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# TIP4P/2005f OpenMM RPMD — one task per (supercell × bead count × replica).
# Tune partition, mem, walltime for your cluster and for 100 ps × 1000 mol × 64 beads.
#
# NYU Greene example (adjust modules):
#   module purge
#   module load anaconda3/2023.9
#   module load cuda/12.2
#   source "$(conda info --base)/etc/profile.d/conda.sh"
#   conda activate openmm

module purge
module load anaconda3/2023.9
module load cuda/12.2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate openmm

export CAMPAIGN_MODEL=tip4p
# Optional: OPENMM_PLUGIN_DIR, CUDA_VISIBLE_DEVICES (SLURM sets when using --gres)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p logs slurm_out
bash run_campaign_task.sh
