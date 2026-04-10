#!/bin/bash
#SBATCH --job-name=rpmd100_ipi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=48GB
#SBATCH --time=48:00:00
#SBATCH --array=0-2499
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# i-PI + LAMMPS UMA RPMD — requires i-pi, LAMMPS (fix ipi), fairchem-lammps in conda env.

module purge
module load anaconda3/2023.9
module load cuda/12.2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate openmm

export CAMPAIGN_MODEL=lammps_uma

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p logs slurm_out
bash run_campaign_task.sh
