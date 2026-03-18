#!/bin/bash
#SBATCH --job-name=rpmd_bench
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=48GB
#SBATCH --time=48:00:00
#SBATCH --array=0-3899
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# NYU Greene: load modules and activate conda env for openmm, openmmml, fairchem, ipi, lammps.
# Adjust module names if needed: module avail anaconda3 cuda
module purge
module load anaconda3/2023.9
module load cuda/12.2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate openmm  # or your env with openmm, openmmml, fairchem, ipi, lammps

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
bash run_rpmd_single.sh
