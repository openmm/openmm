#!/usr/bin/env bash
# LAMMPS (pip) needs MPICH's libmpi.so.12 — not OpenMPI's .40
# Option A: conda install mpich  →  $CONDA_PREFIX/lib/libmpi.so.12
# Option B: point to an env that already has MPICH, e.g. fenics:
export LD_LIBRARY_PATH="${LAMMPS_MPI_LIB:-$CONDA_PREFIX/lib}:/home/mh7373/miniconda3/envs/fenics/lib:${LD_LIBRARY_PATH}"

set -e
cd "$(dirname "$0")"

# Quick 8-molecule ~0.7 ps (default)
DEVICE="${DEVICE:-cuda}"
MOL="${MOL:-8}"
INFILE="${INFILE:-in.ice_uma_quick.lmp}"

echo "LD_LIBRARY_PATH (first dir should provide libmpi.so.12): $LD_LIBRARY_PATH"
python run_lammps_uma_ice.py --molecules "$MOL" --device "$DEVICE" --infile "$INFILE"
echo "Trajectory: $(pwd)/dump.ice_uma_quick.lammpstrj (or dump name in $INFILE)"
