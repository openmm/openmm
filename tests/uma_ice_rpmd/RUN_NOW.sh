#!/usr/bin/env bash
# One place to run OpenMM or LAMMPS UMA ice — copy/paste blocks as needed.

set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "========== OpenMM + UMA (no MPI; use if LAMMPS fails) =========="
echo "CPU (~few minutes for 8 molecules × 0.4 ps):"
echo "  cd $ROOT"
echo "  python test_uma_ice_rpmd.py --input ice.cif --molecules 8 --beads 1 \\"
echo "    --model uma-s-1p1-pythonforce-batch --temperature 243 --dt 1.0 \\"
echo "    --equil 0.1 --prod 0.3 --pressure 0 --pdb-interval 0.1 \\"
echo "    --output run_now --minimal --platform cpu --ml-device cpu"
echo "  → Trajectory: run_now/ice_rpmd_uma_T243_b1.pdb"
echo ""

echo "========== LAMMPS + UMA =========="
echo "Needs: pip install lammps fairchem-lammps"
echo "pip's LAMMPS wants libmpi.so.12 (MPICH) and often libgfortran.so.3."
echo ""
echo "1) libgfortran.so.3 (Ubuntu 24+ has NO apt libgfortran3):"
echo "     conda install -y -c conda-forge 'libgfortran=3.0.0'"
echo "   Or replace old mpich: conda remove mpich --force && conda install -c conda-forge mpich"
echo ""
echo "2) Then:"
echo "     export LD_LIBRARY_PATH=\"\$CONDA_PREFIX/lib:\$HOME/miniconda3/lib/python3.12/site-packages/lammps.libs:\${LD_LIBRARY_PATH}\""
echo "     cd $ROOT/lammps"
echo "     python run_lammps_uma_ice.py --molecules 8 --device cuda --infile in.ice_uma_quick.lmp"
echo "  → Trajectory: lammps/dump.ice_uma_quick.lammpstrj"
echo ""

if [[ "${1:-}" == "openmm" ]]; then
  python test_uma_ice_rpmd.py --input ice.cif --molecules 8 --beads 1 \
    --model uma-s-1p1-pythonforce-batch --temperature 243 --dt 1.0 \
    --equil 0.1 --prod 0.3 --pressure 0 --pdb-interval 0.1 \
    --output run_now --minimal --platform cpu --ml-device cpu
  echo "Done: $ROOT/run_now/ice_rpmd_uma_T243_b1.pdb"
fi
