#!/bin/bash
# Run the extended cavity IR simulation (1 ns @ 300 K)

# Set up environment for installed cavity OpenMM
export LD_LIBRARY_PATH="/home/mh7373/openmm-cavity/lib:${LD_LIBRARY_PATH}"

echo "============================================================"
echo "Cavity OpenMM Extended IR Simulation"
echo "============================================================"
echo "Configuration:"
echo "  Total time: 1 ns"
echo "  Temperature: 100 K (low temp for clear vibrational peaks)"
echo "  Coupling: ON from t=0 (lambda=0.001)"
echo "  Finite-Q displacement: Applied before dynamics"
echo "  Dipole output: Every 1 fs (1 timestep)"
echo ""
echo "Output:"
echo "  Dipole trajectory: cavity_diamer_dipole.npz"
echo "  Log file: cavity_sim_output.log"
echo ""
echo "Note: This will take ~10-30 minutes on CUDA, hours on Reference"
echo "============================================================"
echo ""

# Copy the test script to /tmp to avoid any path conflicts
TEST_SCRIPT="/media/extradrive/Trajectories/openmm/tests/test_cavity_diamer.py"
WORK_DIR="/tmp/cavity_sim_$$"
mkdir -p "${WORK_DIR}"
cp "${TEST_SCRIPT}" "${WORK_DIR}/"

echo "Working directory: ${WORK_DIR}"
echo ""

# Change to work directory and run the simulation with unbuffered output
cd "${WORK_DIR}"
python -u test_cavity_diamer.py 2>&1 | tee cavity_sim_output.log

# Copy results back to test directory
echo ""
echo "Copying results back to test directory..."
cp cavity_diamer_dipole.npz /media/extradrive/Trajectories/openmm/tests/ 2>/dev/null && echo "  ✓ Dipole data copied"
cp cavity_sim_output.log /media/extradrive/Trajectories/openmm/tests/ 2>/dev/null && echo "  ✓ Log file copied"

echo ""
echo "============================================================"
echo "Simulation complete!"
echo "Results saved to: /media/extradrive/Trajectories/openmm/tests/"
echo "============================================================"

# Cleanup
cd /
rm -rf "${WORK_DIR}"
