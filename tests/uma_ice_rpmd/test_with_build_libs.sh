#!/bin/bash
# Script to run test with build directory libraries (before sudo make install)
# This sets up LD_LIBRARY_PATH and PYTHONPATH to use the build directory

BUILD_DIR="/media/extradrive/Trajectories/openmm/build"

# Set library path to find libOpenMM.so and other OpenMM libraries
export LD_LIBRARY_PATH="${BUILD_DIR}:${LD_LIBRARY_PATH}"

# Set Python path to find OpenMM Python bindings
export PYTHONPATH="${BUILD_DIR}/python:${PYTHONPATH}"

# Copy _openmm.so to the right location if needed
if [ ! -f "${BUILD_DIR}/python/openmm/_openmm.so" ] && [ -f "${BUILD_DIR}/python/build/lib.linux-x86_64-cpython-312/openmm/_openmm.cpython-312-x86_64-linux-gnu.so" ]; then
    cp "${BUILD_DIR}/python/build/lib.linux-x86_64-cpython-312/openmm/_openmm.cpython-312-x86_64-linux-gnu.so" \
       "${BUILD_DIR}/python/openmm/_openmm.so" 2>/dev/null || true
fi

cd /media/extradrive/Trajectories/openmm/tests/uma_ice_rpmd

echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
echo "PYTHONPATH=${PYTHONPATH}"
echo "Running test..."
echo ""

python test_uma_ice_rpmd.py --molecules 10 --beads 8 --temperature 243 --dt 1.0 --equil 10 --prod 100 --pressure 0 --model uma-s-1-pythonforce-batch "$@"
