#!/bin/bash
set -e

echo "========================================"
echo "Building OpenMM with RPMD PILE_G support"
echo "========================================"

# Create build directory if it doesn't exist
if [ ! -d "build" ]; then
    mkdir build
fi

cd build

# Configure with CMake
echo ""
echo "[1/3] Configuring with CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON \
    -DOPENMM_BUILD_CUDA_LIB=ON \
    -DOPENMM_BUILD_OPENCL_LIB=ON \
    -DOPENMM_BUILD_RPMD_PLUGIN=ON

# Build (using 4 cores)
echo ""
echo "[2/3] Building..."
make -j4

# Install (may require sudo)
echo ""
echo "[3/3] Installing..."
make install

# Update Python path
export PYTHONPATH=$(pwd)/python:$PYTHONPATH

echo ""
echo "========================================"
echo "Build completed successfully!"
echo "========================================"
echo ""
echo "To test, run:"
echo "  cd .."
echo "  python3 test_rpmd_pile_g.py"
