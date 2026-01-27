#!/bin/bash
# Helper script to set up PYTHONPATH for using OpenMM from build directory
# Usage: source setup_python_path.sh

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_PYTHON="$SCRIPT_DIR/build/python"
BUILD_LIB="$SCRIPT_DIR/build/python/build/lib.linux-x86_64-cpython-312"

# Add both the package directory and the build lib directory
export PYTHONPATH="$BUILD_LIB:$BUILD_PYTHON:$PYTHONPATH"

echo "✓ PYTHONPATH set to:"
echo "  $BUILD_LIB"
echo "  $BUILD_PYTHON"
echo ""
echo "You can now run: python3 examples/cavity/dimer_system/run_simulation_cavity_driven.py"
