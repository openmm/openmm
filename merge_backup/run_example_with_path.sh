#!/bin/bash
# Run the cavity-driven example with proper PYTHONPATH setup

cd "$(dirname "$0")"

# Set up PYTHONPATH to include both build directories
export PYTHONPATH="build/python/build/lib.linux-x86_64-cpython-312:build/python:$PYTHONPATH"

echo "Running cavity-driven simulation example..."
echo "PYTHONPATH: $PYTHONPATH"
echo ""

python3 examples/cavity/dimer_system/run_simulation_cavity_driven.py
