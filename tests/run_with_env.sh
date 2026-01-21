#!/bin/bash
# Run script with custom OpenMM build

export LD_LIBRARY_PATH=/media/extradrive/Trajectories/openmm/build:$LD_LIBRARY_PATH
export OPENMM_PLUGIN_DIR=/media/extradrive/Trajectories/openmm/build
export PYTHONPATH=/media/extradrive/Trajectories/openmm/build/python/build/lib.linux-x86_64-cpython-312:$PYTHONPATH

cd /media/extradrive/Trajectories/openmm/tests
python3 "$@"
