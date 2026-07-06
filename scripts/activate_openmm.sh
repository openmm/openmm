#!/usr/bin/env bash
# Source before Python / Jupyter:  source scripts/activate_openmm.sh
export OPENMM_PLUGIN_DIR="${OPENMM_PLUGIN_DIR:-/usr/local/openmm/lib/plugins}"
export OPENMM_LIB_PATH="${OPENMM_LIB_PATH:-/usr/local/openmm/lib}"
export OPENMM_INCLUDE_PATH="${OPENMM_INCLUDE_PATH:-/usr/local/openmm/include}"
export LD_LIBRARY_PATH="${OPENMM_LIB_PATH}:${LD_LIBRARY_PATH:-}"
