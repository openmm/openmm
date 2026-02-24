#!/bin/bash
# Full install into conda base: build OpenMM from this repo (cmake, make, install
# including Python bindings with integrated openmmml), then install ML runtime
# dependencies (fairchem-core) via pip.
#
# One script, full rebuild. Installs OpenMM (with openmmml baked in), fairchem-core.
#
# Run with base active:
#   conda activate base
#   bash /path/to/scripts/install_openmm_fairchem_base.sh

set -e

# Repo root: parent of scripts/ (portable across machines)
REPO="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="$REPO/build"

# Install prefix: use conda env if active, else default
INSTALL_PREFIX="${CONDA_PREFIX:-/usr/local/openmm}"

# Use python from PATH (conda env provides python; fallback to python3)
PYTHON=$(command -v python 2>/dev/null || command -v python3 2>/dev/null)
if [ -z "$PYTHON" ]; then
    echo "ERROR: python or python3 not found in PATH. Activate your conda env and re-run."
    exit 1
fi

echo "=============================================="
echo "  OpenMM (from repo) + ML dependencies"
echo "=============================================="
echo "Target env: ${CONDA_DEFAULT_ENV:-base}"
echo "Install prefix: $INSTALL_PREFIX"
echo ""

echo "Step 0/6: Prereqs (numpy in Fairchem range, PyTorch via conda)"
echo "----------------------------------------------"
# fairchem-core requires numpy>=2.0,<2.3. OpenMM is built with whatever numpy is present;
# use Fairchem's range up front so we never get NumPy 1.x at build and 2.x at run.
"$PYTHON" -m pip install 'numpy>=2.0,<2.3' -q 2>/dev/null || true
"$PYTHON" -c "
import numpy
v = tuple(map(int, numpy.__version__.split('.')[:2]))
if not (v >= (2, 0) and v < (2, 3)):
    raise SystemExit('numpy must be >=2.0,<2.3 for fairchem-core; got %s. Fix env and re-run.' % numpy.__version__)
" || exit 1
if [ -n "$CONDA_PREFIX" ] && ! "$PYTHON" -c "import torch" 2>/dev/null; then
    echo "PyTorch not found; install via conda so pip does not replace it."
    conda install pytorch -c pytorch -y
elif [ -n "$CONDA_PREFIX" ] && "$PYTHON" -c "import torch" 2>/dev/null; then
    echo "PyTorch already present (prefer conda-managed)."
fi
echo ""

echo "Step 1/6: CMake configure (OpenMM from this repo)"
echo "----------------------------------------------"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake "$REPO" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON \
    -DOPENMM_BUILD_CUDA_LIB=ON \
    -DOPENMM_BUILD_OPENCL_LIB=ON \
    -DOPENMM_BUILD_RPMD_PLUGIN=ON

echo ""
echo "Step 2/6: Build OpenMM (make)"
echo "----------------------------------------------"
make -j"${JOBS:-4}"

echo ""
echo "Step 3/6: Install OpenMM libs and headers (make install)"
echo "----------------------------------------------"
make install

echo ""
echo "Step 4/6: Install OpenMM Python bindings (includes openmmml) into this env"
echo "----------------------------------------------"
# setup.py requires OPENMM_LIB_PATH and OPENMM_INCLUDE_PATH (see wrappers/python/setup.py)
export OPENMM_LIB_PATH="$INSTALL_PREFIX/lib"
export OPENMM_INCLUDE_PATH="$INSTALL_PREFIX/include"
"$PYTHON" -m pip install "$BUILD_DIR/python"

echo ""
echo "Step 5/6: Install ML runtime dependencies via pip"
echo "----------------------------------------------"
# fairchem-core (standard PyPI release for UMA models)
echo "Installing fairchem-core..."
"$PYTHON" -m pip install "fairchem-core>=2.14"


echo "Step 6/6: Verify imports"
echo "----------------------------------------------"
"$PYTHON" -c "
import openmm
from openmm import unit
print('  openmm:', openmm.__version__, '(', openmm.__file__, ')')
import openmm.app
print('  openmm.app: ok')
import openmmml
print('  openmmml: ok (integrated into OpenMM)')
try:
    import fairchem.core
    print('  fairchem.core: ok (pip)')
except ImportError:
    print('  fairchem.core: NOT FOUND (UMA models will not be available)')
print('  Stack OK')
" || {
    echo "One or more required imports failed (openmm, openmm.app, openmmml)."
    echo "ML backend (fairchem-core) is optional runtime dep installed via pip."
    exit 1
}

echo ""
echo "=============================================="
echo "  Install complete"
echo "=============================================="
echo "OpenMM was built from this repo and installed into this env (prefix: $INSTALL_PREFIX)."
echo "openmmml is now integrated into OpenMM (no separate openmm-ml install needed)."
echo "ML runtime dep (fairchem-core) is pip-installed."
echo ""
echo "To reinstall ML deps only:  pip install -r $REPO/requirements-ml.txt"
echo ""
