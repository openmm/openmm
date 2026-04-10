#!/bin/bash
# Full install into conda env: build OpenMM from this repo (cmake, make, install
# including Python bindings with integrated openmmml), then install ML runtime
# dependencies (fairchem-core) via pip.
#
# Run with env active, e.g.:
#   conda activate fairchem
#   bash /path/to/scripts/install_openmm_fairchem_base.sh
#
# -----------------------------------------------------------------------------
# Optional environment (CUDA / driver)
# -----------------------------------------------------------------------------
# OPENMM_CLEAN_BUILD=1     Remove $REPO/build before CMake (no stale CUDA cache).
# CMAKE_CUDA_COMPILER=PATH  Pass -DCMAKE_CUDA_COMPILER=PATH (pin nvcc).
# CUDAToolkit_ROOT=PATH     Pass -DCUDAToolkit_ROOT=PATH if CMake picks wrong toolkit.
# OPENMM_CUDA_ARCHITECTURES  e.g. 89 for RTX 4070 (Ada). If unset, set from nvidia-smi compute_cap.
# OPENMM_SKIP_CUDA_SMOKE=1 Skip Step 7 minimal CUDA Context test (not recommended).
# OPENMM_NO_AUTO_CUDA_PIN=1  Do not auto-pick /usr/local/cuda-12.x for nvcc (even in legacy mode).
#
# Recommended stack (typical fairchem + PyTorch): NVIDIA driver 550+ (e.g. 580, nvidia-smi may
# show “CUDA Version: 13.0” as capability). No extra nvcc/NVRTC tweaks needed.
#
# Legacy (driver series < 550, e.g. 535 on Ubuntu): conda’s libnvrtc + OpenMM JIT can hit PTX 222.
# This script then auto-pins nvcc to /usr/local/cuda-12.{0,1,2} if present and sets LD_PRELOAD
# for that toolkit’s libnvrtc.so.12 during Step 7 (conda python RPATH wins over LD_LIBRARY_PATH).
# OPENMM_FORCE_LEGACY_CUDA=1 forces those workarounds; OPENMM_FORCE_LEGACY_CUDA=0 disables them.
# OPENMM_SKIP_NVRTC_LD_PREPEND=1 skips LD_PRELOAD only (legacy nvcc pin may still apply).
#
# If Step 7 fails with “No compatible CUDA device”, see scripts/smoke_openmm_cuda.py (OOM vs PTX).
# JOBS=N                    Parallel make jobs (default 4).
# -----------------------------------------------------------------------------

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

_nvidia_drv_ver=""
_nvidia_drv_major=""
if command -v nvidia-smi >/dev/null 2>&1; then
    _nvidia_drv_ver=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d '[:space:]')
    if [ -n "$_nvidia_drv_ver" ]; then
        _nvidia_drv_major=$(echo "$_nvidia_drv_ver" | cut -d. -f1 | tr -cd '0-9')
    fi
fi

# Legacy workarounds: old drivers cannot JIT PTX from conda’s NVRTC / nvcc stack.
_OPENMM_LEGACY_CUDA=0
if [ "${OPENMM_FORCE_LEGACY_CUDA:-}" = "1" ]; then
    _OPENMM_LEGACY_CUDA=1
elif [ "${OPENMM_FORCE_LEGACY_CUDA:-}" = "0" ]; then
    _OPENMM_LEGACY_CUDA=0
elif [ -n "$_nvidia_drv_major" ] && [ "$_nvidia_drv_major" -lt 550 ] 2>/dev/null; then
    _OPENMM_LEGACY_CUDA=1
fi

# Pin nvcc to system CUDA 12.0–12.2 only in legacy mode (newer drivers: use CMake / PATH nvcc).
if [ "$_OPENMM_LEGACY_CUDA" = "1" ] && [ -z "${CMAKE_CUDA_COMPILER:-}" ] && [ "${OPENMM_NO_AUTO_CUDA_PIN:-0}" != "1" ]; then
    for cand in /usr/local/cuda-12.0 /usr/local/cuda-12.1 /usr/local/cuda-12.2; do
        if [ -x "$cand/bin/nvcc" ]; then
            export CMAKE_CUDA_COMPILER="$cand/bin/nvcc"
            export CUDAToolkit_ROOT="$cand"
            break
        fi
    done
fi
if [ -z "${OPENMM_CUDA_ARCHITECTURES:-}" ] && command -v nvidia-smi >/dev/null 2>&1; then
    _cap=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d '.')
    if [ -n "$_cap" ]; then
        export OPENMM_CUDA_ARCHITECTURES="$_cap"
    fi
fi

echo "=============================================="
echo "  OpenMM (from repo) + ML dependencies"
echo "=============================================="
echo "Target env: ${CONDA_DEFAULT_ENV:-base}"
echo "Install prefix: $INSTALL_PREFIX"
echo "OPENMM_CLEAN_BUILD: ${OPENMM_CLEAN_BUILD:-0}"
echo "NVIDIA driver: ${_nvidia_drv_ver:-<unknown>}${_nvidia_drv_major:+ (major $_nvidia_drv_major)}"
echo "Legacy CUDA workarounds (CUDA 12.x nvcc pin + Step 7 NVRTC LD_PRELOAD): $([ "$_OPENMM_LEGACY_CUDA" = "1" ] && echo on || echo off)"
echo "CMAKE_CUDA_COMPILER: ${CMAKE_CUDA_COMPILER:-<CMake default from PATH>}"
echo "CUDAToolkit_ROOT: ${CUDAToolkit_ROOT:-<unset>}"
echo "OPENMM_CUDA_ARCHITECTURES: ${OPENMM_CUDA_ARCHITECTURES:-<CMake default>}"
if [ -n "${CMAKE_CUDA_COMPILER:-}" ]; then
    echo "  (nvcc: $($CMAKE_CUDA_COMPILER --version 2>/dev/null | tail -1 || echo '?'))"
fi
echo ""

echo "Step 0/7: Prereqs (numpy in Fairchem range, PyTorch via conda)"
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

echo "Step 1/7: CMake configure (OpenMM from this repo)"
echo "----------------------------------------------"
if [ "${OPENMM_CLEAN_BUILD:-0}" = "1" ]; then
    echo "  Removing build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
fi
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

CMAKE_EXTRA=()
if [ -n "${CMAKE_CUDA_COMPILER:-}" ]; then
    CMAKE_EXTRA+=("-DCMAKE_CUDA_COMPILER=$CMAKE_CUDA_COMPILER")
    echo "  Using CMAKE_CUDA_COMPILER=$CMAKE_CUDA_COMPILER"
fi
if [ -n "${CUDAToolkit_ROOT:-}" ]; then
    CMAKE_EXTRA+=("-DCUDAToolkit_ROOT=$CUDAToolkit_ROOT")
    echo "  Using CUDAToolkit_ROOT=$CUDAToolkit_ROOT"
fi
if [ -n "${OPENMM_CUDA_ARCHITECTURES:-}" ]; then
    CMAKE_EXTRA+=("-DCMAKE_CUDA_ARCHITECTURES=$OPENMM_CUDA_ARCHITECTURES")
    echo "  Using CMAKE_CUDA_ARCHITECTURES=$OPENMM_CUDA_ARCHITECTURES"
fi

cmake "$REPO" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON \
    -DOPENMM_BUILD_CUDA_LIB=ON \
    -DOPENMM_BUILD_OPENCL_LIB=ON \
    -DOPENMM_BUILD_RPMD_PLUGIN=ON \
    "${CMAKE_EXTRA[@]}"

echo ""
echo "  --- CMake CUDA-related cache (from $BUILD_DIR/CMakeCache.txt) ---"
if [ -f CMakeCache.txt ]; then
    grep -E '^(CMAKE_CUDA_COMPILER|CMAKE_CUDA_ARCHITECTURES|CUDAToolkit_ROOT|CUDAToolkit_BIN_DIR):' CMakeCache.txt 2>/dev/null || true
else
    echo "  (CMakeCache.txt not found)"
fi
echo "  --- compare with: $(command -v nvcc 2>/dev/null || echo 'nvcc not in PATH') ---"
if command -v nvcc >/dev/null 2>&1; then
    nvcc --version | head -n 4
fi
echo ""

echo "Step 2/7: Build OpenMM (make)"
echo "----------------------------------------------"
make -j"${JOBS:-4}"

echo ""
echo "Step 3/7: Install OpenMM libs and headers (make install)"
echo "----------------------------------------------"
make install

echo ""
echo "Step 4/7: Install OpenMM Python bindings (includes openmmml) into this env"
echo "----------------------------------------------"
# setup.py requires OPENMM_LIB_PATH and OPENMM_INCLUDE_PATH (see wrappers/python/setup.py)
export OPENMM_LIB_PATH="$INSTALL_PREFIX/lib"
export OPENMM_INCLUDE_PATH="$INSTALL_PREFIX/include"
"$PYTHON" -m pip install "$BUILD_DIR/python"

echo ""
echo "Step 5/7: Install ML runtime dependencies via pip"
echo "----------------------------------------------"
# fairchem-core (standard PyPI release for UMA models)
echo "Installing fairchem-core..."
"$PYTHON" -m pip install "fairchem-core>=2.14"


echo ""
echo "Step 6/7: Verify imports"
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
echo "Step 7/7: CUDA smoke test (minimal Context; catches PTX/driver mismatch)"
echo "----------------------------------------------"
if [ "${OPENMM_SKIP_CUDA_SMOKE:-0}" = "1" ]; then
    echo "  Skipped (OPENMM_SKIP_CUDA_SMOKE=1)."
else
    export OPENMM_PLUGIN_DIR="${INSTALL_PREFIX}/lib/plugins"
    # Legacy only: conda python RPATH → conda's libnvrtc before LD_LIBRARY_PATH; LD_PRELOAD fixes PTX 222.
    if [ "$_OPENMM_LEGACY_CUDA" = "1" ] && [ "${OPENMM_SKIP_NVRTC_LD_PREPEND:-0}" != "1" ]; then
        _nvrtc_so=""
        for cand in /usr/local/cuda-12.0 /usr/local/cuda-12.1 /usr/local/cuda-12.2; do
            if [ -f "$cand/lib64/libnvrtc.so.12" ]; then
                _nvrtc_so="$cand/lib64/libnvrtc.so.12"
                break
            fi
        done
        if [ -z "$_nvrtc_so" ] && [ -n "${CUDAToolkit_ROOT:-}" ] && [ -f "${CUDAToolkit_ROOT}/lib64/libnvrtc.so.12" ]; then
            _nvrtc_so="${CUDAToolkit_ROOT}/lib64/libnvrtc.so.12"
        fi
        if [ -n "$_nvrtc_so" ]; then
            export LD_PRELOAD="$_nvrtc_so${LD_PRELOAD:+:$LD_PRELOAD}"
            echo "  LD_PRELOAD=$_nvrtc_so (legacy driver: NVRTC from system CUDA 12.x)."
        fi
    fi
    if ! "$PYTHON" "$REPO/scripts/smoke_openmm_cuda.py"; then
        echo ""
        echo "CUDA smoke test failed. OpenMM Python is at:"
        "$PYTHON" -c "import openmm; print(' ', openmm.__file__)" || true
        echo "Plugins expected under: $OPENMM_PLUGIN_DIR"
        echo "If PTX 222 / wrong NVRTC: upgrade to NVIDIA driver 550+ (e.g. 580), or on older drivers run:"
        echo "  export LD_PRELOAD=/usr/local/cuda-12.0/lib64/libnvrtc.so.12\${LD_PRELOAD:+:\$LD_PRELOAD}"
        echo "  python $REPO/scripts/smoke_openmm_cuda.py"
        echo "Or OPENMM_FORCE_LEGACY_CUDA=1 bash this script to force CUDA-12.x pin + LD_PRELOAD."
        exit 1
    fi
fi

echo ""
echo "=============================================="
echo "  Install complete"
echo "=============================================="
echo "OpenMM was built from this repo and installed into this env (prefix: $INSTALL_PREFIX)."
echo "openmmml is now integrated into OpenMM (no separate openmm-ml install needed)."
echo "ML runtime dep (fairchem-core) is pip-installed."
echo ""
echo "Runtime check: python -c \"import openmm; print(openmm.__file__)\""
echo "Should live under this prefix / site-packages for the active env."
echo ""
echo "To reinstall ML deps only:  pip install -r $REPO/requirements-ml.txt"
echo ""
echo "Driver <550 + PTX 222: export LD_PRELOAD=/usr/local/cuda-12.0/lib64/libnvrtc.so.12 before Python,"
echo "or upgrade to driver 550+ (recommended with conda PyTorch)."
echo ""
