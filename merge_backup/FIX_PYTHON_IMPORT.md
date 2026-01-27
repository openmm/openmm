# Fixing "No module named 'openmm'" Error

## Problem

When running the example script, you get:
```
Error importing OpenMM: No module named 'openmm'
```

This happens because the Python bindings aren't properly installed or accessible.

## Root Cause

The Python module was built but not installed. The compiled extension (`.so` file) is in `build/python/build/lib.linux-x86_64-cpython-312/` but Python can't find it, or there's a circular import issue in the build directory.

## Solutions

### Option 1: Install OpenMM Properly (Recommended)

```bash
cd /media/extradrive/Trajectories/openmm/build

# Set environment variables needed for installation
export OPENMM_INCLUDE_PATH="/media/extradrive/Trajectories/openmm/include:/media/extradrive/Trajectories/openmm/openmmapi/include"
export OPENMM_LIB_PATH="/home/mh7373/miniconda3/lib"  # or wherever libOpenMM.so is

# Install Python bindings
cd python
python3 setup.py install --user
```

### Option 2: Use Conda OpenMM (Quick Fix)

If you have OpenMM installed via conda:

```bash
# Activate conda environment with OpenMM
conda activate <your-env-with-openmm>

# Or install OpenMM in current environment
conda install -c conda-forge openmm

# Then run the example
python3 examples/cavity/dimer_system/run_simulation_cavity_driven.py
```

### Option 3: Fix Build Directory Installation

The `make PythonInstall` command should work, but it needs proper environment variables:

```bash
cd /media/extradrive/Trajectories/openmm/build

# Make sure libraries are installed first
make install

# Then install Python bindings with proper paths
export OPENMM_INCLUDE_PATH="/media/extradrive/Trajectories/openmm/include"
export OPENMM_LIB_PATH="/home/mh7373/miniconda3/lib"  # Adjust to your install location
make PythonInstall
```

### Option 4: Use System OpenMM (If Available)

If OpenMM is installed system-wide:

```bash
# Check if it's available
python3 -c "import openmm; print(openmm.__file__)"

# If it works, just run the example
python3 examples/cavity/dimer_system/run_simulation_cavity_driven.py
```

## Verification

After installation, verify it works:

```bash
python3 -c "import openmm; print('OpenMM version:', openmm.__version__); print('Location:', openmm.__file__)"
```

## Note

The C++ implementation is complete and all tests pass. The Python import issue is purely an installation/configuration problem, not a code issue.
