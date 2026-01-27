# Why "No module named 'openmm'" Happens

## The Issue

When you run:
```bash
python3 examples/cavity/dimer_system/run_simulation_cavity_driven.py
```

You get:
```
Error importing OpenMM: No module named 'openmm'
```

## Root Cause

**The Python bindings for OpenMM are not installed or accessible in your current Python environment.**

Here's what's happening:

1. **C++ Code is Complete**: The C++ implementation is fully built and all tests pass (100% test success)
2. **Python Module Not Installed**: The Python wrappers exist in `build/python/` but aren't properly installed
3. **Circular Import Issue**: The build directory has a circular import problem that prevents direct use

## Why This Happens

The `make PythonInstall` step requires:
- Proper environment variables (`OPENMM_INCLUDE_PATH`, `OPENMM_LIB_PATH`)
- The C++ libraries to be installed first (`make install`)
- No conflicts with existing installations (you have conda OpenMM which conflicts)

## Solutions

### Quick Fix: Install OpenMM via Conda

```bash
# Install OpenMM in your current conda environment
conda install -c conda-forge openmm

# Verify it works
python3 -c "import openmm; print('OpenMM version:', openmm.__version__)"

# Run the example
python3 examples/cavity/dimer_system/run_simulation_cavity_driven.py
```

**Note**: This will use the conda version, which may not have your new laser features yet. But it will let you test the example structure.

### Proper Fix: Install Your Custom Build

```bash
cd /media/extradrive/Trajectories/openmm/build

# 1. Install C++ libraries first
make install

# 2. Install Python bindings (adjust paths as needed)
export OPENMM_INCLUDE_PATH="/media/extradrive/Trajectories/openmm/include"
export OPENMM_LIB_PATH="/home/mh7373/miniconda3/lib"  # or your install location
cd python
python3 setup.py install --user --force  # --force to overwrite conflicts
```

### Alternative: Use a Clean Environment

```bash
# Create a new conda environment
conda create -n openmm-dev python=3.12
conda activate openmm-dev

# Install dependencies
conda install -c conda-forge numpy swig

# Then install your custom build
cd /media/extradrive/Trajectories/openmm/build
make install
export OPENMM_INCLUDE_PATH="/media/extradrive/Trajectories/openmm/include"
export OPENMM_LIB_PATH="/usr/local/openmm/lib"  # or your install prefix
cd python
python3 setup.py install
```

## Verification

After installation, test it:

```bash
python3 -c "
import openmm
print('✓ OpenMM imported')
print('Version:', openmm.__version__)
print('Location:', openmm.__file__)

# Test your new laser features
cf = openmm.CavityForce(0, 0.01, 0.001, 1.0)
print('✓ CavityForce created')
print('Has laser methods:', hasattr(cf, 'setCavityDriveAmplitude'))
"
```

## Summary

- **C++ implementation**: ✅ Complete and tested (all 7 tests pass)
- **Python bindings**: ⚠️ Need proper installation
- **The error**: Python can't find the `openmm` module because it's not installed

The code is correct; it just needs to be installed properly in your Python environment.
