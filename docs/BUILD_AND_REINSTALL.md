# Build and Reinstall OpenMM (UMA + RPMD)

Full rebuild and install for testing on a new machine or after code changes.

## Prerequisites

- **Conda** (miniconda or anaconda) with Python 3.10–3.12
- **Build tools**: cmake, make, C++ compiler, SWIG, doxygen
- **Optional**: CUDA for GPU (script enables it by default)
- **HuggingFace**: UMA models require a gated repo; run `huggingface-cli login` once before using UMA

### Install build dependencies (Linux)

```bash
# Debian/Ubuntu
sudo apt-get install -y cmake build-essential swig libfftw3-dev libxml2-dev doxygen zlib1g-dev

# Or via conda
conda install -c conda-forge cmake make swig
```

## Quick Rebuild (same machine)

From the repo root:

```bash
conda activate base   # or your env
bash scripts/install_openmm_fairchem_base.sh
```

This builds OpenMM, installs libs and Python bindings, and pip-installs fairchem-core.

## Rebuild on a Different Computer

### 1. Get the repo

```bash
git clone <your-repo-url> openmm
cd openmm
# Or sync/copy the repo if you use rsync or another method
```

### 2. Edit the install script (if repo path differs)

The script uses `REPO="/media/extradrive/Trajectories/openmm"` by default. Either:

- Change `REPO` in `scripts/install_openmm_fairchem_base.sh` to your repo path, or
- Use the portable version below.

### 3. Run the full install

```bash
conda activate base
cd /path/to/openmm
bash scripts/install_openmm_fairchem_base.sh
```

### 4. Set plugin dir (for RPMD)

Plugins (including RPMD) are installed to `$CONDA_PREFIX/lib/plugins`. The test script tries `~/miniconda3/lib/plugins` if `OPENMM_PLUGIN_DIR` is unset. If your env uses a different path:

```bash
export OPENMM_PLUGIN_DIR="$CONDA_PREFIX/lib/plugins"
```

Add to `~/.bashrc` if needed.

### 5. HuggingFace (for UMA)

UMA models are gated. Log in once:

```bash
huggingface-cli login
# Or: python -c "from huggingface_hub import login; login()"
```

Accept the model agreement at https://huggingface.co/facebook/UMA.

## Manual Build (if the script fails)

```bash
cd /path/to/openmm

# 1. Configure
mkdir -p build && cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -DCMAKE_BUILD_TYPE=Release \
  -DOPENMM_BUILD_CUDA_LIB=ON \
  -DOPENMM_BUILD_RPMD_PLUGIN=ON

# 2. Build
make -j$(nproc)

# 3. Install libs
make install

# 4. Install Python bindings
export OPENMM_LIB_PATH="$CONDA_PREFIX/lib"
export OPENMM_INCLUDE_PATH="$CONDA_PREFIX/include"
pip install python/

# 5. ML deps
pip install -r requirements-ml.txt
```

## Python-only changes (no rebuild)

If you only modified Python files (e.g. `umapotential_pythonforce_batch.py`), you can copy them into the installed package instead of rebuilding:

```bash
INSTALL_PATH=$(python -c "import openmmml; import os; print(os.path.dirname(openmmml.__file__))")
cp wrappers/python/openmmml/models/umapotential_pythonforce_batch.py "$INSTALL_PATH/models/"
cp wrappers/python/openmmml/models/umapotential_pythonforce.py "$INSTALL_PATH/models/"
```

## Run the ice RPMD test

```bash
cd tests/uma_ice_rpmd

# CPU (avoids CUDA context conflict)
python test_uma_ice_rpmd.py --molecules 16 --beads 4 --temperature 243 --dt 1.0 --equil 5 --prod 50 --platform cpu

# CUDA (uses lazy-load for single-GPU; add --ml-device cpu if needed)
python test_uma_ice_rpmd.py --molecules 16 --beads 4 --temperature 243 --dt 1.0 --equil 5 --prod 50 --platform cuda
```
