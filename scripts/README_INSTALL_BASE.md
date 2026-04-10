# Install OpenMM (from repo) + Fairchem + OpenMM-ML into base

One script, **full rebuild**. No separate “resume” or “fix” scripts. Run [scripts/install_openmm_fairchem_base.sh](install_openmm_fairchem_base.sh) from scratch; it builds OpenMM from this repo and installs fairchem + openmm-ml into the active conda env (typically base). No conda openmm, no PYTHONPATH.

## What the script does

0. **Prereqs** – Ensures numpy in Fairchem range (`>=2.0,<2.3`) and PyTorch via conda (if in conda and missing), so the build and later runs use one consistent stack.
1. **CMake configure** – OpenMM from this repo, install prefix = `$CONDA_PREFIX` (or `/usr/local/openmm` if not in conda).
2. **make** – Builds OpenMM (core, plugins including RPMD, Python bindings).
3. **make install** – Installs OpenMM libs and headers into that prefix.
4. **pip install build/python** – Installs the OpenMM Python package so `import openmm` works.
5. **pip install -e fairchem/packages/fairchem-core[dev]** – Editable fairchem-core.
6. **pip install -e fairchem/openmm-ml** – Editable openmm-ml.

Verification: **openmm**, **openmm.app**, and **fairchem.core** are required; **openmmml** is tried.

## Run it

```bash
conda activate base
bash scripts/install_openmm_fairchem_base.sh
```

From elsewhere, use the full path to `scripts/install_openmm_fairchem_base.sh`.

## Requirements

- Conda env active (script uses `$CONDA_PREFIX` as install prefix when set).
- Build tools: cmake, make, C++ compiler, SWIG, Python development headers.
- Step 0 will install numpy in Fairchem’s range and, when in conda and PyTorch is missing, run `conda install pytorch -c pytorch`. If you already use PyTorch, prefer a conda-managed install so pip does not replace it in Step 5.
- Optional: CUDA/OpenCL (script enables them by default; edit to turn off if needed).

## Dependency conflicts (OpenMM vs Fairchem)

| Package        | numpy            | ase              |
|----------------|------------------|------------------|
| OpenMM         | any (no version) | —                |
| Fairchem-core  | `>=2.0,<2.3`     | `>=3.26.0`       |

OpenMM does not depend on Fairchem.

## openmm-torch (for UMA/ML in OpenMM)

This script does **not** install openmm or openmm-torch via conda. For UMA/ML forces in OpenMM you need **openmm-torch** built from source against this OpenMM. Use the [openmm-torch](https://github.com/openmm/openmm-torch) repo and point it at this install.

## After install

```bash
python /path/to/openmm/ml-experimental/ml/examples/uma_ice_rpmd/test_uma_ice_rpmd.py -h
```

## Troubleshooting

- **CMake or build fails** – Install build deps (e.g. `conda install -c conda-forge cmake make swig python numpy`). To disable CUDA/OpenCL, edit the script and set `-DOPENMM_BUILD_CUDA_LIB=OFF` / `-DOPENMM_BUILD_OPENCL_LIB=OFF`.
- **“No module named 'openmm'”** – Step 4 installs from `build/python`. Re-run the **full script** (no resume scripts); ensure Step 0 ran so numpy is in Fairchem’s range before the build.
- **Torch broken during or after Step 5** – Use one source for torch (conda or pip). If broken: `pip uninstall torch torchvision torchaudio -y` then `conda install pytorch -c pytorch`. Then re-run the **full install script** from the start.
- **openmmml fails** – Ensure fairchem-core is installed. If import errors persist, try `pip install --force-reinstall scipy`.
- **UMA/ML in OpenMM** – Build and install openmm-torch from source against this OpenMM; do not use conda’s openmm or openmm-torch.
