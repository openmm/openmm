[![GH Actions Status](https://github.com/openmm/openmm/workflows/CI/badge.svg)](https://github.com/openmm/openmm/actions?query=branch%3Amaster+workflow%3ACI)
[![Conda](https://img.shields.io/conda/v/conda-forge/openmm.svg)](https://anaconda.org/conda-forge/openmm)
[![Anaconda Cloud Badge](https://anaconda.org/conda-forge/openmm/badges/downloads.svg)](https://anaconda.org/conda-forge/openmm)

## OpenMM-LM: A High Performance Molecular Dynamics Library with Light-Matter Interactions

### Introduction

[OpenMM](https://openmm.org) is a toolkit for molecular simulation. It can be used either as a stand-alone application for running simulations, or as a library you call from your own code. It
provides a combination of extreme flexibility (through custom forces and integrators), openness, and high performance (especially on recent GPUs) that make it truly unique among simulation codes.
This version adds light-matter interactions, particularly those used in optical cavities, additional modules for nuclear quantum effects (ported from i-PI), and interfaces to SOTA machine learning interatomic potentials, e.g., FAIR Chemistry's UMA and AIMNet2. 

### Getting Help

Need help using OpenMM? There are several places that you can find it:

- [Documentation](https://docs.openmm.org/):
  - [User Manual](https://docs.openmm.org/latest/userguide/)
  - [Python API Reference](https://docs.openmm.org/latest/api-python/)
  - [C++ API Reference](https://docs.openmm.org/latest/api-c++/)
- [Getting Support](SUPPORT.md):
  - [Frequently Asked Questions](https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions)
  - [Discussion Forum](https://github.com/openmm/openmm/discussions)
  - [Issue Tracker](https://github.com/openmm/openmm/issues)
- [Learning Resources](https://openmm.github.io/openmm-cookbook/latest):
  - [OpenMM Tutorials](https://openmm.github.io/openmm-cookbook/latest/tutorials)
  - [OpenMM Cookbook](https://openmm.github.io/openmm-cookbook/latest/cookbook)
  - [OpenMM Examples](examples/README.md)
- [Contributing to OpenMM](CONTRIBUTING.md):
  - [Developer Guide](https://docs.openmm.org/latest/developerguide/)

### License Information

OpenMM is free and open-source software. There are several licenses which cover
different parts of OpenMM, but most of the source code is covered by the MIT
license or the GNU Lesser General Public License (LGPL). Portions copyright
© 2008-2025 Stanford University and the Authors. For more details, see
[Licenses.txt](docs-source/licenses/Licenses.txt).

---

## About this fork (ML / UMA / RPMD)

This repository extends upstream OpenMM with **machine-learning potentials** (bundled as **`openmmml`** in the Python install), **`PythonForce`** with **batched** evaluation for **RPMD**, and related **RPMD plugin** fixes (additive forces, hybrid classical–quantum RPMD, barostat coordination). Details:

| Topic | Where |
|--------|--------|
| RPMD + UMA bug fixes and regression context | [FIXES_SUMMARY.md](FIXES_SUMMARY.md) |
| Hybrid RPMD design | [plugins/rpmd/HYBRID_RPMD.md](plugins/rpmd/HYBRID_RPMD.md) |
| Rebuild / install (CUDA, FairChem, Hugging Face) | [docs/BUILD_AND_REINSTALL.md](docs/BUILD_AND_REINSTALL.md) |
| ML runtime pip deps (FairChem, PyTorch) | [requirements-ml.txt](requirements-ml.txt) |
| Ice Ih RPMD benchmarks and LAMMPS / i-PI parity | [tests/uma_ice_rpmd/README.md](tests/uma_ice_rpmd/README.md) |
| RPMD-focused tests | [tests/rpmd/README.md](tests/rpmd/README.md) |
| Experimental ML architecture notes | [docs/ml-experimental/README.md](docs/ml-experimental/README.md) |

**Support for this fork:** use **this repository’s** issues and discussions for UMA/RPMD/`PythonForce` batching and build questions; use **upstream** OpenMM channels for general API questions (see [SUPPORT.md](SUPPORT.md)).

### Submodule layout (FairChem)

Downstream **FairChem** work may point submodules at a personal fork; no change is required if you already use `https://github.com/muhammadhasyim/fairchem.git` (e.g. branch `les_branch`).
