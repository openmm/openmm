# ml-experimental Overview

This directory contains experimental features for OpenMM developed in the ml-experimental branch.

## Features

For detailed information, see:
- [ml-experimental/README.md](../../ml-experimental/README.md) - Feature overview and quick start
- [ml-experimental/ARCHITECTURE.md](ARCHITECTURE.md) - System architecture and design

## Quick Links

### Cavity Particles
- [Implementation Guide](../../ml-experimental/cavity/docs/CAVITY_IMPLEMENTATION.md)
- [Usage Guide](../../ml-experimental/cavity/docs/CAVITY_USAGE.md)
- [Examples](../../ml-experimental/cavity/examples/)

### Hybrid RPMD
- [Implementation Status](../../ml-experimental/rpmd/docs/HYBRID_RPMD_IMPLEMENTATION_COMPLETE.md)
- [Quick Start](../../ml-experimental/rpmd/docs/HYBRID_RPMD_QUICKSTART.md)
- [Examples](../../ml-experimental/rpmd/examples/)

### ML Potentials
- [UMA Implementation Summary](../../ml-experimental/ml/docs/UMA_IMPLEMENTATION_SUMMARY.md)
- [OMOL25 LES Integration](../../ml-experimental/ml/docs/OMOL25_LES_INTEGRATION_GUIDE.md)
- [Available Models](../../ml-experimental/ml/docs/AVAILABLE_UMA_MODELS.md)
- [Examples](../../ml-experimental/ml/examples/)

## Installation

Build from the ml-experimental branch:

```bash
git checkout ml-experimental
mkdir build && cd build
cmake ..
make -j8
make install
make PythonInstall
```

## Status

All features are production-ready with full GPU support and RPMD compatibility.
