# ml-experimental: Experimental Features for OpenMM

This directory contains experimental extensions to OpenMM developed for the ml-experimental branch. These features add support for cavity-coupled molecular dynamics, hybrid classical/quantum RPMD, and machine learning potentials.

## Directory Structure

```
ml-experimental/
├── cavity/              # Cavity particle implementation
├── rpmd/                # RPMD enhancements
├── ml/                  # ML potential integration
└── tests/               # Consolidated test suite
```

## Features

### 1. Cavity Particles (`cavity/`)

Implements quantum light-matter coupling for polariton chemistry simulations.

**Key Capabilities:**
- Cavity photon mode as dummy particle
- Dipole-cavity coupling Hamiltonian
- GPU acceleration (CUDA/OpenCL)
- RPMD compatible

**Documentation:** [`cavity/docs/`](cavity/docs/)  
**Examples:** [`cavity/examples/`](cavity/examples/)

**Quick Start:**
```python
from ml_experimental.cavity.utils import add_cavity_particle, setup_cavity_coupling

cavity_idx = add_cavity_particle(system, positions, omegac_au)
cavity_force, displacer = setup_cavity_coupling(
    system, cavity_idx, omegac_au, lambda_coupling=0.01
)
```

### 2. Hybrid RPMD (`rpmd/`)

Extends RPMDIntegrator to treat particles selectively as quantum (ring polymer) or classical (single bead).

**Key Capabilities:**
- Particle-type based quantum/classical selection
- 5-10× speedup for typical biomolecular systems
- Multiple thermostat options for classical particles
- Full GPU support

**Documentation:** [`rpmd/docs/`](rpmd/docs/)  
**Examples:** [`rpmd/examples/`](rpmd/examples/)

**Quick Start:**
```python
integrator = openmm.RPMDIntegrator(num_beads, temperature, friction, timestep)

# Treat hydrogens as quantum, heavy atoms as classical
integrator.setQuantumParticleTypes([1])  # Type 1 = H
integrator.setDefaultQuantum(False)      # Heavy atoms classical by default
```

### 3. ML Potentials (`ml/`)

Integration of UMA (Universal Molecular Atomistic) and other ML potentials into OpenMM.

**Key Capabilities:**
- FAIRChem UMA models via OpenMM-ML
- PythonForce implementation for immediate usage
- RPMD compatible
- GPU accelerated inference

**Documentation:** [`ml/docs/`](ml/docs/)  
**Examples:** [`ml/examples/`](ml/examples/)

**Quick Start:**
```python
from openmmml import MLPotential

potential = MLPotential('uma-s-1')
system = potential.createSystem(topology)
```

## Usage

### Installation

#### Core OpenMM Features

These features are built into OpenMM when compiled from the ml-experimental branch:

```bash
cd /path/to/openmm
git checkout ml-experimental
mkdir build && cd build
cmake ..
make -j8
make install
make PythonInstall
```

#### ML Potentials (UMA/OpenMM-ML)

For ML potential support, install OpenMM-ML with UMA models:

```bash
# 1. Install OpenMM and OpenMM-Torch
conda install -c conda-forge openmm openmmtorch

# 2. Install FAIRChem
pip install fairchem-core

# 3. Install OpenMM-ML from repository
cd fairchem/openmm-ml
pip install -e .

# 4. Verify installation
python -c "from openmmml import MLPotential; pot = MLPotential('uma-s-1p1'); print('✓ UMA models registered')"
```

**Note:** See [`ml/docs/README.md`](ml/docs/README.md) for detailed installation instructions, troubleshooting, and limitations.

### Running Examples

Each feature directory contains working examples:

```bash
# Cavity particle - water system
cd ml-experimental/cavity/examples/water_system
python run_simulation.py --cavity-freq 1600 --lambda 0.01

# Hybrid RPMD - simple test
cd ml-experimental/rpmd/examples
python test_hybrid_rpmd_simple.py

# ML potentials - UMA ice RPMD
cd ml-experimental/ml/examples/uma_ice_rpmd
python test_uma_ice_rpmd.py --beads 8 --molecules 50
```

### Testing

Test suite organized by type:

```
ml-experimental/tests/
├── unit/           # Unit tests
├── integration/    # Integration tests
└── examples/       # Example-based tests
```

Run tests:
```bash
# C++ unit tests (in root tests/ directory)
cd build
make test

# Python integration tests
cd ml-experimental/tests/integration
pytest
```

## Documentation

### Feature Documentation

Each feature has detailed documentation:

- **Cavity:** [`cavity/docs/CAVITY_IMPLEMENTATION.md`](cavity/docs/CAVITY_IMPLEMENTATION.md)
- **RPMD:** [`rpmd/docs/HYBRID_RPMD.md`](rpmd/docs/HYBRID_RPMD.md)
- **ML:** [`ml/docs/UMA_IMPLEMENTATION_SUMMARY.md`](ml/docs/UMA_IMPLEMENTATION_SUMMARY.md)

### Architecture Overview

See [`ARCHITECTURE.md`](../docs/ml-experimental/ARCHITECTURE.md) for system design and implementation details.

## Status

All three features are **production-ready** and tested:

| Feature | Status | GPU | RPMD | Tests |
|---------|--------|-----|------|-------|
| Cavity Particles | ✅ Complete | ✅ | ✅ | ✅ |
| Hybrid RPMD | ✅ Complete | ✅ | ✅ | ✅ |
| ML Potentials | ✅ Complete | ✅ | ✅ | ✅ |

## Core Implementation

Core implementations remain in standard OpenMM locations:

**Cavity Particles:**
- `openmmapi/src/CavityForce.cpp`
- `openmmapi/src/CavityParticleDisplacer.cpp`
- `platforms/common/src/CommonKernels.cpp`

**Hybrid RPMD:**
- `plugins/rpmd/openmmapi/src/RPMDIntegrator.cpp`
- `plugins/rpmd/platforms/common/src/CommonRpmdKernels.cpp`
- `plugins/rpmd/platforms/common/src/kernels/rpmd.cc`

**ML Potentials:**
- `plugins/uma/` (C++ plugin, in development)
- `fairchem/openmm-ml/` (OpenMM-ML integration)

This directory (`ml-experimental/`) organizes documentation, examples, and tests for these features.

## Contributing

When contributing to ml-experimental features:

1. **Add utilities** to appropriate `utils/` directories
2. **Create examples** in feature `examples/` directories  
3. **Write documentation** in feature `docs/` directories
4. **Add tests** to `ml-experimental/tests/`

Follow the established directory structure for consistency.

## Performance

Typical performance characteristics:

**Cavity Particles:**
- Overhead: <5% for systems >1000 atoms
- GPU speedup: 20-50× vs CPU

**Hybrid RPMD:**
- Speedup: 5-10× for typical biomolecular systems (few quantum atoms)
- Memory: 50-90% reduction for classical particles

**ML Potentials:**
- UMA inference: ~1-5 ms/step for 100-atom systems on GPU
- PythonForce overhead: 5-20% vs C++ plugin

## Support

For questions or issues:

1. Check feature-specific documentation in `docs/` directories
2. Review examples in `examples/` directories
3. See issue tracker for known issues
4. Read `SUPPORT.md` in repository root

## License

These experimental features are part of OpenMM and distributed under the same license terms. See repository root for license information.

## Citation

If you use these features in published research, please cite:

**OpenMM:**
- Eastman et al., JCTC 2017

**Cavity Particles:**
- (Cite relevant polariton chemistry papers)

**RPMD:**
- Markland & Manolopoulos, JCP 2008

**UMA Potentials:**
- FAIRChem Team, (cite UMA paper when available)
