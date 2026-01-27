# Cavity Particle Implementation

## Overview

This directory contains the cavity particle implementation for OpenMM, which enables quantum light-matter coupling in molecular dynamics simulations. The cavity particle represents a quantized electromagnetic field mode that couples to the molecular dipole moment.

## Key Concepts

### Physical Model

The cavity-matter coupling Hamiltonian is:

```
H = H_matter + ℏωc a†a + λ(a + a†) · μ
```

where:
- `H_matter`: Molecular Hamiltonian
- `ωc`: Cavity frequency
- `a†, a`: Photon creation/annihilation operators
- `λ`: Light-matter coupling strength
- `μ`: Molecular dipole moment

### Implementation

The cavity photon mode is represented as a dummy particle with:
- Mass: ~1 a.u. (very light, decoupled from molecular dynamics)
- Position: Represents the cavity field amplitude
- No nonbonded interactions with matter

Two custom forces implement the coupling:

1. **CavityForce**: Dipole-cavity coupling term λ(q_cavity) · μ
2. **CavityParticleDisplacer**: Finite-q corrections for cavity field displacement

## Directory Structure

```
cavity/
├── docs/                           # Cavity-specific documentation
│   ├── CAVITY_IMPLEMENTATION.md   # This file
│   ├── CAVITY_USAGE.md            # Usage examples
│   └── CAVITY_THEORY.md           # Theoretical background
├── examples/                       # Example systems
│   ├── dimer_system/              # Simple toy model
│   ├── water_system/              # Realistic water simulations
│   └── protein_system/            # Protein cavity dynamics
└── utils/                          # Shared utilities
    ├── cavity_helpers.py          # Common cavity functions
    └── __init__.py
```

## Core Implementation Files

The core C++ implementation is located in:
- `openmmapi/src/CavityForce.cpp` - CavityForce implementation
- `openmmapi/src/CavityParticleDisplacer.cpp` - Displacer implementation
- `openmmapi/include/openmm/CavityForce.h` - CavityForce header
- `openmmapi/include/openmm/CavityParticleDisplacer.h` - Displacer header

Platform-specific implementations:
- `platforms/reference/src/ReferenceCavityKernels.cpp` - CPU reference
- `platforms/common/src/CommonKernels.cpp` - GPU kernels (CUDA/OpenCL)

## Usage

### Basic Setup

```python
from ml_experimental.cavity.utils import (
    add_cavity_particle, 
    setup_cavity_coupling,
    wavenumber_to_hartree
)

# Add cavity particle
omegac_au = wavenumber_to_hartree(3663.0)  # OH stretch frequency
cavity_idx = add_cavity_particle(system, positions, omegac_au)

# Setup coupling
cavity_force, displacer = setup_cavity_coupling(
    system, cavity_idx, omegac_au, lambda_coupling=0.01
)
```

### Complete Example

See `examples/water_system/run_simulation.py` for a complete working example.

## Key Features

1. **GPU Acceleration**: Full CUDA/OpenCL support
2. **RPMD Compatible**: Works with ring polymer molecular dynamics
3. **Flexible Integration**: Compatible with standard OpenMM integrators
4. **Efficient**: Minimal overhead (<5% for typical systems)

## Physical Parameters

### Typical Cavity Frequencies

- **Water HOH bending**: ~1600 cm⁻¹ (0.00729 Hartree)
- **Water OH stretch**: ~3663 cm⁻¹ (0.01669 Hartree)
- **Carbonyl stretch**: ~1700 cm⁻¹ (0.00775 Hartree)
- **Amide I**: ~1650 cm⁻¹ (0.00752 Hartree)

### Coupling Strengths

Typical values for observing Rabi splitting:
- **Weak coupling**: λ = 0.001-0.01
- **Strong coupling**: λ = 0.01-0.1
- **Ultra-strong**: λ > 0.1

### Photon Mass

Standard choice: `1.0 / 1822.888` amu = 1 atomic unit of mass

This ensures the photon particle moves much faster than nuclei, effectively averaging over molecular motion.

## Observables

### Rabi Splitting

The signature of strong light-matter coupling is the appearance of two polariton peaks in the vibrational spectrum, split by the Rabi frequency:

```
ΩR ∝ λ√N
```

where N is the number of molecules coupled to the cavity.

### Cavity Displacement

Monitor the cavity particle position to track the field amplitude:

```python
from ml_experimental.cavity.utils import get_cavity_displacement

state = context.getState(getPositions=True)
disp, mag = get_cavity_displacement(state, cavity_idx)
print(f"Cavity displacement: {mag:.4f} nm")
```

## Testing

Unit tests:
- C++: `tests/TestCavityForce.h`
- Platform-specific: `platforms/cuda/tests/TestCudaCavityForce.cpp`

Integration tests:
- Simple: `examples/dimer_system/`
- Realistic: `examples/water_system/`
- Complex: `examples/protein_system/`

## References

1. **Polariton Chemistry**: Ebbesen, T. W. (2016). Acc. Chem. Res. 49, 2403.
2. **Cavity QED**: Flick, J. et al. (2017). Proc. Natl. Acad. Sci. 114, 3026.
3. **Implementation**: See `CAVITY_THEORY.md` for detailed derivations.

## Status

**Implementation**: Complete and tested
- ✅ CPU (Reference platform)
- ✅ GPU (CUDA/OpenCL)
- ✅ RPMD integration
- ✅ Python bindings

**Validated Against**:
- Analytical results for harmonic oscillators
- Literature benchmark systems
- HOOMD-blue cavity implementation

## Contributing

When adding new cavity-related code:
1. Add common utilities to `utils/cavity_helpers.py`
2. Create example systems in `examples/`
3. Add documentation to `docs/`
4. Write tests for new functionality

## Known Limitations

1. **Single cavity mode**: Currently limited to one cavity mode per simulation
2. **Dipole approximation**: Uses long-wavelength approximation (valid for molecular systems)
3. **Classical field**: Cavity treated semiclassically (quantum matter + classical field)

Future work may address these limitations.
