# Hybrid Classical/Quantum RPMD Implementation

## Overview

This implementation extends OpenMM's `RPMDIntegrator` to support hybrid simulations where some particles are treated quantum mechanically (full RPMD with ring polymer dynamics) while others are treated classically (standard molecular dynamics).

## Key Design Principles

Following guidance from i-PI's architecture:

1. **Uniform Memory Layout**: All particles stored with all beads for simplicity and robustness
2. **All-GPU Execution**: Hybrid logic implemented entirely on GPU for maximum performance  
3. **Per-Particle Flags**: GPU array `isQuantum[numParticles]` marks each particle as quantum (1) or classical (0)
4. **Minimal Memory Overhead**: Forces keep full layout for OpenMM compatibility; positions/velocities could be sparse in future optimization

## Memory Considerations

Currently uses uniform layout (all particles × all beads). For a 100 O + 200 H system with 16 beads:
- Position storage: 300 × 16 × 12 bytes = 57.6 KB (single precision)
- This is typically small compared to other simulation data (nonbonded lists, etc.)

Future optimization could implement sparse storage for larger systems.

## API

### Particle Type Assignment

```python
from openmm import RPMDIntegrator

integrator = RPMDIntegrator(numCopies=16, temperature=300, friction=1.0, stepSize=0.0005)

# Assign particle types (optional, for grouping)
integrator.setParticleType(0, 1)  # Particle 0 is type 1 (e.g., O)
integrator.setParticleType(1, 2)  # Particle 1 is type 2 (e.g., H)
integrator.setParticleType(2, 2)  # Particle 2 is type 2 (e.g., H)

# Mark which types are quantum
integrator.setQuantumParticleTypes({2})  # Type 2 (H) is quantum

# Set default behavior for untyped particles
integrator.setDefaultQuantum(False)  # Default is classical
```

### Classical Thermostat Options

```python
# Classical particles can use different thermostats:
integrator.setClassicalThermostat(RPMDIntegrator.BussiClassical)   # Bussi velocity rescaling
integrator.setClassicalThermostat(RPMDIntegrator.LangevinClassical)  # Langevin dynamics
integrator.setClassicalThermostat(RPMDIntegrator.NoneClassical)    # No thermostat (NVE)
```

## Implementation Details

### Memory Layout

All particles are stored with all beads:
- `positions[numCopies * paddedNumAtoms]`
- `velocities[numCopies * paddedNumAtoms]`
- `forces[numCopies * paddedNumAtoms * 3]`

### GPU Arrays

- `isQuantum[paddedNumAtoms]`: Per-particle flag (1=quantum, 0=classical)
- `classicalKE[numWorkGroups]`: For Bussi thermostat kinetic energy reduction

### GPU Kernels

#### Core Kernels (modified for hybrid)
- `integrateStep`: Handles both quantum (FFT ring polymer) and classical (velocity Verlet) particles
- `advanceVelocities`: Velocity update for all particles
- `applyPileThermostat`: PILE thermostat (quantum particles only)

#### Hybrid-Specific Kernels
- `syncClassicalBeads`: Copies bead 0 to all other beads for classical particles
- `applyClassicalThermostat`: Langevin thermostat for classical particles
- `computeClassicalKEHybrid`: Kinetic energy reduction for Bussi thermostat
- `applyBussiClassicalHybrid`: Velocity scaling for Bussi thermostat

### Integration Algorithm

1. **Thermostat (First Half)**
   - Quantum particles: PILE/PILE_G thermostat (internal modes + optionally centroid)
   - Classical particles: Langevin or Bussi thermostat (bead 0 only)

2. **Integration**
   - All particles: Standard RPMD `integrateStep` kernel
   - Classical particles: `syncClassicalBeads` to couple all beads

3. **Force Computation**
   - Forces computed on all beads for all particles
   - Classical particles naturally get identical forces since all beads are at same position

4. **Velocity Update**
   - Standard `advanceVelocities` kernel
   - Classical particles: `syncClassicalBeads` again

5. **Thermostat (Second Half)**
   - Same as first half

## Usage Example

```python
import openmm as mm
from openmm import unit, RPMDIntegrator

# Create system (e.g., water)
system = mm.System()
system.addParticle(16.0 * unit.amu)  # O
system.addParticle(1.0 * unit.amu)   # H
system.addParticle(1.0 * unit.amu)   # H

# Add forces (bonds, etc.)
harmonic = mm.HarmonicBondForce()
harmonic.addBond(0, 1, 0.1 * unit.nanometer, 1000 * unit.kJ_per_mole / unit.nm**2)
harmonic.addBond(0, 2, 0.1 * unit.nanometer, 1000 * unit.kJ_per_mole / unit.nm**2)
system.addForce(harmonic)

# Create hybrid RPMD integrator
integrator = RPMDIntegrator(16, 300, 1.0 / unit.picosecond, 0.5 * unit.femtosecond)

# Configure: O is classical, H is quantum
integrator.setParticleType(0, 1)  # O
integrator.setParticleType(1, 2)  # H
integrator.setParticleType(2, 2)  # H
integrator.setQuantumParticleTypes({2})
integrator.setDefaultQuantum(False)
integrator.setClassicalThermostat(RPMDIntegrator.LangevinClassical)

# Create context and run
platform = mm.Platform.getPlatformByName('CUDA')
context = mm.Context(system, integrator, platform)
context.setPositions(positions)
integrator.step(1000)
```

## Performance

Tested on CUDA with 100 water molecules (300 particles), 16 beads:
- ~7000 steps/second
- ~300 ns/day

The hybrid mode has minimal overhead compared to full RPMD since:
1. Memory layout is unchanged
2. GPU kernels check `isQuantum` flag inline
3. Classical bead sync is a simple copy operation

## Files Modified

- `plugins/rpmd/openmmapi/include/openmm/RPMDIntegrator.h`: API additions
- `plugins/rpmd/openmmapi/src/RPMDIntegrator.cpp`: API implementation
- `plugins/rpmd/platforms/common/src/CommonRpmdKernels.h`: Kernel declarations
- `plugins/rpmd/platforms/common/src/CommonRpmdKernels.cpp`: Kernel logic
- `plugins/rpmd/platforms/common/src/kernels/rpmd.cc`: GPU kernels
- `wrappers/python/src/swig_doxygen/swigInputConfig.py`: Python bindings

## Testing

Run the test scripts:
```bash
python3 test_hybrid_rpmd_v2.py     # Basic tests (Reference platform)
python3 test_hybrid_rpmd_gpu.py   # GPU tests (CUDA/OpenCL)
```
