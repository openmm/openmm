# Hybrid Classical/Quantum RPMD Implementation

## Overview

This document describes the implementation of hybrid classical/quantum simulations in RPMDIntegrator, where specific particles can receive quantum (RPMD with multiple beads) treatment while others are treated classically (single copy).

## Scientific Motivation

Quantum nuclear effects are most important for:
- **Light atoms**: Hydrogen (H), deuterium (D) have significant zero-point energy and tunneling
- **Low temperatures**: Below ~200 K, quantum effects become important
- **High-frequency modes**: O-H, N-H stretches (>3000 cm⁻¹)

For large biomolecular systems, applying RPMD to all atoms is computationally expensive. Hybrid classical/quantum RPMD allows:
- **Quantum treatment** for protons and active site atoms
- **Classical treatment** for protein backbone and solvent
- **6-10× speedup** for typical systems while preserving quantum accuracy where it matters

## API Design

### Type-Based Particle Selection

Following the QTBIntegrator/DPDIntegrator pattern:

```cpp
// C++ API
RPMDIntegrator integrator(8, 300.0, 1.0, 0.001);  // 8 beads, 300K, 1 ps⁻¹, 1 fs

// Assign particle types
integrator.setParticleType(0, 1);   // Particle 0 is type 1
integrator.setParticleType(1, 1);   // Particle 1 is type 1
integrator.setParticleType(2, 0);   // Particle 2 is type 0 (default)

// Set which types are quantum
std::set<int> quantumTypes = {1};
integrator.setQuantumParticleTypes(quantumTypes);

// Set default behavior for type 0
integrator.setDefaultQuantum(false);  // Type 0 is classical

// Configure classical thermostat
integrator.setClassicalThermostat(RPMDIntegrator::BussiClassical);
```

```python
# Python API
integrator = openmm.RPMDIntegrator(8, 300*unit.kelvin, 1.0/unit.picosecond, 
                                    0.001*unit.picosecond)

# Mark hydrogen atoms as quantum
for i, atom in enumerate(topology.atoms()):
    if atom.element.symbol == 'H':
        integrator.setParticleType(i, 1)
    else:
        integrator.setParticleType(i, 0)

# Set type 1 as quantum
integrator.setQuantumParticleTypes({1})
integrator.setDefaultQuantum(False)  # Type 0 = classical

# Configure classical thermostat
integrator.setClassicalThermostat(openmm.RPMDIntegrator.BussiClassical)
```

### ClassicalThermostatType Enum

Three options for classical particle thermostats:

1. **BussiClassical (default)**: Bussi stochastic velocity rescaling
   - Maintains canonical ensemble efficiently
   - No friction parameter needed
   - Best for NVT sampling

2. **LangevinClassical**: Langevin dynamics
   - Uses same friction coefficient as RPMD particles
   - Consistent friction across quantum/classical sectors
   - Good for dynamics

3. **NoneClassical**: No thermostat (NVE)
   - Microcanonical ensemble for classical particles
   - Useful for testing energy conservation

## Memory Layout

### Current (All-Quantum) Layout

```
positions[numParticles × numBeads × 4]
```

For 500 atoms, 8 beads:
- Storage: 500 × 8 × 4 = 16,000 floats
- Memory: 64 KB (single precision) or 128 KB (double precision)

### Hybrid (Sparse) Layout

```
positions[
    0 ... numQuantumParticles × numBeads × 4 - 1:       quantum particles (all beads)
    numQuantumParticles × numBeads × 4 ... end:         classical particles (bead 0 only)
]
```

For 500 atoms (20 H quantum, 480 heavy classical), 8 beads:
- Quantum storage: 20 × 8 × 4 = 640 floats
- Classical storage: 480 × 1 × 4 = 1,920 floats
- Total: 2,560 floats vs 16,000 (6.25× reduction!)
- Memory: 10 KB vs 64 KB (single precision)

### Index Mapping

GPU arrays store mappings:
- `quantumIndices[numQuantumParticles]`: Maps quantum particle local index → global particle index
- `classicalIndices[numClassicalParticles]`: Maps classical particle local index → global particle index

## Integration Algorithm

### Quantum Particles (RPMD)

1. **Thermostat (first half-step)**:
   - PILE_G mode: Apply Bussi to centroid, PILE to internal modes
   - PILE mode: Apply PILE to all modes

2. **Ring polymer evolution**:
   - FFT to normal modes
   - Free ring polymer propagation in normal mode space
   - Inverse FFT to bead space

3. **Force computation**:
   - Compute forces on all beads

4. **Velocity advance**:
   - Half-step velocity update from forces

5. **Thermostat (second half-step)**:
   - Same as step 1

### Classical Particles (Standard MD)

1. **Thermostat (first half-step)**:
   - Bussi: Velocity rescaling based on kinetic energy
   - Langevin: Random force + friction
   - None: Skip

2. **Velocity Verlet**:
   ```
   v(t+dt/2) = v(t) + f(t) × dt/(2m)
   r(t+dt) = r(t) + v(t+dt/2) × dt
   ```

3. **Force computation**:
   - Compute forces on single copy (bead 0)

4. **Velocity advance**:
   ```
   v(t+dt) = v(t+dt/2) + f(t+dt) × dt/(2m)
   ```

5. **Thermostat (second half-step)**:
   - Same as step 1

### Force Computation Strategy

**Challenge**: Forces involve interactions between quantum and classical particles.

**Solution**:
- For each quantum bead `i`:
  - Set quantum particle positions to bead `i` positions
  - Set classical particle positions to bead 0 (always)
  - Compute forces
- Classical particles feel **averaged force** from all quantum beads
- Quantum particles feel force from **static** classical particle positions

This ensures:
- Quantum particles sample quantum distribution
- Classical particles act as mean-field bath
- Energy is properly conserved (total system energy)

## Kernel Modifications

### Phase 3: Kernel Modifications (rpmd.cc)

#### 3a. Modify `applyPileThermostat`

**Current signature**:
```cpp
KERNEL void applyPileThermostat(GLOBAL mixed4* velm, GLOBAL float4* random, 
                                unsigned int randomIndex, mixed dt, mixed kT, 
                                mixed friction, int applyToCentroid)
```

**New signature**:
```cpp
KERNEL void applyPileThermostat(GLOBAL mixed4* velm, GLOBAL float4* random, 
                                GLOBAL const int* restrict quantumIndices,
                                int numQuantumParticles,
                                unsigned int randomIndex, mixed dt, mixed kT, 
                                mixed friction, int applyToCentroid)
```

**Key changes**:
- Add `quantumIndices` array parameter
- Add `numQuantumParticles` parameter
- Change loop to iterate over `numQuantumParticles` instead of `NUM_ATOMS`
- Map particle index: `int globalParticle = quantumIndices[particle]`
- Access velocities using sparse layout

#### 3b. Modify `integrateStep`

Similar modifications for quantum particle indexing.

#### 3c. Add `integrateClassical` Kernel

```cpp
KERNEL void integrateClassical(
    GLOBAL mixed4* restrict posq,
    GLOBAL mixed4* restrict velm,
    GLOBAL const mm_long* restrict force,
    GLOBAL const int* restrict classicalIndices,
    mixed dt
) {
    int particle = GLOBAL_ID;
    if (particle >= NUM_CLASSICAL_PARTICLES) return;
    
    int idx = classicalIndices[particle];
    mixed forceScale = 1/(mixed) 0x100000000;
    
    // Load state
    mixed4 velocity = velm[idx];
    mixed4 position = posq[idx];
    mixed invMass = velocity.w;
    
    if (invMass == 0.0f) return;
    
    // Velocity half-step: v(t+dt/2) = v(t) + f(t)*dt/(2m)
    velocity.x += force[3*idx] * (dt * 0.5f) * forceScale * invMass;
    velocity.y += force[3*idx+1] * (dt * 0.5f) * forceScale * invMass;
    velocity.z += force[3*idx+2] * (dt * 0.5f) * forceScale * invMass;
    
    // Position update: r(t+dt) = r(t) + v(t+dt/2)*dt
    position.x += velocity.x * dt;
    position.y += velocity.y * dt;
    position.z += velocity.z * dt;
    
    // Store updated state
    posq[idx] = position;
    velm[idx] = velocity;
}
```

#### 3d. Add `advanceClassicalVelocities` Kernel

```cpp
KERNEL void advanceClassicalVelocities(
    GLOBAL mixed4* restrict velm,
    GLOBAL const mm_long* restrict force,
    GLOBAL const int* restrict classicalIndices,
    mixed dt
) {
    int particle = GLOBAL_ID;
    if (particle >= NUM_CLASSICAL_PARTICLES) return;
    
    int idx = classicalIndices[particle];
    mixed forceScale = 1/(mixed) 0x100000000;
    
    mixed4 velocity = velm[idx];
    mixed invMass = velocity.w;
    
    if (invMass == 0.0f) return;
    
    // Velocity advance: v(t+dt) = v(t+dt/2) + f(t+dt)*dt/(2m)
    velocity.x += force[3*idx] * (dt * 0.5f) * forceScale * invMass;
    velocity.y += force[3*idx+1] * (dt * 0.5f) * forceScale * invMass;
    velocity.z += force[3*idx+2] * (dt * 0.5f) * forceScale * invMass;
    
    velm[idx] = velocity;
}
```

#### 3e. Add `applyBussiClassical` Kernel

```cpp
KERNEL void applyBussiClassical(
    GLOBAL mixed4* restrict velm,
    GLOBAL const int* restrict classicalIndices,
    mixed scalingFactor
) {
    int particle = GLOBAL_ID;
    if (particle >= NUM_CLASSICAL_PARTICLES) return;
    
    int idx = classicalIndices[particle];
    mixed4 velocity = velm[idx];
    
    if (velocity.w != 0.0f) {
        velocity.x *= scalingFactor;
        velocity.y *= scalingFactor;
        velocity.z *= scalingFactor;
        velm[idx] = velocity;
    }
}
```

#### 3f. Add `applyLangevinClassical` Kernel

```cpp
KERNEL void applyLangevinClassical(
    GLOBAL mixed4* restrict velm,
    GLOBAL const float4* restrict random,
    GLOBAL const int* restrict classicalIndices,
    unsigned int randomIndex,
    mixed dt,
    mixed kT,
    mixed friction
) {
    int particle = GLOBAL_ID;
    if (particle >= NUM_CLASSICAL_PARTICLES) return;
    
    int idx = classicalIndices[particle];
    mixed4 velocity = velm[idx];
    mixed invMass = velocity.w;
    
    if (invMass == 0.0f) return;
    
    mixed c1 = exp(-friction * dt);
    mixed c2 = sqrt(1.0f - c1*c1);
    mixed c3 = c2 * sqrt(kT * invMass);
    
    float4 rand = random[randomIndex + particle];
    
    velocity.x = c1 * velocity.x + c3 * rand.x;
    velocity.y = c1 * velocity.y + c3 * rand.y;
    velocity.z = c1 * velocity.z + c3 * rand.z;
    
    velm[idx] = velocity;
}
```

#### 3g. Modify Data Copy Kernels

Update `copyDataToContext` and `copyDataFromContext` to handle sparse storage:

```cpp
KERNEL void copyDataToContext(
    GLOBAL const mixed4* restrict rpmdPosq,
    GLOBAL mixed4* restrict contextPosq,
    GLOBAL const int* restrict particleIndices,
    int numParticlesInSet,
    int copy
) {
    for (int i = GLOBAL_ID; i < numParticlesInSet; i += GLOBAL_SIZE) {
        int globalIdx = particleIndices[i];
        int storageIdx = i + copy * numParticlesInSet;  // Sparse layout
        contextPosq[globalIdx] = rpmdPosq[storageIdx];
    }
}
```

## Implementation Status

### ✅ Phase 1: API Extension
- [x] Add ClassicalThermostatType enum to RPMDIntegrator.h
- [x] Add particle type management methods
- [x] Add quantum particle type selection methods
- [x] Add classical thermostat configuration methods
- [x] Add private member variables
- [x] Implement getter/setter methods in RPMDIntegrator.cpp
- [x] Initialize new members in constructors

### ⚠️ Phase 2: Particle Classification (Partial)
- [x] Add member variables to CommonRpmdKernels.h
- [x] Implement particle classification in initialize()
- [x] Allocate sparse bead storage
- [ ] **TODO**: Update initialization of positions/velocities for sparse layout
- [ ] **TODO**: Upload quantum/classical index arrays to GPU

### ⏳ Phase 3: Kernel Modifications (TODO)
- [ ] Modify applyPileThermostat kernel signature and implementation
- [ ] Modify integrateStep kernel
- [ ] Add integrateClassical kernel
- [ ] Add advanceClassicalVelocities kernel
- [ ] Add applyBussiClassical kernel
- [ ] Add applyLangevinClassical kernel
- [ ] Modify copyDataToContext/FromContext for sparse storage

### ⏳ Phase 4: Integration Loop (TODO)
- [ ] Restructure execute() in CommonRpmdKernels.cpp
- [ ] Separate quantum and classical integration paths
- [ ] Add classical thermostat application
- [ ] Handle kernel launches for different particle sets

### ⏳ Phase 5: Force Computation (TODO)
- [ ] Modify computeForces() for hybrid system
- [ ] Handle quantum beads + classical bead 0 interactions
- [ ] Update copyToContext/FromContext kernel calls

### ⏳ Phases 6-9: Remaining Work (TODO)
- [ ] Reference platform implementation
- [ ] Serialization support
- [ ] Comprehensive test suite
- [ ] Complete documentation

## Test Plan

### Unit Tests

1. **testHybridClassicalQuantum()**
   - System: 10 atoms (5 H quantum, 5 C classical), 8 beads
   - Verify: Quantum particles have 8 copies, classical have 1
   - Check: Energy conservation (NVE mode)

2. **testClassicalThermostatTypes()**
   - Test: Bussi, Langevin, None thermostats
   - Verify: Temperature distributions match theoretical

3. **testAllQuantum()**
   - Backward compatibility: All particles quantum
   - Should match existing RPMD behavior exactly

4. **testAllClassical()**
   - Edge case: All particles classical
   - Should reduce to standard MD

5. **testMemoryEfficiency()**
   - Large system: 1000 atoms (10 quantum, 990 classical)
   - Verify: Memory usage proportional to (10×8 + 990×1)

### Integration Tests

1. **testHydrogenOnlyQuantum()**
   - Water box: 32 H2O (64 H quantum, 32 O classical)
   - Compare IR spectrum with all-quantum RPMD
   - Verify: O-H stretch peak identical

2. **testProtonTransfer()**
   - Proton transfer reaction in enzyme active site
   - Quantum: H+ and nearby residues
   - Classical: Protein backbone and solvent
   - Verify: Reaction rate matches all-quantum

## Performance Expectations

### Benchmark System
- 250 dimers (500 atoms)
- 20 H atoms → quantum (8 beads)
- 480 heavy atoms → classical (1 copy)

### Expected Performance

| Metric | All-Quantum | Hybrid | Speedup |
|--------|-------------|--------|---------|
| Memory | 128 KB | 20.5 KB | **6.2×** |
| Force calls/step | 8 | 8 | 1× |
| Integration ops | 500×8 | 20×8 + 480×1 | **5.9×** |
| Total speedup | — | — | **~5×** |

### Scaling

For system with `N` total atoms, `Q` quantum atoms, `B` beads:

- **Memory**: `O(Q×B + (N-Q))` vs `O(N×B)`
- **Computation**: `O(Q×B + N-Q)` for integration, `O(N×B)` for forces
- **Speedup**: Approximately `N×B / (Q×B + N-Q)` for integration

## Best Practices

### Choosing Quantum Particles

1. **Light atoms** (H, D): Always quantum below 300 K
2. **High-frequency modes**: Bonds to H (O-H, N-H, C-H)
3. **Active sites**: Atoms directly involved in chemistry
4. **Convergence testing**: Start with more quantum, reduce gradually

### Performance vs Accuracy Trade-offs

- **Conservative**: All H atoms quantum → ~5× slower than all-classical
- **Balanced**: Active site + nearby H quantum → ~10× speedup vs all-quantum
- **Aggressive**: Only transferring H quantum → ~20× speedup (check convergence!)

### Classical Thermostat Choice

- **Bussi**: Best for equilibrium sampling (NVT), canonical ensemble
- **Langevin**: Better for dynamics, friction matches quantum sector
- **None**: Use only for testing energy conservation

## References

1. Habershon et al., "Ring-polymer molecular dynamics: Quantum effects in chemical dynamics from classical trajectories", *Annu. Rev. Phys. Chem.* **64**, 387 (2013)

2. Ceriotti et al., "Efficient stochastic thermostatting of path integral molecular dynamics", *J. Chem. Phys.* **133**, 124104 (2010)

3. Bussi et al., "Canonical sampling through velocity rescaling", *J. Chem. Phys.* **126**, 014101 (2007)

4. Markland & Manolopoulos, "An efficient ring polymer contraction scheme for imaginary time path integral simulations", *J. Chem. Phys.* **129**, 024105 (2008)

5. Cheng & Ceriotti, "Direct path integral estimators for isotope fractionation ratios", *J. Chem. Phys.* **141**, 244112 (2014)

## Example Usage

### Python Example: Water with Quantum Hydrogens

```python
import openmm
from openmm import app, unit
import sys

# Load water system
pdb = app.PDBFile('water_box.pdb')
forcefield = app.ForceField('tip3p.xml')
system = forcefield.createSystem(pdb.topology, 
                                 nonbondedCutoff=1.0*unit.nanometer,
                                 constraints=None)  # No constraints for RPMD

# Create hybrid RPMD integrator
integrator = openmm.RPMDIntegrator(
    numBeads=8,
    temperature=300*unit.kelvin,
    frictionCoeff=1.0/unit.picosecond,
    stepSize=0.0005*unit.picosecond  # 0.5 fs
)

# Set PILE_G thermostat (Bussi on centroid)
integrator.setThermostatType(openmm.RPMDIntegrator.PileG)

# Mark hydrogen atoms as quantum
for i, atom in enumerate(pdb.topology.atoms()):
    if atom.element.symbol == 'H':
        integrator.setParticleType(i, 1)  # Type 1 = H
    else:
        integrator.setParticleType(i, 0)  # Type 0 = O

# Set type 1 (H) as quantum, type 0 (O) as classical
integrator.setQuantumParticleTypes({1})
integrator.setDefaultQuantum(False)  # Type 0 is classical

# Use Bussi thermostat for classical particles
integrator.setClassicalThermostat(openmm.RPMDIntegrator.BussiClassical)

# Create context
platform = openmm.Platform.getPlatformByName('CUDA')
context = openmm.Context(system, integrator, platform)

# Set initial positions
context.setPositions(pdb.positions)

# Initialize velocities
context.setVelocitiesToTemperature(300*unit.kelvin)

# Equilibration
print("Equilibrating...")
integrator.step(10000)  # 5 ps

# Production
print("Production run...")
for i in range(1000):
    integrator.step(1000)  # 0.5 ps per iteration
    if i % 100 == 0:
        state = integrator.getState(0, getEnergy=True)
        print(f"Step {i*1000}: E = {state.getPotentialEnergy()}")

print("Done!")
```

### Performance Monitoring

```python
import time

start = time.time()
integrator.step(10000)
elapsed = time.time() - start

ns_per_day = (10000 * 0.0005e-3) / elapsed * 86400  # ns/day
print(f"Performance: {ns_per_day:.1f} ns/day")
print(f"Time per step: {elapsed/10000*1000:.2f} ms")
```

## Validation Checklist

- [ ] Energy conservation in NVE mode (< 0.1% drift per ns)
- [ ] Temperature distribution matches Boltzmann (quantum and classical)
- [ ] Centroid theorem: Quantum centroids sample classical distribution
- [ ] IR spectrum matches all-quantum RPMD for quantum modes
- [ ] Performance scales linearly with number of quantum particles
- [ ] Memory usage matches theoretical prediction
- [ ] Backward compatible: All-quantum mode identical to current implementation
- [ ] Serialization preserves hybrid simulation state
