# GPU-Optimized Bussi Thermostat Implementation for RPMD

## Overview

This document describes the GPU-optimized implementation of the Bussi (Stochastic Velocity Rescaling) thermostat for the centroid mode in Ring Polymer Molecular Dynamics (RPMD), following the i-PI PILE_G convention.

## Architecture

The implementation keeps all velocity data on the GPU throughout the thermostatting process, minimizing CPU-GPU data transfers. Only the reduced centroid kinetic energy value (a single scalar) is transferred to the CPU for computing the scaling factor.

### Key Components

#### 1. GPU Kernels (`rpmd.cc`)

**`computeCentroidKE`**: Computes the total kinetic energy of the centroid mode
- Each thread processes one particle across all beads
- Computes centroid velocity: `v_centroid = (1/N) * Σ v_bead[i]`
- Calculates kinetic energy: `KE = 0.5 * m * v_centroid²`
- Performs parallel reduction within work groups
- Output: per-work-group kinetic energy contributions

**`applyBussiScaling`**: Applies the Bussi scaling factor to velocities
- Takes the computed scaling factor `alpha` as input
- For each particle, computes centroid velocity on GPU
- Applies scaling: `v_new[bead] = v_old[bead] + (alpha - 1) * v_centroid`
- This preserves internal mode velocities while scaling only the centroid component

#### 2. Host-Side Orchestration (`CommonRpmdKernels.cpp`)

**`applyBussiCentroidThermostat`**: Coordinates the thermostatting process
1. **GPU Phase 1**: Launch `computeCentroidKEKernel` to compute centroid KE
2. **CPU Phase**: 
   - Download reduced KE values (one per work group)
   - Sum to get total centroid kinetic energy
   - Generate Gaussian random numbers for Bussi algorithm
   - Compute scaling factor `alpha` using Bussi formula
3. **GPU Phase 2**: Launch `applyBussiScalingKernel` with computed `alpha`

### Bussi Algorithm

The scaling factor is computed using:

```
ratio = K_target / K_current
alpha² = c1 + ratio*(1-c1)*(R1² + R_gamma) + 2*R1*sqrt(ratio*c1*(1-c1))
alpha = sqrt(alpha²)
```

Where:
- `K_target = 0.5 * ndof * kT` (target kinetic energy)
- `c1 = exp(-γ * dt)` (friction damping factor)
- `R1` ~ N(0,1) (Gaussian random number)
- `R_gamma` = Σ[i=1 to ndof-1] Ri² (chi-squared distribution)

## Performance Benefits

### Compared to CPU-Based Implementation

| Operation | CPU-Based | GPU-Optimized |
|-----------|-----------|---------------|
| Data Transfer | Full velocity array | Single scalar (KE) |
| Centroid Calculation | O(N*M) serial | O(N*M) parallel |
| Scaling Application | O(N*M) serial | O(N*M) parallel |

Where N = number of particles, M = number of beads

### Transfer Volume Reduction

For a typical system:
- N = 10,000 particles
- M = 8 beads
- Precision = double (8 bytes)

**CPU-based**: 
- Download: 10,000 × 8 × 3 × 8 = 1.92 MB
- Upload: 1.92 MB
- **Total: 3.84 MB per thermostat application**

**GPU-optimized**:
- Download: ~100 work groups × 8 bytes = 800 bytes
- Upload: 0 bytes
- **Total: 800 bytes per thermostat application**

**Improvement: ~5000× reduction in data transfer!**

### Computational Speedup

The GPU kernels parallelize both the centroid velocity computation and the scaling application, providing significant speedup on modern GPUs with thousands of CUDA cores.

## Integration with RPMD

The Bussi thermostat is applied only in `PILE_G` mode:

```cpp
if (applyThermostat && thermostatType == RPMDIntegrator::PileG) {
    applyBussiCentroidThermostat(system, integrator, dt*0.5);
}
pileKernel->execute(numParticles*numCopies, workgroupSize);  // Internal modes
```

This ensures:
1. Bussi thermostat acts on centroid mode only
2. PILE Langevin thermostat acts on internal modes
3. The `applyToCentroid` flag in the PILE kernel prevents double-thermostatting

## Usage from Python

```python
import openmm

# Create RPMD integrator with PILE_G thermostat
integrator = openmm.RPMDIntegrator(
    num_beads=8,
    temperature=300.0 * unit.kelvin,
    friction=1.0 / unit.picosecond,
    dt=0.5 * unit.femtosecond
)

# Enable PILE_G mode (Bussi on centroid + PILE on internal modes)
integrator.setThermostatType(openmm.RPMDIntegrator.PileG)

# Optional: set different friction for centroid
integrator.setCentroidFriction(0.5)  # ps^-1

# Create context (automatically uses GPU if available)
context = openmm.Context(system, integrator, platform)
```

## References

- Bussi, G., Donadio, D., & Parrinello, M. (2007). Canonical sampling through velocity rescaling. *J. Chem. Phys.*, 126, 014101.
- Ceriotti, M., et al. (2010). Efficient stochastic thermostatting of path integral molecular dynamics. *J. Chem. Phys.*, 133, 124104.
- i-PI documentation: https://ipi-code.org/

## Technical Notes

- The parallel reduction in `computeCentroidKE` uses shared memory for efficiency
- Work group size is set to match the number of beads for optimal FFT performance
- The implementation handles both single and double precision automatically
- Random number generation for the Bussi algorithm uses OpenMM's built-in RNG
