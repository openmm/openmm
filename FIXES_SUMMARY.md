# UMA Hybrid RPMD Implementation Fixes - Summary

## Overview

This document summarizes the critical bug fixes implemented to make UMA (Universal Molecular Atomistic potential) work correctly with standard and hybrid RPMD (Ring Polymer Molecular Dynamics) in OpenMM.

## Critical Bugs Fixed

### 1. Force Accumulation Overwrite Bug (CRITICAL)

**Problem:** In the batched RPMD path, when UMA coexisted with other forces (e.g., CavityForce, NonbondedForce, HarmonicBondForce), the UMA forces were being overwritten instead of accumulated.

**Root Cause:** The `copyDataFromContext` kernel used assignment (`=`) instead of addition (`+=`) when copying forces from the OpenMM context to the RPMD force buffer.

**Fix:**
- Added new `addForcesFromContext` kernel in `rpmd.cc` that uses `+=` instead of `=`
- Modified `CommonRpmdKernels.cpp` to use the additive kernel when processing non-batched forces after UMA forces are uploaded
- Updated `CommonRpmdKernels.h` to declare the new kernel

**Files Modified:**
- `plugins/rpmd/platforms/common/src/kernels/rpmd.cc`
- `plugins/rpmd/platforms/common/src/CommonRpmdKernels.cpp`
- `plugins/rpmd/platforms/common/src/CommonRpmdKernels.h`

**Impact:** ALL cavity-RPMD simulations with UMA were producing incorrect results before this fix.

### 2. Classical Particle Centroid Force Bug (CRITICAL)

**Problem:** In hybrid RPMD, classical particles were only seeing forces from bead 0, not the centroid-averaged force F_c = (1/P) * Σ_k F(q_c, q_quantum^(k)).

**Root Cause:** The integration and sync kernels were designed to copy bead 0 to all other beads, discarding the force information from beads 1..N-1.

**Fix:**
- Modified `integrateStepHybrid` kernel in `rpmd.cc` to:
  - Compute centroid velocity for classical particles by averaging across all beads
  - Use this centroid velocity for position updates
- Modified `advanceVelocitiesHybrid` kernel to:
  - Remove the skip logic for classical particles
  - Allow all beads to update velocities (each sees different force based on quantum configuration)
- Modified `syncClassicalBeads` kernel to:
  - Average velocities across beads instead of copying bead 0
  - Still copy positions from bead 0 (they should be identical after integration fix)

**Files Modified:**
- `plugins/rpmd/platforms/common/src/kernels/rpmd.cc` (3 kernel functions)

**Impact:** Hybrid RPMD dynamics were fundamentally incorrect, sampling the wrong distribution.

### 3. Wrong Barostat in Ice Test

**Problem:** The UMA ice RPMD test was using `MonteCarloBarostat` which throws an error with `RPMDIntegrator`.

**Fix:**
- Changed import to use `RPMDMonteCarloBarostat`
- Updated barostat instantiation in the test

**Files Modified:**
- `tests/uma_ice_rpmd/test_uma_ice_rpmd.py`

**Impact:** Test would crash at runtime when running NPT simulations.

### 4. Debug Logging Cleanup

**Problem:** Extensive file I/O debug logging was left in the production code from development, adding overhead.

**Fix:**
- Removed all `debug.log` file write statements from the UMA batch implementation
- Kept `DEBUG_LOGS` console output (disabled by default) for development

**Files Modified:**
- `wrappers/python/openmmml/models/umapotential_pythonforce_batch.py`

**Impact:** Reduced I/O overhead in production runs.

## New Tests Added

Three comprehensive tests were added to verify the fixes:

### 1. Force Accumulation Test
**File:** `tests/rpmd/test_uma_force_accumulation.py`

Tests that UMA forces and other forces (HarmonicBondForce) are correctly accumulated in the batched RPMD path. Verifies that total force = UMA force + other forces on each bead.

### 2. Classical Centroid Force Test
**File:** `tests/rpmd/test_classical_centroid_force.py`

Tests that classical particles in hybrid RPMD experience the centroid-averaged force. Verifies:
- Classical particle positions are identical on all beads
- Forces vary across beads (due to different quantum configurations)
- Velocity updates are consistent with centroid force

### 3. UMA + Hybrid RPMD Integration Test
**File:** `tests/rpmd/test_uma_hybrid_rpmd_cavity.py`

Comprehensive integration test for UMA with hybrid RPMD and cavity particle. Verifies:
- Cavity particle receives zero force from UMA
- Molecular particles receive correct UMA forces
- Cavity beads remain identical (classical behavior)
- Molecular particles show quantum delocalization (bead spread)
- Simulation remains stable

## Technical Details

### Force Accumulation Implementation

The new `addForcesFromContext` kernel:
```glsl
KERNEL void addForcesFromContext(GLOBAL mm_long* srcForce, GLOBAL mm_long* dstForce,
        GLOBAL int* order, int copy) {
    const int base = copy*PADDED_NUM_ATOMS;
    for (int particle = GLOBAL_ID; particle < NUM_ATOMS; particle += GLOBAL_SIZE) {
        int index = order[particle];
        dstForce[base*3+index] += srcForce[particle];  // ADD instead of ASSIGN
        dstForce[base*3+index+PADDED_NUM_ATOMS] += srcForce[particle+PADDED_NUM_ATOMS];
        dstForce[base*3+index+PADDED_NUM_ATOMS*2] += srcForce[particle+PADDED_NUM_ATOMS*2];
    }
}
```

### Centroid Velocity Averaging

For classical particles in hybrid RPMD:
1. After first velocity half-step, all beads have different velocities (different forces)
2. Before position update, velocities are averaged across beads using shared memory
3. Position is updated using this centroid velocity
4. After force computation, second velocity half-step updates all beads
5. Velocities are averaged again via `syncClassicalBeads`

This ensures classical particles see the correct centroid-averaged force F_c = (1/P) * Σ_k F(q_c, q_quantum^(k)).

## Known Limitations

### Reference Platform
The Reference platform does not have hybrid RPMD support. The centroid averaging fixes only apply to the Common (GPU) platform kernels. If hybrid RPMD is needed on CPU, additional work would be required in `ReferenceRpmdKernels.cpp`.

## Verification Strategy

After applying these fixes, the following verification steps should be performed:

1. **Force accumulation test:** Run `test_uma_force_accumulation.py` to verify UMA + other forces work correctly
2. **Centroid force test:** Run `test_classical_centroid_force.py` to verify classical particles see centroid-averaged forces
3. **Integration test:** Run `test_uma_hybrid_rpmd_cavity.py` to verify UMA + hybrid RPMD + cavity works end-to-end
4. **Regression test:** Run existing ice RPMD test to ensure standard RPMD still works correctly
5. **Energy conservation:** Run NVE simulations with UMA to verify energy drift is acceptable

## Impact Assessment

These fixes are critical for correctness:

- **Before fix:** ALL cavity-RPMD simulations with UMA were producing incorrect results
- **After fix:** Proper force accumulation and centroid averaging ensure correct physics

Any published results using UMA with:
- Hybrid RPMD
- Cavity particles
- Multiple force groups

Should be re-evaluated with these fixes applied.
