# Critical Fixes Applied to Dimer Simulation

## Date: 2026-01-21

## Issues Identified and Fixed

### 1. Force Constants Were WRONG (√2 Too High)
**Problem:** Frequencies were 1.414× too high due to incorrect force constant conversion.
- **Observed:** O-O at ~2200 cm⁻¹, N-N at ~3300 cm⁻¹  
- **Expected:** O-O at ~1560 cm⁻¹, N-N at ~2325 cm⁻¹

**Root Cause:** I incorrectly assumed HOOMD's harmonic bond potential didn't have a 1/2 factor, so I doubled the force constants.

**Truth:** Both HOOMD and OpenMM use `E = 0.5 * k * (r - r0)²`, so force constants should be IDENTICAL.

**Fix Applied:**
```python
# BEFORE (WRONG - frequencies √2 too high):
k_OO_au = 2 * 0.73204  # = 1.46408 → ω ≈ 2200 cm⁻¹
k_NN_au = 2 * 1.4325   # = 2.865  → ω ≈ 3300 cm⁻¹

# AFTER (CORRECT):
k_OO_au = 0.73204   # → ω ≈ 1560 cm⁻¹ ✓
k_NN_au = 1.4325    # → ω ≈ 2325 cm⁻¹ ✓
```

### 2. Cavity Coupling Was NOT Active from t=0
**Problem:** No Rabi splitting observed because coupling was OFF during equilibration.

**Evidence from terminal output:**
```
Initial lambda: 0.0
Coupling will turn on at step 100000 with lambda=0.001
```

This meant the first 100 ps had λ=0 (no coupling), so molecules evolved independently!

**Fix Applied:**
```python
# BEFORE (WRONG - coupling OFF initially):
cavity_force = openmm.CavityForce(cavity_index, omegac_au, 0.0, photon_mass)
cavity_force.setCouplingOnStep(equilibration_steps, lambda_coupling)

# AFTER (CORRECT - coupling ON from t=0):
cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass)
# No setCouplingOnStep() - coupling active immediately!
```

### 3. Missing Finite-Q Displacement
**Problem:** Cavity particle started at rest with no initial excitation.

**Fix Applied:**
1. Added `CavityParticleDisplacer`:
   ```python
   displacer = openmm.CavityParticleDisplacer(cavity_index, photon_mass)
   displacer.setSwitchOnLambda(lambda_coupling)
   system.addForce(displacer)
   ```

2. Applied displacement after minimization, before dynamics:
   ```python
   # Displace cavity particle by 0.01 nm along x-axis
   positions[cavity_index] = [0.01, 0.0, 0.0] * unit.nanometer
   ```

## Verification Output

From the corrected simulation start:
```
--- Adding Cavity Force ---
  Lambda coupling: 0.001 (ACTIVE from t=0)  ✓
  CavityParticleDisplacer added (finite-Q mode)  ✓

--- Applying Finite-Q Displacement ---
  Cavity position before: (0.005747, -0.000262, 0.000000) nm
  Cavity position after:  (0.010000, 0.000000, 0.000000) nm  ✓
  Displacement: 0.01 nm along x-axis

--- Running Full Simulation (Coupling ON from t=0) ---
  Running 200000 steps (1 ns)...
```

## Expected Results (After Correction)

### Vibrational Frequencies
- **O-O stretch:** ~1560 cm⁻¹ (dominant peak, 80% of dimers)
- **N-N stretch:** ~2325 cm⁻¹ (smaller peak, 20% of dimers)

### Rabi Splitting (with λ=0.001, ωc=1560 cm⁻¹)
The O-O peak should split into:
- **Lower polariton:** ωc - ΔE/2
- **Upper polariton:** ωc + ΔE/2

Where Rabi splitting: ΔE ≈ 2λ√(N_OO) ωc

For 200 O-O dimers:
- ΔE ≈ 2 × 0.001 × √200 × 1560 ≈ 44 cm⁻¹

Expected polariton peaks:
- Lower: ~1538 cm⁻¹
- Upper: ~1582 cm⁻¹

## Files Modified

1. `tests/dimer_system/run_simulation.py`
   - Corrected force constants (k_OO, k_NN)
   - Fixed CavityForce initialization (λ ≠ 0 from t=0)
   - Added CavityParticleDisplacer
   - Added finite-Q displacement before dynamics

2. `tests/dimer_system/README.md`
   - Updated with correct frequencies and parameters

3. `tests/FORCE_CONSTANT_FIX.md`
   - Documented the force constant correction

## Current Status

✅ Force constants corrected  
✅ Coupling active from t=0  
✅ Finite-Q displacement applied  
🔄 Simulation running (ETA: ~8.5 minutes)  
⏳ Awaiting spectrum analysis to confirm Rabi splitting

## Next Steps

1. Wait for simulation to complete (~200 ps at 1 fs timestep)
2. Run: `python analyze_spectrum.py cavity_diamer_lambda0.0010.npz`
3. Verify:
   - Peaks at correct frequencies (1560, 2325 cm⁻¹)
   - O-O peak shows Rabi splitting (~44 cm⁻¹)
   - N-N peak remains sharp (off-resonant)

## References

- cav-hoomd source: `src/cavitymd/simulation/core.py` (lines 2204-2205)
- cav-hoomd docs: `docs/part2_theory/introduction.rst` (line 349: "O₂ - 1580 cm⁻¹")
- HOOMD potential: E = 0.5 * k * (r - r0)²
- OpenMM potential: E = 0.5 * k * (r - r0)² (IDENTICAL!)
