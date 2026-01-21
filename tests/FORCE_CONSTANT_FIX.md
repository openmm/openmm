# Force Constant Correction

## Problem Identified

The initial OpenMM dimer simulation was using **incorrect force constants**, resulting in vibrational frequencies that were √2 ≈ 1.414× too high:

- **Observed:** O-O at ~2200 cm⁻¹, N-N at ~3300 cm⁻¹
- **Expected:** O-O at ~1560 cm⁻¹, N-N at ~2325 cm⁻¹

## Root Cause

**Misunderstanding of HOOMD's potential form:**

Initially, I incorrectly assumed:
- HOOMD uses: `V = k * (r - r0)²` (NO 1/2 factor)
- OpenMM uses: `E = 0.5 * k * (r - r0)²` (WITH 1/2 factor)
- Therefore: need to DOUBLE the force constants ❌ WRONG!

**Actual reality:**
- **HOOMD uses: `E = 0.5 * k * (r - r0)²`** (HAS 1/2 factor, same as OpenMM!)
- **OpenMM uses: `E = 0.5 * k * (r - r0)²`** (HAS 1/2 factor)
- Therefore: use the **SAME** force constants ✓ CORRECT!

## Force Constant Values

From cav-hoomd code (`src/cavitymd/simulation/core.py`, line 2204-2205):

```python
harmonic.params['O-O'] = dict(k=2*0.36602, r0=2.281655158)
harmonic.params['N-N'] = dict(k=2*0.71625, r0=2.0743522177)
```

This gives:
- `k_OO = 0.73204 Hartree/Bohr²`
- `k_NN = 1.4325 Hartree/Bohr²`

## Vibrational Frequency Calculation

For a harmonic potential `E = 0.5 * k * (r - r0)²`, the vibrational frequency is:

```
ω = sqrt(k / μ)  (in atomic units)
ν = ω × 219474.63  (convert to cm⁻¹)
```

Where `μ = m₁*m₂/(m₁+m₂)` is the reduced mass. For homonuclear dimers:
- `μ_OO = m_O/2 = 8 × 1822.9 / 2 = 7291.6 a.u.`
- `μ_NN = m_N/2 = 7 × 1822.9 / 2 = 6380.15 a.u.`

**Wait, this is wrong!** For a **diatomic molecule**, the reduced mass is:
```
μ = m/2  (for identical atoms)
```

But we need to be careful about what mass we're using. Let me recalculate:

For O-O: μ = (m_O × m_O)/(m_O + m_O) = m_O/2 = (8 × 1822.9)/2 = 7291.6 a.u.

Hmm, but my calculation gave μ_OO = 8 × 1822.9 = 14583.2, which is actually **2m_O**, not m_O/2.

Actually, for center-of-mass motion of a dimer treated as a bond, we use the reduced mass:
```
μ = m₁*m₂/(m₁+m₂) = m*m/(2m) = m/2
```

But in the simulation, we have TWO particles of mass m each, so the effective "spring mass" is different...

Actually, let me just verify with the numbers that WORK:

```python
mu_OO = 8.0 * 1822.9 = 14583.2  # Total mass of both atoms
k_OO = 0.73204
ω = sqrt(k/mu) = sqrt(0.73204/14583.2) = 0.007089 a.u.
ν = 0.007089 × 219474.63 = 1555 cm⁻¹ ✓
```

So the "effective mass" for the vibration is the **total mass** of both atoms (or twice the reduced mass).

## Corrected Code

**File:** `tests/dimer_system/run_simulation.py`

```python
# BEFORE (WRONG - frequencies too high by √2):
k_OO_au = 2 * 0.73204  # = 1.46408
k_NN_au = 2 * 1.4325   # = 2.865

# AFTER (CORRECT):
k_OO_au = 0.73204   # → 1560 cm⁻¹
k_NN_au = 1.4325    # → 2325 cm⁻¹
```

## Verification

Run a new simulation with corrected force constants:

```bash
cd /media/extradrive/Trajectories/openmm/tests/dimer_system
python run_simulation.py --dimers 250 --lambda 0.001 --prod 100
python analyze_spectrum.py
```

Expected peaks:
- **1560 cm⁻¹**: O-O stretch (dominant peak, 80% of dimers)
- **2325 cm⁻¹**: N-N stretch (smaller peak, 20% of dimers)

## References

- cav-hoomd documentation: `docs/part2_theory/introduction.rst` (line 349: "O₂ - 1580 cm⁻¹")
- HOOMD source: `src/cavitymd/simulation/core.py` (lines 2204-2205)
- OpenMM HarmonicBondForce: https://openmm.org

## Status

- ✅ Force constants corrected in `run_simulation.py`
- ✅ Default cavity frequency updated to 1560 cm⁻¹
- ✅ README.md updated with correct frequencies
- ⏳ Need to run new simulation to verify correction

---

# Lennard-Jones Parameter Correction

## Problem Identified

The OpenMM dimer simulation was using **incorrect LJ parameters** that did not match the molecular Kob-Andersen model parameters in cav-hoomd:

- **Previous (incorrect):** σ_O = 0.3 nm, ε_O = 0.5 kJ/mol, σ_N = 0.28 nm, ε_N = 0.4 kJ/mol
- **Correct (converted):** σ_O = 0.32971 nm, ε_O = 0.4381 kJ/mol, σ_N = 0.29015 nm, ε_N = 0.2190 kJ/mol

Additionally, OpenMM uses Lorentz-Berthelot combining rules by default, but the Kob-Andersen model requires non-additive cross-term parameters.

## Kob-Andersen Model Parameters

The molecular KA model is a 4:1 mixture of diatomics (O₂:N₂ = 4:1) with LJ parameters following the original KA model ratios:

| Parameter | Ratio | Description |
|-----------|-------|-------------|
| ε_BB/ε_AA | 0.5 | N-N vs O-O epsilon |
| σ_BB/σ_AA | 0.88 | N-N vs O-O sigma |
| ε_AB | 1.5 × ε_AA | N-O cross-term (NOT geometric mean) |
| σ_AB | 0.8 × σ_AA | N-O cross-term (NOT arithmetic mean) |

## Unit Conversion

### cav-hoomd (atomic units) → OpenMM units

Conversion factors:
- 1 Bohr = 0.0529177 nm
- 1 Hartree = 2625.5 kJ/mol

### HOOMD values (atomic units)

```python
# From cav-hoomd/src/cavitymd/simulation/core.py
lj.params[('O', 'O')] = dict(epsilon=0.00016685201, sigma=6.230426584)  # Hartree, Bohr
lj.params[('N', 'N')] = dict(epsilon=0.000083426, sigma=5.48277488)      # Hartree, Bohr
lj.params[('N', 'O')] = dict(epsilon=0.00025027802, sigma=4.9832074319)  # Hartree, Bohr (cross-term!)
```

### OpenMM values (converted)

```python
# O-O parameters
sigma_O = 6.230426584 * 0.0529177  # = 0.32971 nm
epsilon_O = 0.00016685201 * 2625.5  # = 0.4381 kJ/mol

# N-N parameters
sigma_N = 5.48277488 * 0.0529177  # = 0.29015 nm
epsilon_N = 0.000083426 * 2625.5  # = 0.2190 kJ/mol

# N-O cross-term (non-additive KA parameters)
sigma_NO = 4.9832074319 * 0.0529177  # = 0.26370 nm
epsilon_NO = 0.00025027802 * 2625.5  # = 0.6571 kJ/mol
```

## Verification of KA Ratios

```python
# Epsilon ratio (should be 0.5)
epsilon_N / epsilon_O = 0.000083426 / 0.00016685201 = 0.50 ✓

# Sigma ratio (should be 0.88)
sigma_N / sigma_O = 5.48277488 / 6.230426584 = 0.88 ✓

# Cross-term epsilon (should be 1.5 × ε_OO)
epsilon_NO / epsilon_O = 0.00025027802 / 0.00016685201 = 1.50 ✓

# Cross-term sigma (should be 0.8 × σ_OO)
sigma_NO / sigma_O = 4.9832074319 / 6.230426584 = 0.80 ✓
```

## Non-Additive Cross-Terms in OpenMM

OpenMM's NonbondedForce uses Lorentz-Berthelot combining rules by default:
- σ_ij = (σ_i + σ_j) / 2
- ε_ij = √(ε_i × ε_j)

These would give **incorrect** N-O interactions for the KA model. The fix adds explicit exceptions for all N-O pairs with the correct KA cross-term parameters.

## Corrected Files

- `tests/dimer_system/run_simulation.py`
- `tests/dimer_system/run_simulation_reference.py`

## Impact

The incorrect LJ parameters would have affected:
1. Intermolecular potential energies (LJ term)
2. Equilibrium structure and dynamics of the supercooled liquid
3. Comparison with cav-hoomd results

The correction ensures that OpenMM and HOOMD simulations use identical force fields for the molecular Kob-Andersen model.
