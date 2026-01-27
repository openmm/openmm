# IR Spectrum Calculation - Critical Fixes Applied

## Issues Found and Fixed

### 1. ✅ Wrong Autocorrelation Method (FIXED)
**Problem**: Used direct O(N²) computation instead of FFT-based method
**Impact**: Slower computation, slightly different normalization
**Fix**: Switched to FFT-based autocorrelation matching `tests/water_system/analyze_spectrum.py`

### 2. ✅ Wrong Spectral Weight (FIXED)
**Problem**: Used `ω × |FT[ACF]|²` instead of `ω² × |FT[ACF]|`
**Impact**: Incorrect frequency dependence, wrong peak intensities
**Fix**: Changed to `omega**2 * np.abs(fft_result)` (line ~108)

### 3. ✅ Missing Mean Subtraction (FIXED)
**Problem**: Didn't subtract mean dipole before computing ACF
**Impact**: DC component not removed, affects low-frequency spectrum
**Fix**: Added `subtract_mean=True` to `compute_acf()` function

### 4. ✅ Wrong ACF Normalization (FIXED)
**Problem**: Didn't normalize by decreasing number of samples
**Impact**: ACF doesn't properly account for finite trajectory length
**Fix**: Added `norm = np.arange(N, 0, -1); C_mumu /= norm`

### 5. ✅ **CRITICAL BUG: Wrapped Positions for Dipole** (FIXED)
**Problem**: Used `enforcePeriodicBox=True` when getting positions
**Impact**: **Molecules split across boundaries create artificial huge dipoles!**
- Example: O at x=0.1 nm, H at x=2.9 nm in 3 nm box
- Creates ~2.8 nm artificial separation → huge wrong dipole
- **Completely ruins IR spectrum**

**Fix**: Changed to `state.getState(getPositions=True)` without `enforcePeriodicBox`
This gives **unwrapped positions** where molecules stay intact.

## Verification Against Reference Implementation

All changes now match `tests/water_system/analyze_spectrum.py`:

```python
# Reference implementation (water_system)
def compute_dipole_autocorrelation(dipole_nm, subtract_mean=True):
    dipole = dipole_nm.copy()
    if subtract_mean:
        dipole = dipole - dipole.mean(axis=0)  # ✅ Mean subtraction
    
    N = len(dipole)
    C_mumu = np.zeros(N)
    
    for alpha in range(3):  # Sum over x,y,z
        mu_alpha = dipole[:, alpha]
        nfft = 2 * N
        mu_fft = np.fft.fft(mu_alpha, n=nfft)  # ✅ FFT method
        acf = np.fft.ifft(mu_fft * np.conj(mu_fft))
        C_mumu += acf[:N].real
    
    norm = np.arange(N, 0, -1)  # ✅ Normalization
    C_mumu /= norm
    return C_mumu

def compute_ir_spectrum(...):
    fft_result = rfft(C_windowed)
    omega_rad = 2 * np.pi * freq_hz
    spectrum = np.abs(fft_result) * omega_rad**2  # ✅ ω² weighting
```

## Impact Assessment

### Before Fixes:
- ❌ Wrapped positions: Artificial dipoles from split molecules
- ❌ Wrong spectral weight: Peaks at wrong relative intensities
- ❌ DC component: Low frequency noise
- ❌ Wrong normalization: ACF decay incorrect

### After Fixes:
- ✅ Unwrapped positions: Correct dipole moments
- ✅ Correct spectral weight: Proper ω² dependence
- ✅ Mean subtracted: Clean ACF without DC
- ✅ Proper normalization: Correct ACF decay
- ✅ FFT method: Faster computation (O(N log N) vs O(N²))

## Dipole Moment Calculation

Verified to be correct in both implementations:
```python
# water_system/run_simulation.py
positions = state.getPositions(asNumpy=True)  # ← No enforcePeriodicBox!
dipole += charges[i] * np.array(pos)

# cace-lr_water/cace_dipole_calculator.py
dipole = np.sum(charges_e[:, None] * positions, axis=0)
```

Both compute: **μ = Σ qᵢ · rᵢ** (in e·Å or e·nm)

## Testing Recommendation

Rerun simulations with the fixed code:
```bash
# Fresh simulation with unwrapped positions
python run_water_cace_lr.py --molecules 15 --prod 10.0 --dt 0.0005 --output fixed_run

# Compute spectrum
python compute_ir_realtime.py fixed_run_dipoles_realtime.txt --output fixed_spectrum
```

The spectrum should now be physically meaningful and match expected water vibrational modes.

## Key Takeaway

**Always use unwrapped positions for dipole moment calculations in periodic systems!**

The `enforcePeriodicBox` option is useful for visualization and some analyses, but **not for computing extensive properties like dipole moments** where molecules must remain intact.
