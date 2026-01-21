# Cavity OpenMM IR Spectrum Simulation

## Overview

The extended test script `test_cavity_diamer.py` now runs a **1 nanosecond** simulation at **300 K (room temperature)** with dipole moment output for IR spectrum calculation.

## Simulation Parameters

- **Temperature**: 300 K (room temperature)
- **Total Time**: 1 ns (1000 ps)
  - **Equilibration**: 100 ps (coupling OFF, lambda=0)
  - **Production**: 900 ps (coupling ON, lambda=0.001)
- **Timestep**: 1 fs (0.001 ps)
- **Dipole Output Interval**: 10 fs (0.01 ps)
  - This gives ~100,000 data points for good frequency resolution

## System Details

- 50 diatomic molecules (O-O and N-N dimers)
- 80% O-O, 20% N-N
- Harmonic bonds, LJ interactions, Coulomb interactions
- 1 cavity photon particle
- Periodic boundary conditions (2.5 nm box)

## Output File

The simulation saves dipole trajectory to: **`cavity_diamer_dipole.npz`**

### Contents:
```python
data = np.load('cavity_diamer_dipole.npz')
time_ps = data['time_ps']        # (N,) array of times in picoseconds
dipole_nm = data['dipole_nm']    # (N, 3) array of dipole in e*nm units
metadata = data['metadata']       # dict with simulation parameters
```

## IR Spectrum Calculation

### Step 1: Load and Cut Data

```python
import numpy as np

# Load data
data = np.load('cavity_diamer_dipole.npz')
time = data['time_ps']
dipole = data['dipole_nm']

# Cut off first 20 ps as requested
idx_cutoff = np.where(time >= 20.0)[0][0]
time_cut = time[idx_cutoff:]
dipole_cut = dipole[idx_cutoff:]

print(f"Data points after cutoff: {len(time_cut)}")
print(f"Time range: {time_cut[0]:.2f} - {time_cut[-1]:.2f} ps")
```

### Step 2: Compute Dipole Autocorrelation

```python
from scipy import signal

# Total dipole magnitude at each timestep
dipole_total = np.linalg.norm(dipole_cut, axis=1)

# Autocorrelation using FFT (efficient for large datasets)
autocorr = signal.correlate(dipole_total, dipole_total, mode='full', method='fft')
autocorr = autocorr[len(autocorr)//2:]  # Keep only positive lags
autocorr = autocorr / autocorr[0]  # Normalize

# Time axis for autocorrelation
dt_ps = time_cut[1] - time_cut[0]
time_autocorr = np.arange(len(autocorr)) * dt_ps
```

### Step 3: Fourier Transform to Get IR Spectrum

```python
from scipy.fft import rfft, rfftfreq

# Apply window function to reduce spectral leakage
window = np.hanning(len(autocorr))
autocorr_windowed = autocorr * window

# FFT to get spectrum
spectrum = np.abs(rfft(autocorr_windowed))
frequencies_THz = rfftfreq(len(autocorr_windowed), d=dt_ps * 1e-12)  # Convert ps to s

# Convert THz to cm^-1
frequencies_cm = frequencies_THz * 33.356  # 1 THz = 33.356 cm^-1

# Plot
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.plot(frequencies_cm[:10000], spectrum[:10000])  # Plot up to ~3000 cm^-1
plt.xlabel('Frequency (cm$^{-1}$)')
plt.ylabel('Intensity (arb. units)')
plt.title('IR Spectrum from Dipole Autocorrelation')
plt.grid(True)
plt.savefig('ir_spectrum.png', dpi=300)
plt.show()
```

## Running the Simulation

### Using the Run Script (Recommended)

```bash
cd /media/extradrive/Trajectories/openmm/tests
./run_ir_simulation.sh
```

### Manual Execution

```bash
cd /media/extradrive/Trajectories/openmm/tests
export OPENMM_BUILD_DIR="/media/extradrive/Trajectories/openmm/build"
export PYTHONPATH="${OPENMM_BUILD_DIR}/python:${PYTHONPATH}"
export LD_LIBRARY_PATH="${OPENMM_BUILD_DIR}:${LD_LIBRARY_PATH}"
python test_cavity_diamer.py
```

**Note**: The simulation will take considerable time (~hours) as it runs 1 million steps on the Reference platform. For production runs, consider using CUDA platform (requires fixing the CUDA kernel atomic operations first).

### Progress Monitoring

The simulation provides detailed progress updates every **10 ps (10,000 steps)** including:

- **Progress percentage**: Overall completion (0-100%)
- **Step counter**: Current step / total steps
- **Simulation time**: Current simulation time in picoseconds
- **Performance**: Steps per second
- **ETA**: Estimated time remaining in hours
- **Energy**: Potential energy and cavity energy components (during production)
- **Cavity position**: x, y, z coordinates of the cavity particle

**Example output:**
```
  [ 15.0%] Step  150000/1000000 | Sim time:  150.0 ps | Speed:  245.3 steps/s | ETA:   0.96 hr
         PE: -1245.67 kJ/mol | Cavity: H= 12.34, C=-5.67, D= 2.11
         Cavity: ( 0.0123, -0.0045,  0.0000) nm
```

### Output Files

- **`cavity_diamer_dipole.npz`**: Main output with dipole trajectory
- **`cavity_sim_output.log`**: Complete simulation log (when using run script)

## Key Features

1. **Time-varying coupling**: Coupling is OFF during equilibration (0-100 ps) and ON during production (100-1000 ps)
2. **Room temperature**: 300 K for realistic conditions
3. **High-frequency sampling**: 10 fs intervals provide frequency resolution up to ~5000 cm^-1
4. **Complete trajectory**: All molecular dipole moments saved for post-processing

## Expected Results

The IR spectrum should show:
- Molecular vibrations (bond stretches, typically 1000-3000 cm^-1)
- Cavity mode influence on molecular frequencies
- Possible collective effects from cavity-molecule coupling

The 20 ps cutoff ensures the system is equilibrated before data collection for spectrum analysis.
