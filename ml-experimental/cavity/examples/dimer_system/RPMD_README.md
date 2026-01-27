# RPMD Cavity Dimer Simulation

## Overview

The `run_simulation_rpmd.py` script performs Ring Polymer Molecular Dynamics (RPMD) simulations of 250 O-O and N-N dimers coupled to a cavity photon. It uses the PILE_G thermostat (Bussi for centroid + Langevin for internal modes) following i-PI's convention.

## Key Features

- **Quantum Nuclear Effects**: RPMD with configurable number of beads (default: 8)
- **PILE_G Thermostat**: 
  - Bussi stochastic velocity rescaling for centroid mode
  - Langevin thermostat for internal modes
  - GPU-optimized implementation
- **Cavity Coupling**: Light-matter coupling with lambda parameter
- **IR Spectrum**: Computed from centroid dipole moment autocorrelation
- **GPU Acceleration**: CUDA-optimized for high performance

## Usage

### Basic Usage (250 dimers, 8 beads, lambda=0.0700)
```bash
cd /media/extradrive/Trajectories/openmm/tests/dimer_system
python3 run_simulation_rpmd.py
```

### Custom Parameters
```bash
python3 run_simulation_rpmd.py \
    --dimers 250 \
    --beads 8 \
    --lambda 0.0700 \
    --temp 100 \
    --dt 0.001 \
    --equil 100 \
    --prod 900 \
    --cavity-freq 1560
```

### Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `--dimers` | 250 | Number of dimer molecules |
| `--beads` | 8 | Number of RPMD beads (quantum nuclei) |
| `--lambda` | 0.0700 | Cavity coupling strength |
| `--temp` | 100 | Temperature in Kelvin |
| `--dt` | 0.001 | Timestep in ps (1 fs) |
| `--equil` | 100 | Equilibration time in ps |
| `--prod` | 900 | Production time in ps |
| `--cavity-freq` | 1560 | Cavity frequency in cm⁻¹ |
| `--no-dipole` | False | Disable dipole output (benchmark mode) |

## Output

The simulation generates an NPZ file: `cavity_diamer_lambda{lambda:.4f}.npz`

### Data Structure
```python
import numpy as np

data = np.load('cavity_diamer_lambda0.0700.npz', allow_pickle=True)
time_ps = data['time_ps']           # Time in picoseconds
dipole_nm = data['dipole_nm']       # Centroid dipole moment (e*nm)
metadata = data['metadata'].item()  # Simulation parameters
```

### Metadata Fields
- `temperature_K`: Simulation temperature
- `num_molecules`: Number of dimers
- `num_beads`: Number of RPMD beads
- `cavity_freq_cm`: Cavity frequency
- `lambda_coupling`: Coupling strength
- `thermostat`: 'PILE_G'
- `centroid_friction`: Bussi thermostat friction
- `dt_ps`: Timestep
- `total_steps`: Total simulation steps
- `elapsed_time_s`: Wall-clock time

## Computing the IR Spectrum

### Step 1: Load Data
```python
import numpy as np
from scipy import signal

data = np.load('cavity_diamer_lambda0.0700.npz', allow_pickle=True)
time = data['time_ps']
dipole = data['dipole_nm']  # Shape: (n_timesteps, 3)
metadata = data['metadata'].item()
dt = metadata['dt_ps']
```

### Step 2: Calculate Dipole Autocorrelation
```python
# Use only production phase (skip equilibration)
equilibration_ps = 100.0  # ps
idx_start = np.where(time >= equilibration_ps)[0][0]
time_prod = time[idx_start:]
dipole_prod = dipole[idx_start:]

# Compute autocorrelation for each component
n_points = len(dipole_prod)
autocorr = np.zeros(n_points)

for i in range(3):  # x, y, z components
    d = dipole_prod[:, i]
    d = d - d.mean()  # Remove DC offset
    
    # FFT-based autocorrelation
    fft_d = np.fft.fft(d, n=2*n_points)
    power = np.abs(fft_d)**2
    autocorr_component = np.fft.ifft(power)[:n_points].real
    autocorr += autocorr_component

autocorr /= autocorr[0]  # Normalize
```

### Step 3: Fourier Transform to IR Spectrum
```python
# Apply window function
window = signal.windows.blackman(len(autocorr))
autocorr_windowed = autocorr * window

# FFT to frequency domain
spectrum = np.fft.rfft(autocorr_windowed)
freqs_ps = np.fft.rfftfreq(len(autocorr), d=dt)

# Convert to cm⁻¹
c_cm_per_ps = 33.356  # Speed of light in cm/ps
freqs_cm = freqs_ps * c_cm_per_ps

# IR intensity
intensity = np.abs(spectrum)**2
```

### Step 4: Plot
```python
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(freqs_cm, intensity)
plt.xlim(0, 3000)
plt.xlabel('Frequency (cm$^{-1}$)')
plt.ylabel('IR Intensity (arb. units)')
plt.title(f'IR Spectrum (RPMD, {metadata["num_beads"]} beads, λ={metadata["lambda_coupling"]:.4f})')
plt.grid(True, alpha=0.3)
plt.savefig('ir_spectrum_rpmd.png', dpi=300, bbox_inches='tight')
plt.show()
```

## RPMD Theory

### Ring Polymer Path Integral
RPMD represents quantum nuclei as classical "ring polymers" with P beads connected by harmonic springs. The centroid represents the average position:

$$\mathbf{R}_{\text{centroid}} = \frac{1}{P} \sum_{k=1}^P \mathbf{R}_k$$

### PILE_G Thermostat
Following i-PI convention:
- **Centroid mode (k=0)**: Bussi stochastic velocity rescaling
- **Internal modes (k≠0)**: PILE Langevin thermostat

This ensures proper quantum statistics while maintaining computational efficiency.

### Quantum Effects
At 100 K with O-O stretching at 1560 cm⁻¹:
- Classical limit: ℏω/kT → 0
- Quantum regime: ℏω/kT ≈ 2.2
- Zero-point energy: ½ℏω = 780 cm⁻¹

RPMD captures:
- Zero-point motion
- Quantum tunneling
- Nuclear delocalization
- Proper quantum statistics

## Performance

### Expected Runtime (CUDA GPU)
- **250 dimers (500 atoms), 8 beads**: ~2-3 hours for 1 ns
- **GPU acceleration**: ~100-150 ns/day
- **Memory**: ~2-3 GB GPU memory

### Scaling with Beads
| Beads | Computational Cost | Quantum Accuracy |
|-------|-------------------|------------------|
| 1 | 1× (classical) | No quantum effects |
| 4 | 4× | Basic quantum |
| 8 | 8× | Good quantum |
| 16 | 16× | Converged quantum |
| 32 | 32× | Overkill for 100 K |

## Comparison with Classical MD

### Classical MD (`run_simulation.py`)
- Single trajectory
- No quantum effects
- Bussi thermostat for all particles
- Faster (~300 ns/day)

### RPMD (`run_simulation_rpmd.py`)
- Multiple bead trajectories
- Quantum nuclear effects
- PILE_G thermostat (centroid + internal)
- Slower (~100-150 ns/day with 8 beads)
- More accurate for light atoms at low T

## Troubleshooting

### Out of Memory
Reduce number of beads: `--beads 4`

### Slow Performance
1. Check GPU usage: `nvidia-smi`
2. Reduce number of dimers: `--dimers 100`
3. Disable dipole output: `--no-dipole` (benchmark mode)

### Numerical Instability
Reduce timestep: `--dt 0.0005` (0.5 fs)

## References

1. **RPMD Theory**: Habershon et al., Annu. Rev. Phys. Chem. 64, 387 (2013)
2. **PILE Thermostat**: Ceriotti et al., J. Chem. Phys. 133, 124104 (2010)
3. **Bussi Thermostat**: Bussi et al., J. Chem. Phys. 126, 014101 (2007)
4. **i-PI Software**: Kapil et al., Comput. Phys. Commun. 236, 214 (2019)
