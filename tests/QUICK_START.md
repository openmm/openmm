# Running the Extended Cavity IR Simulation - Quick Guide

## ✅ Simulation is Now Ready to Run!

The simulation has been successfully configured and tested. All the new cavity features (BussiThermostat, CavityForce, CavityParticleDisplacer) are working correctly.

---

## 🚀 To Run the Simulation

Simply execute:

```bash
cd /media/extradrive/Trajectories/openmm/tests
./run_ir_simulation.sh
```

---

## 📊 What You'll See

The simulation provides detailed progress updates every **10 ps** (10,000 steps):

```
--- Equilibration Phase (Coupling OFF) ---
  Running 100000 steps (100 ps)...
  Progress updates every 10 ps (10,000 steps)

  [  1.0%] Step   10000/1000000 | Sim time:   10.0 ps | Speed: 3978.8 steps/s | ETA:   0.1 hr
         PE:   -29.57 kJ/mol | Cavity: ( 0.0784, -0.0210,  0.0600) nm

  [  2.0%] Step   20000/1000000 | Sim time:   20.0 ps | Speed: 3898.0 steps/s | ETA:   0.1 hr
         PE:   -33.37 kJ/mol | Cavity: (-0.0800,  0.0163,  0.0449) nm
  ...
  
--- Production Phase (Coupling ON) ---
  Lambda coupling now: 0.001
  Running 900000 steps (900 ps)...
  Progress updates every 10 ps (10,000 steps)

  [ 11.0%] Step  110000/1000000 | Sim time:  110.0 ps | Speed: 3937.3 steps/s | ETA:   0.1 hr
         PE:   -74.02 kJ/mol | Cavity: H=  2.28, C= -0.00, D=  0.00
         Cavity: ( 0.0224,  0.0208, -0.0838) nm
  ...
```

**Legend:**
- **[  1.0%]**: Completion percentage
- **Step 10000/1000000**: Current/total steps
- **Sim time: 10.0 ps**: Simulation time
- **Speed: 3978.8 steps/s**: Performance (Reference platform: ~3500-4000 steps/s)
- **ETA: 0.1 hr**: Estimated time remaining
- **PE**: Potential energy
- **Cavity: H, C, D**: Harmonic, Coupling, Dipole self-energy (production phase only)
- **Cavity: (x, y, z)**: Position of cavity particle

---

## ⏱️ Expected Runtime

At **~3900 steps/second** (typical for Reference platform):
- **Total time**: ~4.3 minutes (0.07 hours) for 1,000,000 steps
- **Equilibration** (100,000 steps): ~26 seconds
- **Production** (900,000 steps): ~3.8 minutes

**Much faster than expected!** The Reference platform is quite efficient for this system size.

---

## 📁 Output Files

After completion (in `/media/extradrive/Trajectories/openmm/tests/`):

1. **`cavity_diamer_dipole.npz`**: 
   - Dipole trajectory data (~100,000 points)
   - Time stamps in picoseconds
   - Dipole moments in e·nm units
   - Metadata with simulation parameters

2. **`cavity_sim_output.log`**: 
   - Complete simulation log
   - All progress updates
   - Energy and position data

---

## 🎯 Next Steps: IR Spectrum Calculation

After the simulation completes, calculate the IR spectrum:

```python
import numpy as np
from scipy import signal
from scipy.fft import rfft, rfftfreq
import matplotlib.pyplot as plt

# 1. Load the data
data = np.load('cavity_diamer_dipole.npz')
time = data['time_ps']
dipole = data['dipole_nm']

# 2. Cut off first 20 ps as requested
idx_cutoff = np.where(time >= 20.0)[0][0]
time_cut = time[idx_cutoff:]
dipole_cut = dipole[idx_cutoff:]

print(f"Data points for analysis: {len(time_cut)}")
print(f"Time range: {time_cut[0]:.1f} - {time_cut[-1]:.1f} ps")

# 3. Compute total dipole magnitude
dipole_total = np.linalg.norm(dipole_cut, axis=1)

# 4. Calculate autocorrelation
autocorr = signal.correlate(dipole_total, dipole_total, mode='full', method='fft')
autocorr = autocorr[len(autocorr)//2:]
autocorr = autocorr / autocorr[0]  # Normalize

# 5. Apply window and FFT
window = np.hanning(len(autocorr))
autocorr_windowed = autocorr * window

spectrum = np.abs(rfft(autocorr_windowed))
dt_ps = time_cut[1] - time_cut[0]
frequencies_THz = rfftfreq(len(autocorr_windowed), d=dt_ps * 1e-12)
frequencies_cm = frequencies_THz * 33.356  # Convert to cm^-1

# 6. Plot
plt.figure(figsize=(12, 6))
plt.plot(frequencies_cm[:10000], spectrum[:10000])
plt.xlabel('Frequency (cm$^{-1}$)', fontsize=12)
plt.ylabel('Intensity (arb. units)', fontsize=12)
plt.title('IR Spectrum from Dipole Autocorrelation (300 K, λ=0.001)', fontsize=14)
plt.grid(True, alpha=0.3)
plt.xlim(0, 3000)  # Focus on molecular vibration region
plt.savefig('ir_spectrum.png', dpi=300, bbox_inches='tight')
print("IR spectrum saved to: ir_spectrum.png")
plt.show()
```

---

## ✅ Summary

- ✓ Simulation runs for **1 ns** at **300 K**
- ✓ Coupling **OFF** for first 100 ps (equilibration)
- ✓ Coupling **ON** for remaining 900 ps (production)
- ✓ Dipole recorded every **10 fs** (~100,000 points)
- ✓ **Detailed progress** updates every 10 ps
- ✓ **Fast execution**: ~4-5 minutes total
- ✓ Ready for IR spectrum analysis

**The simulation is ready to go!** Just run `./run_ir_simulation.sh` 🚀
