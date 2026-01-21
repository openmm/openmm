# How to Run the Extended Cavity IR Simulation

## Quick Start

```bash
cd /media/extradrive/Trajectories/openmm/tests
./run_ir_simulation.sh
```

This will run the **1 nanosecond simulation at 300 K** with detailed progress reporting.

---

## What the Simulation Does

1. **Creates a diamer system**: 50 diatomic molecules (O-O and N-N dimers) in a periodic box
2. **Adds a cavity photon**: Fictitious particle representing the electromagnetic field mode
3. **Equilibrates for 100 ps**: Coupling is OFF (lambda=0) to let the system reach thermal equilibrium
4. **Production run for 900 ps**: Coupling is ON (lambda=0.001) to observe cavity effects
5. **Records dipole moment**: Every 10 fs for IR spectrum calculation (~100,000 data points)

---

## Progress Output

The simulation reports progress **every 10 ps (10,000 steps)** with the following information:

```
  [ 15.0%] Step  150000/1000000 | Sim time:  150.0 ps | Speed:  245.3 steps/s | ETA:   0.96 hr
         PE: -1245.67 kJ/mol | Cavity: H= 12.34, C=-5.67, D= 2.11
         Cavity: ( 0.0123, -0.0045,  0.0000) nm
```

**Progress indicators:**
- **[ 15.0%]**: Overall completion percentage
- **Step 150000/1000000**: Current step out of total steps
- **Sim time: 150.0 ps**: Simulation time elapsed
- **Speed: 245.3 steps/s**: Performance metric
- **ETA: 0.96 hr**: Estimated time remaining
- **PE**: Potential energy in kJ/mol
- **Cavity: H, C, D**: Harmonic, Coupling, and Dipole self-energy (production phase only)
- **Cavity: (x, y, z)**: Position of cavity particle in nm

---

## Expected Runtime

On the **Reference platform** (CPU):
- **Speed**: ~200-500 steps/second (depends on CPU)
- **Total time**: 0.5-2.5 hours for the full 1 ns simulation

For faster execution, rebuild with CUDA support (requires fixing atomic operations in CUDA kernels).

---

## Output Files

After completion, you'll find:

1. **`cavity_diamer_dipole.npz`**: Main data file containing:
   - `time_ps`: Time points in picoseconds
   - `dipole_nm`: Dipole moment vectors (N × 3 array) in e·nm units
   - `metadata`: Dictionary with simulation parameters

2. **`cavity_sim_output.log`**: Complete log of the simulation

---

## Monitoring a Running Simulation

If you run in the background, monitor progress with:

```bash
tail -f /media/extradrive/Trajectories/openmm/tests/cavity_sim_output.log
```

---

## Next Steps: Calculating IR Spectrum

After the simulation completes, see **`README_IR_SIMULATION.md`** for detailed instructions on:
1. Loading the dipole trajectory
2. Cutting off the first 20 ps
3. Computing the dipole autocorrelation function
4. Fourier transforming to get the IR spectrum

Quick preview:
```python
import numpy as np

# Load data
data = np.load('cavity_diamer_dipole.npz')
time = data['time_ps']
dipole = data['dipole_nm']

# Cut off first 20 ps as requested
idx_cutoff = np.where(time >= 20.0)[0][0]
dipole_cut = dipole[idx_cutoff:]

print(f"Data ready for IR analysis: {len(dipole_cut)} points")
```

---

## Troubleshooting

**If you get import errors:**
The script uses the locally built OpenMM. Make sure:
1. OpenMM was built successfully in `/media/extradrive/Trajectories/openmm/build`
2. The environment variables are set correctly (the run script does this automatically)

**If simulation is too slow:**
The Reference platform is slow but guaranteed to work. For faster execution:
1. Fix the CUDA kernel atomic operations (in progress)
2. Use the CUDA platform instead

**If you need to stop and restart:**
Currently, the simulation doesn't support checkpointing. You would need to start from the beginning. Consider implementing checkpointing for very long runs.
