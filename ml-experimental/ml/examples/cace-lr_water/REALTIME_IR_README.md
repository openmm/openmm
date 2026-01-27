# Real-time IR Spectrum Computation

This setup allows you to compute IR spectra while the MD simulation is still running.

## How it works

1. **During simulation**: The `run_water_cace_lr.py` script now writes dipole moments to a text file (`*_dipoles_realtime.txt`) in real-time, flushing to disk every 50 steps.

2. **Compute spectrum on-the-fly**: Use `compute_ir_realtime.py` to read the dipole file and compute the IR spectrum, even while the simulation is running.

## Usage

### Start the simulation:
```bash
python run_water_cace_lr.py --molecules 15 --prod 10.0 --dt 0.0005 --output water_test
```

This will create:
- `water_test_dipoles_realtime.txt` - Real-time dipole data (updated every 50 steps)
- `water_test_checkpoint.npz` - Checkpoint files (updated every 500 steps)
- `water_test.npz` - Final output (created at the end)
- `water_test.pdb` - PDB trajectory

### Compute IR spectrum in real-time (one-shot):
```bash
python compute_ir_realtime.py water_test_dipoles_realtime.txt --output water_test_ir
```

### Monitor and auto-update (while simulation runs):
```bash
python compute_ir_realtime.py water_test_dipoles_realtime.txt --monitor --interval 60 --output water_test_ir
```

This will:
- Check for new data every 60 seconds
- Automatically recompute the IR spectrum when new data appears
- Save numbered output files: `water_test_ir_iter001_spectrum.png`, `water_test_ir_iter002_spectrum.png`, etc.

### Options:
- `--window 5.0` - Use only the last 5 ps of data (useful for convergence checking)
- `--max-freq 5000` - Set maximum frequency to plot (default: 4500 cm⁻¹)
- `--interval 30` - Update every 30 seconds in monitor mode

## Files created:

### Real-time dipole file format:
```
# Real-time dipole moment data
# time_ps dipole_x dipole_y dipole_z
0.000000 2.345678 -1.234567 0.987654
0.000500 2.356789 -1.245678 0.998765
...
```

### Output spectrum files:
- `*_spectrum.png` - Plot with ACF and IR spectrum
- `*_spectrum.npz` - Numerical data (frequency, intensity, ACF)

## Tips:

1. **Start monitoring early**: You can start the real-time monitor immediately after starting the simulation, even if only a few data points exist.

2. **Check convergence**: Run with `--window` using different time windows (e.g., last 2 ps, last 5 ps, last 10 ps) to see if the spectrum is converging.

3. **Multiple terminals**: Run the simulation in one terminal and the monitor in another.

## Example workflow:

Terminal 1 (simulation):
```bash
python run_water_cace_lr.py --molecules 15 --prod 20.0 --dt 0.0005 --output long_run
```

Terminal 2 (monitoring):
```bash
python compute_ir_realtime.py long_run_dipoles_realtime.txt --monitor --interval 60
```

The monitor will automatically update the spectrum every minute as new data arrives!
