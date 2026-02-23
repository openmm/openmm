# OpenMM Tests

Test suite for verifying OpenMM functionality and correctness.

## Directory Structure

```
tests/
├── dimer_system/      # Two-component toy model for cavity coupling
├── water_system/      # Flexible TIP4P-Ew water
├── cace-lr_water/    # CACE-LR ML water
├── rpmd/              # UMA/RPMD tests (see rpmd/README.md)
├── protein_system/   # Protein cavity tests
├── uma_ice_rpmd/     # Ice RPMD
├── _archive/         # Legacy outputs
├── docs/             # Fix history, technical notes (docs/fix-history.md)
└── Test*.h, Test*.cpp # C++ unit tests for OpenMM core
```

## Running Tests

### Verification Tests
```bash
cd dimer_system
python test_cavity_coupling.py
```

### C++ Unit Tests
The C++ unit tests (Test*.h, Test*.cpp) are compiled and run as part of the OpenMM test suite.

## Simulation Scripts

**Note:** Simulation scripts (run_simulation.py, analyze_spectrum.py, etc.) have been moved to `examples/cavity/` for clarity. This directory now contains only verification and unit tests.

For actual simulation examples, see:
- `examples/cavity/dimer_system/` - Dimer system simulations
- `examples/cavity/water_system/` - Water system simulations
- `examples/cavity/protein_system/` - Protein system data

## System Requirements

- **OpenMM** with custom BussiThermostat, CavityForce, and CavityParticleDisplacer
- **Python packages**: numpy, scipy, matplotlib, memspectrum
- **GPU recommended**: CUDA or OpenCL for fast water simulations
- **Memory**: ~2 GB for 1000 water molecules

## Key Features

✓ **Dimer System**
- Simple test case for debugging
- Fast execution (~minutes)
- Clear vibrational peaks

✓ **Water System**
- Flexible TIP4P-Ew water model
- OH stretch and HOH bending modes
- Demonstrates Rabi splitting
- Comparison with baseline simulations

✓ **Analysis Tools**
- Maximum Entropy Spectral Analysis (MESA)
- Automatic peak detection
- Rabi splitting quantification
- High-quality publication plots

## Output Files

Each simulation produces:
- `*.npz` - NumPy archive with dipole trajectory, metadata
- `*.png` - IR spectrum plot with DACF and time series

## Documentation

- `dimer_system/README.md`, `water_system/README.md` - Per-system details
- `rpmd/README.md` - UMA/RPMD test suite
- `docs/fix-history.md` - Force constants, cavity coupling, IR fixes

## Notes

- All simulations use **LangevinMiddleIntegrator** with low friction to preserve vibrational dynamics
- **Flexible bonds** are essential for observing molecular vibrations in IR spectra
- Cavity frequency should be tuned to molecular resonances (use MESA on baseline first)
- Production runs should be ≥100 ps for clean spectra (1000 ps for publication quality)

## Troubleshooting

**No vibrational peaks?**
- Check that water is flexible (not rigid)
- Reduce temperature to minimize thermal broadening
- Increase production time for better statistics

**No Rabi splitting?**
- Verify cavity coupling is activated (λ > 0)
- Check cavity frequency matches molecular resonance
- Increase coupling strength or number of molecules

**Slow performance?**
- Use CUDA platform: `--platform CUDA`
- Reduce system size for testing
- Check GPU driver installation
