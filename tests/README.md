# OpenMM Cavity-Coupled MD Tests

Clean, organized test suite for cavity-coupled molecular dynamics simulations.

## Directory Structure

```
tests/
├── dimer_system/          # Toy two-component dimer system
│   ├── run_simulation.py  # Main simulation
│   ├── analyze_spectrum.py # IR analysis with MESA
│   ├── benchmark_speed.py  # Performance testing
│   └── README.md
│
├── water_system/          # Flexible TIP4P-Ew water with cavity
│   ├── run_simulation.py  # Cavity-coupled simulation
│   ├── run_baseline.py    # Baseline (no cavity) for comparison
│   ├── analyze_spectrum.py # Rabi splitting analysis
│   ├── benchmark_speed.py  # Performance testing
│   ├── tip4pew_flexible.xml # Force field parameters
│   └── README.md
│
├── _archive/              # Legacy files and old outputs
└── Test*.h, Test*.cpp     # C++ unit tests for OpenMM core
```

## Quick Start

### Dimer System (Fast, ~5 minutes)
```bash
cd dimer_system
python run_simulation.py
python analyze_spectrum.py cavity_diamer_dipole.npz
```

### Water System (Slower, ~30 minutes for 100 ps)
```bash
cd water_system

# Run baseline (no cavity)
python run_baseline.py --molecules 1000 --prod 100

# Run with cavity coupling
python run_simulation.py --lambda 0.010 --molecules 1000 --prod 100

# Analyze and compare
python analyze_spectrum.py *.npz --xlim 1400 1800
```

### Performance Benchmarking
```bash
cd water_system  # or dimer_system
python benchmark_speed.py
```

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

- `dimer_system/README.md` - Dimer system details
- `water_system/README.md` - Water system details
- `TESTS_ORGANIZATION.md` - Organization notes

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
