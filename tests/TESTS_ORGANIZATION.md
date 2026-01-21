# Tests Directory Organization

This directory contains two main test systems for cavity-coupled molecular dynamics:

## 1. Dimer System (`dimer_system/`)
Simple two-component toy model for testing cavity coupling implementation.
- Fast to run (~minutes)
- Good for debugging and development
- See `dimer_system/README.md` for details

## 2. Water System (`water_system/`)
Realistic flexible TIP4P-Ew water with cavity coupling.
- Production-quality simulations
- Demonstrates Rabi splitting
- See `water_system/README.md` for details

## C++ Unit Tests
The `.h` and `.cpp` files (TestBussiThermostat.h, TestCavityForce.h, etc.) are C++ unit tests for the OpenMM core library. These are compiled and run separately from the Python integration tests.

## Legacy Files
The root `tests/` directory contains legacy test scripts and output files. These are kept for reference but new work should use the organized subdirectories.

## Quick Start

### Dimer System:
```bash
cd dimer_system
python run_simulation.py         # Run simulation
python analyze_spectrum.py *.npz  # Analyze results
python benchmark_speed.py         # Test performance
```

### Water System:
```bash
cd water_system
python run_baseline.py --molecules 1000 --prod 100     # Baseline
python run_simulation.py --lambda 0.010 --prod 100     # Cavity-coupled
python analyze_spectrum.py *.npz                       # Analyze & compare
python benchmark_speed.py                              # Test performance
```
