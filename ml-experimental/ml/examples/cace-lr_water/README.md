# UMA Water System - Machine Learning Potential

Water simulation using FAIRChem UMA (Universal Machine Learning Architecture) potential.

## System Description
- UMA-S-1P1 machine learning potential from FAIRChem
- Flexible water (ML potential captures all interactions including bonds)
- No explicit force field parameters - learned from quantum mechanical data
- GPU-accelerated via PyTorch CUDA
- Optional cavity coupling for polariton chemistry

## Files
- `run_simulation.py` - Main UMA water simulation with cavity coupling support
- `run_baseline.py` - Baseline UMA simulation for IR spectrum validation
- `analyze_spectrum.py` - IR spectrum analysis
- `benchmark_speed.py` - Performance benchmarking

## Quick Start

### Run baseline (no cavity):
```bash
python run_simulation.py --molecules 10 --prod 10 --lambda 0
```

### Run with cavity coupling (OH stretch):
```bash
python run_simulation.py --molecules 10 --lambda 0.01 --cavity-freq 3500 --prod 10
```

### Run with cavity coupling (HOH bend):
```bash
python run_simulation.py --molecules 10 --lambda 0.01 --cavity-freq 1600 --prod 10
```

### Run longer production for IR spectrum:
```bash
python run_baseline.py --molecules 64 --equil 10 --prod 100
```

### Analyze IR spectrum:
```bash
python analyze_spectrum.py uma_water_cavity_lambda0.0100.npz --xlim 0 4500
```

### Benchmark performance:
```bash
python benchmark_speed.py
```

## Parameters
- Temperature: 300 K
- Timestep: 0.5 fs (for flexible bonds)
- Model: uma-s-1p1 (small UMA model)
- Task: omol (organic molecules)
- Cavity frequencies:
  - OH stretch: ~3500 cm⁻¹
  - HOH bending: ~1600 cm⁻¹

## Performance
- ~30-35 steps/s on RTX 4070 (10 water molecules)
- ~1.3-1.5 ns/day with 0.5 fs timestep
- Main bottleneck: UMA neural network inference (~36 ms/step)

## Differences from Classical Water (TIP4P-Ew)
- UMA is an ML potential - no force field parameters needed
- All interactions (bonds, angles, nonbonded) captured by neural network
- More accurate but ~2000x slower than classical force fields
- Uses PythonForce callback mechanism for OpenMM integration

## Expected Results
- OH stretch: ~3200-3700 cm⁻¹ (should be captured by ML potential)
- HOH bending: ~1645 cm⁻¹
- Librational modes: 400-900 cm⁻¹
- With cavity: Rabi splitting around cavity frequency
