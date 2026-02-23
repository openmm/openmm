# CACE-LR / UMA Water Systems

Water simulations using ML potentials: CACE-LR (CACE Long Range) and FAIRChem UMA.

## System Description
- UMA-S-1P1 machine learning potential from FAIRChem
- Flexible water (ML potential captures all interactions including bonds)
- No explicit force field parameters - learned from quantum mechanical data
- GPU-accelerated via PyTorch CUDA
- Optional cavity coupling for polariton chemistry

## Files
- `run_water_cace_lr.py` - CACE-LR water with real-time dipole output
- `run_simulation.py` - UMA water with cavity coupling
- `run_baseline.py` - Baseline UMA for IR validation
- `compute_ir_realtime.py` - IR from dipole file while simulation runs
- `analyze_spectrum.py` - IR analysis

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

## Real-time IR (CACE-LR)

`run_water_cace_lr.py` writes dipoles to `*_dipoles_realtime.txt` every 50 steps. Compute spectrum while simulation runs:
```bash
python run_water_cace_lr.py --molecules 15 --prod 10 --output water_test
python compute_ir_realtime.py water_test_dipoles_realtime.txt --output ir
# Or monitor mode: --monitor --interval 60
```

## IR Spectrum Notes

Use unwrapped positions for dipoles (no `enforcePeriodicBox`); molecules across box edges must stay intact. ACF: subtract mean, FFT-based, ω² spectral weight. See `tests/docs/fix-history.md`.

## Expected Results
- OH stretch: ~3200-3700 cm⁻¹
- HOH bending: ~1645 cm⁻¹
- Librational: 400-900 cm⁻¹
- With cavity: Rabi splitting around cavity frequency
