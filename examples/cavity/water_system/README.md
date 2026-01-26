# Water System - Cavity-Coupled Flexible TIP4P-Ew

Realistic water simulation with cavity coupling to study Rabi splitting and vibrational strong coupling.

## System Description
- Flexible TIP4P-Ew water model (allows OH stretch vibrations)
- 1000 water molecules (3000 atoms)
- Cavity resonance tuned to HOH bending (~1600 cm⁻¹) or OH stretch (~3663 cm⁻¹)
- Variable coupling strength λ to observe Rabi splitting

## Files
- `run_simulation.py` - Main cavity-coupled simulation
- `run_baseline.py` - Baseline (uncoupled) simulation for comparison
- `analyze_spectrum.py` - IR spectrum analysis and Rabi splitting detection
- `benchmark_speed.py` - Performance benchmarking
- `tip4pew_flexible.xml` - Flexible TIP4P-Ew force field parameters

## Quick Start

### Run baseline (no cavity):
```bash
python run_baseline.py --molecules 1000 --prod 100
```

### Run cavity-coupled simulation:
```bash
# Weak coupling (HOH bending mode)
python run_simulation.py --lambda 0.001 --cavity-freq 1600 --molecules 1000 --prod 100

# Medium coupling (HOH bending mode)
python run_simulation.py --lambda 0.010 --cavity-freq 1600 --molecules 1000 --prod 100

# Strong coupling (HOH bending mode)
python run_simulation.py --lambda 0.050 --cavity-freq 1600 --molecules 1000 --prod 100

# Target OH stretch instead (higher frequency)
python run_simulation.py --lambda 0.050 --cavity-freq 3663 --molecules 1000 --prod 100
```

### Analyze Rabi splitting:
```bash
python analyze_spectrum.py water_flexible_lambda*.npz --xlim 1400 1800
```

### Benchmark performance:
```bash
python benchmark_speed.py
```

## Parameters
- Temperature: 300 K
- Timestep: 0.5 fs (for flexible bonds)
- Equilibration: 100 ps
- Production: 100 ps (for quick tests) or 1000 ps (for publication-quality)
- Cavity frequencies:
  - HOH bending: 1600-1660 cm⁻¹
  - OH stretch: 3660-3670 cm⁻¹

## Expected Results
- Baseline: Single peak at vibrational frequency
- With cavity: Rabi splitting into upper and lower polariton branches
- Splitting increases with √(N×λ)
