# Tests Directory - Final Organization Summary

## ✅ Completed Tasks

### 1. Removed Duplicate Test Scripts
All test scripts moved to `_archive/` to avoid confusion:
- `test_cavity_diamer.py` → in `dimer_system/run_simulation.py`
- `test_cavity_water.py` → in `water_system/run_simulation.py`
- `test_water_baseline.py` → in `water_system/run_baseline.py`
- `analyze_water_ir.py` → in `water_system/analyze_spectrum.py`
- `tip4pew_flexible.xml` → in `water_system/tip4pew_flexible.xml`

### 2. Added Cavity Frequency Flag
**Water System now supports `--cavity-freq` flag:**

```bash
# HOH bending mode (1600 cm⁻¹, default)
python run_simulation.py --lambda 0.010 --cavity-freq 1600 --prod 100

# OH stretch mode (3663 cm⁻¹)
python run_simulation.py --lambda 0.050 --cavity-freq 3663 --prod 100
```

**Dimer System:**
- Cavity frequency is hardcoded to 0.00913 Hartree (~2000 cm⁻¹)
- Can be changed by editing line 286 in `dimer_system/run_simulation.py`

### 3. Critical Bug Fixed
**Cavity coupling was never activated in previous versions!**
- OLD: `CavityForce(..., 0.0, ...)` ← Always zero coupling
- NEW: `CavityForce(..., lambda_coupling, ...)` ← Uses actual λ value

This explains why previous λ variations showed no effect on IR spectra.

## 📁 Final Directory Structure

```
tests/
├── dimer_system/              [4 files - Toy model]
│   ├── run_simulation.py      
│   ├── analyze_spectrum.py    
│   ├── benchmark_speed.py     
│   └── README.md
│
├── water_system/              [6 files - Production]
│   ├── run_simulation.py      ← --cavity-freq flag added
│   ├── run_baseline.py        
│   ├── analyze_spectrum.py    
│   ├── benchmark_speed.py     
│   ├── tip4pew_flexible.xml   
│   └── README.md
│
├── _archive/                  [100+ old files]
├── Test*.h, Test*.cpp         [C++ unit tests for build]
└── README.md                  [Main documentation]
```

## 🚀 Quick Start Examples

### Water System - HOH Bending Mode
```bash
cd tests/water_system

# Run baseline (no cavity)
python run_baseline.py --molecules 1000 --prod 100

# Run with different coupling strengths
python run_simulation.py --cavity-freq 1600 --lambda 0.001 --prod 100
python run_simulation.py --cavity-freq 1600 --lambda 0.010 --prod 100
python run_simulation.py --cavity-freq 1600 --lambda 0.050 --prod 100

# Analyze Rabi splitting
python analyze_spectrum.py water_*.npz --xlim 1400 1800
```

### Water System - OH Stretch Mode
```bash
cd tests/water_system

# Need new baseline at 3663 cm⁻¹ or use existing baseline
python run_simulation.py --cavity-freq 3663 --lambda 0.001 --prod 100
python run_simulation.py --cavity-freq 3663 --lambda 0.010 --prod 100
python run_simulation.py --cavity-freq 3663 --lambda 0.050 --prod 100

# Analyze Rabi splitting
python analyze_spectrum.py water_*.npz --xlim 3400 3900
```

### Dimer System
```bash
cd tests/dimer_system

# Run simulation (cavity freq hardcoded to ~2000 cm⁻¹)
python run_simulation.py

# Analyze spectrum
python analyze_spectrum.py cavity_diamer_dipole.npz

# Benchmark performance
python benchmark_speed.py
```

## 🔍 Cavity Frequency Reference

### Water Modes (from MESA analysis)
- **HOH bending**: 1654 cm⁻¹ (default: 1600 cm⁻¹)
  - Collective mode, easier to couple
  - Better for observing Rabi splitting
  - Recommended for initial studies

- **OH stretch**: 3663 cm⁻¹
  - Higher frequency, stronger coupling needed
  - More challenging but demonstrates versatility
  - Requires larger λ for observable splitting

### Dimer Modes
- **O-O stretch**: ~1740 cm⁻¹
- **N-N stretch**: ~2280 cm⁻¹
- Current cavity: ~2000 cm⁻¹ (0.00913 Hartree)

## 📊 Command-Line Flags

### run_simulation.py (Water System)
```
--lambda         Coupling strength (default: 0.01)
--cavity-freq    Cavity frequency in cm⁻¹ (default: 1600)
--molecules      Number of water molecules (default: 1000)
--box            Box size in nm (default: 3.0)
--temp           Temperature in K (default: 300)
--dt             Timestep in ps (default: 0.0005 = 0.5 fs)
--equil          Equilibration time in ps (default: 100)
--prod           Production time in ps (default: 900)
--output         Output prefix (default: water_cavity)
--test           Quick test run (216 molecules, 10 ps)
```

### run_baseline.py (Water System)
```
--molecules      Number of water molecules (default: 1000)
--box            Box size in nm (default: 3.0)
--temp           Temperature in K (default: 300)
--dt             Timestep in ps (default: 0.0005)
--equil          Equilibration time in ps (default: 100)
--prod           Production time in ps (default: 100)
--output         Output prefix (default: water_baseline)
```

## 💡 Recommendations

1. **Start with HOH bending** (1600 cm⁻¹) - easier to observe Rabi splitting
2. **Use default parameters** for initial tests (1000 molecules, 100 ps)
3. **Run baseline first** to identify actual peak frequencies
4. **Test different λ values** (0.001, 0.010, 0.050) to see scaling
5. **Increase production time** to 1000 ps for publication-quality spectra
6. **Use CUDA platform** for faster simulations on GPU

## ⚠️ Important Notes

- All old output files are in `_archive/` and can be deleted
- The root `tests/` directory should stay clean (only organized subdirs + C++ tests)
- Flexible water is essential - rigid water won't show vibrational modes
- Timestep must be ≤0.5 fs for stable flexible bond integration
- Cavity coupling is now ON from t=0 (not gradually activated)
