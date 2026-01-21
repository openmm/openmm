# Tests Directory Cleanup - Summary

## ✅ Completed Actions

### 1. Created Organized Structure
```
tests/
├── dimer_system/          # Two-component toy model
│   ├── run_simulation.py
│   ├── analyze_spectrum.py
│   ├── benchmark_speed.py
│   └── README.md
│
├── water_system/          # Flexible TIP4P-Ew water
│   ├── run_simulation.py  (from test_cavity_water.py)
│   ├── run_baseline.py    (from test_water_baseline.py)
│   ├── analyze_spectrum.py (from analyze_water_ir.py)
│   ├── benchmark_speed.py (NEW)
│   ├── tip4pew_flexible.xml
│   └── README.md
│
├── _archive/              # Old outputs and legacy scripts
└── README.md              # Main documentation
```

### 2. Archived Legacy Files
Moved to `_archive/`:
- All `.npz`, `.png`, `.log` files (old simulation outputs)
- Old analysis scripts (`calculate_ir_spectrum*.py`, old benchmarks)
- Diagnostic scripts (`diagnose_bonds.py`, `check_bond_vibrations.py`)
- Legacy test scripts (`run_rabi_*.py`, `test_baseline_*.py`)
- Old documentation (`.md` files except main README)
- Shell scripts (`run_ir_simulation.sh`, `run_with_env.sh`)

### 3. Kept in Root Directory
- C++ unit tests (`Test*.h`, `Test*.cpp`) - needed for OpenMM build
- Main simulation scripts referenced by paths
- `CMakeLists.txt` - build configuration
- Data files (`*.dat` - reference data for unit tests)
- Main `README.md` and `TESTS_ORGANIZATION.md`

### 4. Fixed Critical Bug
**Fixed CavityForce initialization** - coupling was never activated!
- OLD: `CavityForce(..., 0.0, ...)` ← Always zero!
- NEW: `CavityForce(..., lambda_coupling, ...)` ← Actual coupling!

### 5. Created New Features
- **Benchmark scripts** for both systems with multi-platform testing
- **Comprehensive READMEs** with quick-start guides
- **Clean separation** between toy model and production simulations

## 🎯 How to Use

### Dimer System (Fast Testing)
```bash
cd tests/dimer_system
python run_simulation.py
python analyze_spectrum.py cavity_diamer_dipole.npz
python benchmark_speed.py
```

### Water System (Production)
```bash
cd tests/water_system

# Step 1: Run baseline
python run_baseline.py --molecules 1000 --prod 100

# Step 2: Run cavity simulations with different λ
python run_simulation.py --lambda 0.001 --molecules 1000 --prod 100
python run_simulation.py --lambda 0.010 --molecules 1000 --prod 100  
python run_simulation.py --lambda 0.050 --molecules 1000 --prod 100

# Step 3: Analyze Rabi splitting
python analyze_spectrum.py water_*.npz --xlim 1400 1800

# Step 4: Benchmark performance
python benchmark_speed.py
```

## 📊 What Changed

### File Consolidation
- **Before**: ~150 files scattered in tests/
- **After**: Clean structure with 2 organized subdirectories

### Script Naming
- **Before**: `test_cavity_water.py`, `test_water_baseline.py`, `analyze_water_ir.py`
- **After**: `run_simulation.py`, `run_baseline.py`, `analyze_spectrum.py`

### Each System Now Has
1. **One run script** - main simulation
2. **One benchmark script** - performance testing  
3. **One analysis script** - IR spectrum with MESA
4. (Water also has baseline script for comparison)

## 🔧 Fixed Issues

1. **Coupling was never activated** - Critical bug fixed
2. **Scattered files** - Now organized by system type
3. **Inconsistent naming** - Standardized to `run_*` and `analyze_*`
4. **Missing benchmarks** - Added for both systems
5. **Poor documentation** - READMEs in each directory

## 📝 Next Steps

The simulations are ready to run with the corrected code. You should:

1. Run baseline water simulation first (for reference)
2. Run 3 cavity simulations at different λ values
3. Analyze to observe Rabi splitting
4. Benchmark to optimize performance

The cavity frequency is currently set to **1600 cm⁻¹** (as you modified).
This targets the HOH bending mode. For OH stretch, change to **3663 cm⁻¹**.
