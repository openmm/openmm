# Water Cavity Molecular Dynamics Simulation

## Overview

This directory contains scripts for simulating **nonthermal aging** in liquid water using strong light-matter coupling in an optical cavity. The cavity mode is tuned to resonate with the OH stretch frequency (~3656 cm⁻¹), creating polaritons that selectively pump vibrational energy and drive structural relaxation.

## Scientific Background

From the paper "Nonthermal Aging of Supercooled Liquids in Optical Cavities":

> When a supercooled liquid is strongly coupled to an optical cavity, the cavity can act as a tunable *non-thermal* energy source, selectively pumping energy into dipole-active vibrations and thereby inducing what we refer to as non-thermal aging, as it drives aging without changing the temperature.

**Key Physics:**
- **OH stretches** = intramolecular vibrations (cavity couples here)
- **Hydrogen bond network** = structural degrees of freedom (compensates via deeper minima)
- **Result**: Structural fictive temperature T_s < T even at constant kinetic temperature

## Files

### Simulation Scripts

1. **`test_cavity_water.py`** - Main cavity MD simulation
   - Creates TIP4P-Ew water box (~1000 molecules)
   - Adds cavity photon mode resonant with OH stretch (3656 cm⁻¹)
   - Runs equilibration (λ=0) then activates coupling
   - Outputs streaming NPZ file with dipole/cavity trajectories

2. **`analyze_water_ir.py`** - IR spectrum analysis
   - Computes dipole autocorrelation function C_μμ(t)
   - Calculates IR absorption spectrum via FFT
   - Identifies polariton peaks (upper/lower branches)
   - Measures Rabi splitting Ω_R vs coupling strength λ

### Data Format

NPZ files contain:
- `time_ps`: Time points in picoseconds
- `dipole_nm`: Dipole moment vectors (N×3 array, units: e·nm)
- `cavity_nm`: Cavity position vectors (N×3 array, units: nm)
- `metadata`: Simulation parameters and status

## Quick Start

### Test Run (Quick Validation)

Run a small test system to verify everything works:

```bash
python test_cavity_water.py --test
```

This runs:
- 216 water molecules
- 1.8 nm box
- 5 ps equilibration + 5 ps production
- Should complete in ~5-10 minutes on GPU

### Production Runs

Run full simulation at different coupling strengths:

```bash
# Weak coupling
python test_cavity_water.py --lambda 0.001 --molecules 1000 --prod 900

# Medium coupling  
python test_cavity_water.py --lambda 0.01 --molecules 1000 --prod 900

# Strong coupling
python test_cavity_water.py --lambda 0.05 --molecules 1000 --prod 900
```

Each run produces: `water_cavity_lambda{value}.npz`

### Analyze IR Spectra

Analyze a single trajectory:

```bash
python analyze_water_ir.py water_cavity_lambda0.0100.npz
```

Compare multiple coupling strengths:

```bash
python analyze_water_ir.py --compare water_cavity_lambda*.npz
```

This produces:
- Individual IR spectrum plots
- Rabi splitting Ω_R vs λ plot (should be linear)

## Expected Results

### Equilibrium Spectrum (λ=0)

Water IR spectrum shows:
- **OH symmetric stretch**: ~3656 cm⁻¹
- **OH antisymmetric stretch**: ~3756 cm⁻¹  
- **HOH bending mode**: ~1645 cm⁻¹

### Cavity-Coupled Spectrum (λ>0)

When cavity couples to OH stretch:
- **Polariton splitting**: OH peak splits into upper/lower branches
- **Rabi splitting**: Ω_R = ω_upper - ω_lower increases linearly with λ
- **Strong coupling regime**: Ω_R >> cavity linewidth

Example for λ = 0.01:
- Lower polariton: ~3500 cm⁻¹
- Upper polariton: ~3800 cm⁻¹
- Rabi splitting: ~300 cm⁻¹

### Nonthermal Aging Dynamics

Time-dependent observables show:
- **Structural slowdown**: Relaxation times increase despite constant T
- **Fictive temperature**: T_s drops below bath temperature T
- **Energy redistribution**: Vibrational → structural (deeper H-bond minima)

## System Parameters

### Default Settings

| Parameter | Value | Description |
|-----------|-------|-------------|
| Water model | TIP4P-Ew | 4-site rigid water |
| N molecules | 1000 | ~3000 atoms total |
| Box size | 3.0 nm | Cubic periodic box |
| Temperature | 300 K | Room temperature |
| ω_c | 3656 cm⁻¹ | Cavity frequency (OH stretch) |
| λ | 0.001-0.05 | Coupling strength |
| Timestep | 1 fs | Integration timestep |
| Equilibration | 100 ps | λ=0 phase |
| Production | 900 ps | λ≠0 phase |

### Cavity Parameters

- **Photon mass**: 1 a.u. electron mass = 5.49×10⁻⁴ amu
- **Cavity frequency**: ω_c = 3656 cm⁻¹ = 0.0169 Hartree
- **Coupling schedule**: 
  - t < 100 ps: λ = 0 (equilibration)
  - t ≥ 100 ps: λ = target value (production)

### Thermostats

- **Molecular thermostat**: BussiThermostat (τ = 1.0 ps)
  - Applied to all water atoms only
  - Preserves microscopic dynamics
- **Cavity thermostat**: Langevin (friction = 0.5 ps⁻¹)
  - Applied to cavity photon mode
  - Models finite polariton lifetime

## Command-Line Options

### test_cavity_water.py

```
usage: test_cavity_water.py [options]

Options:
  --lambda FLOAT      Coupling strength (default: 0.01)
  --molecules INT     Number of water molecules (default: 1000)
  --box FLOAT         Box size in nm (default: 3.0)
  --temp FLOAT        Temperature in K (default: 300)
  --dt FLOAT          Timestep in ps (default: 0.001)
  --equil FLOAT       Equilibration time in ps (default: 100)
  --prod FLOAT        Production time in ps (default: 900)
  --output STR        Output prefix (default: water_cavity)
  --test              Quick test mode (216 molecules, 10 ps)
```

### analyze_water_ir.py

```
usage: analyze_water_ir.py [files...] [options]

Positional:
  files               NPZ file(s) to analyze

Options:
  --compare           Compare multiple files, plot Ω_R vs λ
  --xlim MIN MAX      Frequency range for plot (default: 2500 4500)
  --output FILE       Output filename for plot
```

## Performance

Approximate runtimes on NVIDIA A100 GPU:

| System Size | Timestep | Simulation Time | Wall Time |
|-------------|----------|-----------------|-----------|
| 216 mol (648 atoms) | 1 fs | 10 ps | ~5 min |
| 1000 mol (3000 atoms) | 1 fs | 100 ps | ~30 min |
| 1000 mol (3000 atoms) | 1 fs | 1000 ps | ~5 hr |

*Note: Times include equilibration + production phases*

## Troubleshooting

### CavityForce not available

If you get `AttributeError: CavityForce not available`:
- Make sure you compiled OpenMM with the cavity force implementation
- Check that `libOpenMM.so` includes `CavityForce` and `BussiThermostat`

### BussiThermostat not available

The script will fall back to Langevin-only thermalization if BussiThermostat is not compiled. This is acceptable but may affect long-time dynamics slightly.

### Memory issues with large systems

For very large systems (>2000 molecules):
- Reduce production time or increase output frequency
- Use single precision (`Precision: 'single'`)
- Save trajectory data less frequently

### Slow on CPU

If GPU is not available:
- Simulation will run on CPU (Reference platform)
- Expect ~100× slower performance
- Consider reducing system size for testing

## Scientific Interpretation

### Polariton Formation

Strong coupling creates hybrid light-matter states:
- **Lower polariton**: ω_- = (ω_c - Ω_R/2)
- **Upper polariton**: ω_+ = (ω_c + Ω_R/2)
- **Rabi splitting**: Ω_R ∝ λ√N (collective enhancement)

### Nonthermal Aging Mechanism

1. **Cavity activation** (t = 100 ps): λ switches from 0 → target value
2. **Energy pumping**: Cavity deposits energy into OH stretches
3. **Vibrational heating**: T_v increases to ~500-800 K transiently
4. **Structural compensation**: H-bond network relaxes to deeper minima
5. **Fictive temperature drop**: T_s decreases below T_bath
6. **Aging back**: System slowly returns to equilibrium over ~2-5 ns

### Comparison to Paper Results

The diatomic model in the paper showed:
- Slowdown factor: ~4.75× at λ = 0.141 (T = 100 K)
- Fictive temperature: T_s ≈ 70 K (from T = 100 K)
- Memory timescale: ~2.5 ns (~50× equilibrium relaxation)

Water system should show qualitatively similar behavior:
- Polariton splitting confirms strong coupling regime
- Structural slowdown despite constant T
- Deeper H-bond configurations (lower T_s)

## References

1. Hasyim, M. R., Damiani, A., & Hoffmann, N. M. "Nonthermal Aging of Supercooled Liquids in Optical Cavities" (2026)

2. Li, T. E., Nitzan, A., & Subotnik, J. E. "Cavity molecular dynamics simulations of liquid water under vibrational ultrastrong coupling" *PNAS* **117**, 18324 (2020)

3. Horn, H. W. et al. "Development of an improved four-site water model for biomolecular simulations: TIP4P-Ew" *J. Chem. Phys.* **120**, 9665 (2004)

## Contact

For questions or issues:
- Muhammad R. Hasyim (mh7373@nyu.edu)
- GitHub: [OpenMM Cavity MD Repository]

---

*Last updated: January 2026*
