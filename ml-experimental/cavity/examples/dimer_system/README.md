# Dimer System Simulations

This directory contains simulations of a two-component diatomic dimer system (O-O and N-N dimers) coupled to an optical cavity.

## System Parameters

**Molecular frequencies:**
- O-O stretch: **1560 cm⁻¹** (most abundant, 80% of dimers)
- N-N stretch: **2325 cm⁻¹** (less abundant, 20% of dimers)

**Default simulation settings:**
- Number of dimers: 250 (200 O-O + 50 N-N)
- Temperature: 100 K
- Timestep: 1.0 fs
- Equilibration: 20 ps
- Production: 100 ps
- Cavity frequency: 1560 cm⁻¹ (resonant with O-O stretch)
- Coupling strength (λ): 0.001 a.u.

## Scripts

### `run_simulation.py`
Run the cavity-coupled dimer simulation with custom parameters.

```bash
# Basic run (default parameters)
python run_simulation.py

# Custom parameters
python run_simulation.py --dimers 500 --lambda 0.005 --cavity-freq 2325 --prod 200

# Speed benchmark (no dipole output)
python run_simulation.py --no-dipole --prod 10
```

**Flags:**
- `--dimers N`: Number of dimers (default: 250)
- `--lambda λ`: Coupling strength in a.u. (default: 0.001)
- `--temp T`: Temperature in K (default: 100)
- `--dt Δt`: Timestep in ps (default: 0.001 = 1 fs)
- `--equil t`: Equilibration time in ps (default: 20)
- `--prod t`: Production time in ps (default: 100)
- `--cavity-freq ω`: Cavity frequency in cm⁻¹ (default: 1560)
- `--no-dipole`: Disable dipole output for speed benchmarking

### `analyze_spectrum.py`
Calculate and plot the IR spectrum using Maximum Entropy Spectral Analysis (MESA).

```bash
# Analyze most recent simulation
python analyze_spectrum.py

# Analyze specific file
python analyze_spectrum.py cavity_diamer_lambda0.0010.npz
```

**Outputs:**
- `cavity_ir_spectrum_lambda{λ}_full.png`: 3-panel plot (autocorrelation, spectrum, dipole trajectory)
- `cavity_ir_spectrum_lambda{λ}_hires.png`: High-resolution spectrum (0-4000 cm⁻¹)
- `cavity_ir_spectrum_lambda{λ}_data.npz`: Processed data

### `benchmark_speed.py`
Benchmark simulation speed without dipole output overhead.

```bash
python benchmark_speed.py
```

## Expected Results

**Without cavity coupling (λ=0):**
- Two peaks at 1560 cm⁻¹ (O-O) and 2325 cm⁻¹ (N-N)
- Peak heights reflect 80/20 mixture ratio

**With cavity coupling (λ > 0, resonant with O-O):**
- O-O peak splits into upper and lower polaritons (Rabi splitting)
- N-N peak remains largely unchanged (off-resonant)
- Splitting magnitude: ΔE ≈ 2λ√(N_OO) ωc

## Force Field Parameters

**Bond parameters** (from cav-hoomd):
```python
# HOOMD and OpenMM both use: E = 0.5 * k * (r - r0)²
k_OO = 0.73204 Hartree/Bohr²  →  1560 cm⁻¹
k_NN = 1.4325 Hartree/Bohr²   →  2325 cm⁻¹

r0_OO = 2.281655158 Bohr
r0_NN = 2.0743522177 Bohr
```

**LJ parameters:**
```python
σ_O = 0.3 nm, ε_O = 0.5 kJ/mol
σ_N = 0.25 nm, ε_N = 0.3 kJ/mol
```

**Charges:**
```python
q = ±0.3 e
```

## References

- cav-hoomd repository: Force field and system definition
- Maximum Entropy Spectral Analysis: https://maximum-entropy-spectrum.readthedocs.io/
