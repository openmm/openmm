# Dimer System Simulations

This directory contains simulations of a two-component diatomic dimer system (O-O and N-N dimers) coupled to an optical cavity.

## Install (OpenMM from this repo)

To run `run_simulation.py` you need OpenMM built and installed from this repository (with CavityForce and Bussi thermostat). From the repo root:

```bash
# Optional: use a conda env
conda activate base
# Full build + install (OpenMM, optional Fairchem/openmm-ml); see scripts/README_INSTALL_BASE.md
bash scripts/install_openmm_fairchem_base.sh
```

If you only need OpenMM (no Fairchem): configure and build OpenMM, then `make install` and install the Python package from `build/python` (e.g. `pip install build/python`). The script uses CUDA when available and falls back to Reference (CPU) otherwise.

## System Parameters

**Molecular frequencies:**
- O-O stretch: **1560 cm⁻¹** (most abundant, 80% of dimers)
- N-N stretch: **2325 cm⁻¹** (less abundant, 20% of dimers)

**Default simulation settings (cav-hoomd parity: 100 K, 40 Bohr box, 250 dimers):**
- Number of dimers: 250 (200 O-O + 50 N-N)
- Box: 2.117 nm (40 Bohr)
- Temperature: 100 K
- Timestep: 1.0 fs
- Equilibration: 100 ps
- Production: 900 ps
- Cavity frequency: 1560 cm⁻¹ (resonant with O-O stretch)
- Coupling strength (λ): 0.001 a.u. (default g = λ√N ≈ 0.0158 for N=250)

## Scripts

### `run_simulation.py`
Run the cavity-coupled dimer simulation with custom parameters. Default box (40 Bohr), temperature (100 K), and coupling (λ=0.001) match the cav-hoomd "100 K, 40 Bohr, 250 dimers" regime; use `--lambda 0.001` or the default `--g` to compare with cav-hoomd `--coupling 1e-3`.

```bash
# Basic run (default parameters)
python run_simulation.py

# Custom parameters
python run_simulation.py --dimers 500 --lambda 0.005 --cavity-freq 2325 --prod 200

# Speed benchmark (no dipole output)
python run_simulation.py --no-dipole --prod 10

# F(k,t) + PDB, constant density, no cavity coupling (g=0), no dipole output
python run_simulation.py --dimers 250 --g 0.0 --temp 100.0 --dt 0.001 --equil 200 --prod 3000 \
  --cavity-freq 1555 --report-interval 10000 --pdb-interval 10000 --enable-fkt \
  --fkt-output-period-ps 1.0 --fkt-ref-interval-ps 100.0 --constant-density --no-dipole --pdb ./test.pdb
```

**Flags:**
- `--dimers N`: Number of dimers (default: 250)
- `--lambda λ`: Coupling λ in a.u.; overrides `--g` when set (e.g. 0.001 for cav-hoomd parity)
- `--g g`: Coupling g = λ√N when `--lambda` not set (default: 0.0158 ⇒ λ≈0.001 for N=250)
- `--box-size L`: Box edge in nm (default: 2.117 = 40 Bohr). Use `--constant-density` to scale with N.
- `--temp T`: Temperature in K (default: 100)
- `--dt Δt`: Timestep in ps (default: 0.001 = 1 fs)
- `--equil t`: Equilibration time in ps (default: 100)
- `--prod t`: Production time in ps (default: 900)
- `--cavity-freq ω`: Cavity frequency in cm⁻¹ (default: 1560)
- `--no-dipole`: Disable dipole output for speed benchmarking

### Force comparison: OpenMM vs cav-hoomd (reference)

To run OpenMM cavity and cav-hoomd **side by side** and compare per-particle forces (cav-hoomd is the reference; OpenMM must match within a chosen tolerance):

1. **Shared configuration:** Use a GSD snapshot (e.g. from cav-hoomd `initlattice_equilibrium.py`) with the same N, box, and particle order.

2. **OpenMM forces:** From this directory:
   ```bash
   python compare_openmm_cav_hoomd_forces.py init-0.gsd --lambda 0.001 --output-forces openmm_f.npz
   ```
   This builds the OpenMM system from the GSD (positions and box in Bohr → nm), adds the cavity force, and writes forces (N×3) in kJ/(mol·nm) to `openmm_f.npz`.

3. **cav-hoomd forces:** In an environment where cav-hoomd is installed, run a force-dump that loads the same GSD, runs one step, and writes forces (N×3) in **Hartree/Bohr** to an NPZ (e.g. `hoomd_forces.npz`). The reference forces **must be in GSD particle (tag) order** so the comparison script can align them to OpenMM bond order. Use `dump_hoomd_forces_from_gsd.py` in this directory; it uses HOOMD's `cpu_local_snapshot` + `rtag` to write forces in GSD (tag) order. If you use another script, ensure it outputs tag order.

4. **Compare:** Pass the cav-hoomd NPZ to the comparison script. All forces are converted to Hartree/Bohr for comparison; default tolerance is max |F_openmm − F_hoomd| ≤ 1e−5 Ha/Bohr.
   ```bash
   python compare_openmm_cav_hoomd_forces.py init-0.gsd --lambda 0.001 --hoomd-forces-npz hoomd_forces.npz [--tol 1e-5]
   ```
   Exit code 0 means OpenMM matches the reference within tolerance; exit code 1 means the test fails.

**Scripts:**
- `compare_openmm_cav_hoomd_forces.py`: Loads GSD, builds OpenMM system from it, computes forces; optionally compares to `--hoomd-forces-npz` and asserts tolerance. Use `--diagnose-order` to try alignment variants and see which minimizes error.
- `dump_hoomd_forces_from_gsd.py`: Dumps cav-hoomd forces to NPZ in **GSD particle (tag) order** for use with the comparison script (requires cav-hoomd installed).

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
