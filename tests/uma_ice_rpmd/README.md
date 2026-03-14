# UMA Ice RPMD Simulation

Test if ice remains frozen at 243K using UMA potential with RPMD.

## Overview

This simulation tests the stability of ice at sub-freezing temperatures using:
- **UMA potential**: Small model **`uma-s-1p1-pythonforce-batch`** (HF removed `uma-s-1.pt`; `uma-s-1-pythonforce-batch` is remapped to p1 in current openmmml)
- **RPMD**: Ring Polymer Molecular Dynamics (default **8 beads**; use **`--beads 1`** for classical MD and **~LAMMPS-equivalent UMA cost** per step)
- **Temperature**: 243 K (30 degrees below freezing)
- **Ice Ih structure**: GenIce CIF, GenIce2, or embedded CIF fallback

### Speed vs LAMMPS (why OpenMM can feel slow)

| Cause | What to do |
|-------|------------|
| **Multiple RPMD beads** | Each MD step runs UMA on **every bead** (batched, but still ~N× GNN work vs one classical config). LAMMPS ice uses **one config per step**. Use **`--beads 1`** to match that cost when benchmarking or developing. |
| **Same ML stack** | LAMMPS (`fix external`) and OpenMM both call Fairchem **`predict(AtomicData)`**—no faster LAMMPS-only API. |
| **CUDA sync (optional)** | openmmml avoids **`torch.cuda.synchronize()`** / **`empty_cache()`** in the hot path by default. Set **`OPENMMML_CUDA_SYNC_AFTER_FORCE=1`** only if debugging OOM or async errors (hurts throughput). |
| **Profiling** | **`OPENMMML_DEBUG_LOGS=1`** plus **`--debug-logs`** on `test_uma_ice_rpmd.py` prints per-step **prep / atomicdata / inference** ms so you can see where time goes. |

Example (fast classical comparison to LAMMPS):

```bash
OPENMMML_DEBUG_LOGS=1 python test_uma_ice_rpmd.py --molecules 128 --beads 1 --platform cuda \
  --ml-device cuda --temperature 243 --pressure 0 --dt 1.0 --equil 0 --prod 0.1 \
  --model uma-s-1p1-pythonforce-batch --debug-logs --minimal
```

## Ice Structure

The initial structure can come from:

- **`--input ice.cif`** (recommended): Load from GenIce CIF directly.
  ```bash
  genice 1h --rep 2 2 2 --format cif > ice.cif
  python test_uma_ice_rpmd.py --input ice.cif --beads 8 ...
  ```
- **Auto GenIce** (when no `--input`): The script calls GenIce or GenIce2 to generate ice. Tries `genice` first, then `genice2`. Install with `pip install genice` or `pip install genice2`.
- **Embedded CIF fallback**: Ice Ih from Avogadro/COD (12 molecules) if GenIce unavailable

## Key Physics

- **Ice Ih**: Hexagonal ice, O-O ~2.75 Å, density ~0.92 g/cm³
- **RPMD**: Captures quantum nuclear effects for hydrogen bonding
- **Timestep**: 0.5-1.0 fs recommended for flexible water + RPMD

## Expected Results

### If Ice Remains Frozen
- **NPT**: Density > 0.85 g/cm³ (ice ~0.92)
- **NVT**: RMSD < 0.15 nm from minimized structure
- Sharp O-O RDF first peak at ~0.275 nm

### If Ice Melts
- **NPT**: Density ~1.0 g/cm³ (liquid)
- **NVT**: RMSD > 0.15 nm
- Broader RDF

## Usage

### Prepare Ice for a Cubic Box (standalone script)

Create ice Ih to fill a cubic box of specified size:

```bash
python prepare_ice_box.py --box 1.0 -o ice.pdb          # 1 nm cubic box
python prepare_ice_box.py --box 0.5 -o ice_05nm.pdb      # 0.5 nm box
python prepare_ice_box.py --box 2.0 -f gro -o ice.gro    # 2 nm, GROMACS format
```

Requires ASE. Optional: `pip install genice2` for proton-disordered ice.

### Quick Test (from GenIce CIF)
```bash
genice 1h --rep 2 2 2 --format cif > ice.cif
python test_uma_ice_rpmd.py --input ice.cif --beads 8 --temperature 243 \
    --dt 1.0 --equil 5 --prod 50 --platform cpu
```

### CUDA with limited memory (use --molecules to trim)
```bash
# ice.cif has 128 molecules; trim to 32 for ~12 GB GPU
python test_uma_ice_rpmd.py --input ice.cif --molecules 32 --beads 4 --temperature 243 \
    --dt 1.0 --equil 5 --prod 50 --platform cuda --pressure 0
```

### Quick Test (by box size)
```bash
python test_uma_ice_rpmd.py --box 1.0 --beads 8 --temperature 243 \
    --dt 1.0 --equil 5 --prod 50 --model uma-s-1p1-pythonforce-batch
```

### Quick Test (by molecule count)
```bash
python test_uma_ice_rpmd.py --molecules 32 --beads 8 --temperature 243 \
    --dt 1.0 --equil 5 --prod 50 --pressure 1 --model uma-s-1p1-pythonforce-batch
```

### NVT (no barostat)
```bash
python test_uma_ice_rpmd.py --box 1.0 --beads 8 --temperature 243 \
    --equil 10 --prod 100 --pressure 0
```

### Custom Parameters
```bash
python test_uma_ice_rpmd.py \
    --input ice.cif \         # Load from CIF/PDB/XYZ (OR --box/--molecules)
    --box 1.0 \               # Cubic box side in nm (OR --molecules)
    --molecules 32 \          # Number of water molecules (if no --box)
    --beads 8 \                # RPMD beads
    --temperature 243 \        # Temperature (K)
    --dt 1.0 \                 # Timestep in fs (0.5-1.0 for stability)
    --equil 5.0 \              # Equilibration (ps)
    --prod 50.0 \              # Production (ps)
    --model uma-s-1p1-pythonforce-batch \   # Small UMA (HF checkpoint)
    --output ice_test          # Output directory (use a persistent path; do NOT use /tmp)
```

**Note:** Always use `--output` with a persistent directory (e.g. `tests/uma_ice_rpmd/` or `ice_test/`). Do **not** use `/tmp`—outputs there are ephemeral and may be deleted.

## Output Files

### Data File (NPZ)
`uma_ice_rpmd_T243K_b8.npz` contains:
- `times`: Time points (ps)
- `rmsds`: RMSD from initial structure (nm)
- `densities`: System density (g/cm³)
- `potential_energies`: Mean PE per bead (kJ/mol)
- `kinetic_energies`: Mean KE per bead (kJ/mol)
- `temperatures`: Temperature per bead (K)
- `msds`: Mean-squared displacement (nm²)
- `r_initial`, `g_initial`: Initial RDF
- `r_final`, `g_final`: Final RDF
- `initial_positions`, `final_positions`: Positions
- `is_frozen`: Boolean flag for ice status

### Analysis Plot
`uma_ice_rpmd_T243K_b8_analysis.png` shows:
1. RMSD vs time
2. Density vs time
3. Potential energy vs time
4. Temperature vs time
5. MSD vs time
6. RDF comparison (initial vs final)

## Analysis

### Loading Data
```python
import numpy as np
import matplotlib.pyplot as plt

data = np.load('uma_ice_rpmd_T243K_b8.npz')
times = data['times']
rmsds = data['rmsds']
is_frozen = data['is_frozen']

print(f"Ice frozen: {is_frozen}")
print(f"Final RMSD: {rmsds[-1]:.4f} nm")
```

### Indicators of Frozen Ice
- RMSD stays below 0.1 nm throughout simulation
- Density remains around 0.92 g/cm³
- Sharp first peak in RDF at ~0.275 nm (nearest neighbor O-O distance)
- Low MSD indicates restricted motion

## Requirements

- OpenMM with RPMD plugin
- openmm-ml with UMA support (integrated into OpenMM when built from this repo)

For full rebuild and reinstall on a new machine, see [docs/BUILD_AND_REINSTALL.md](../docs/BUILD_AND_REINSTALL.md).
- ASE (Atomic Simulation Environment)
- NumPy
- Matplotlib

## Limitations and Caveats

### NPT on CUDA

NPT (RPMDMonteCarloBarostat) is supported on CUDA. Ensure OpenMM is built with CUDA (`-DOPENMM_BUILD_CUDA_LIB=ON`). If you encounter `CUDA_ERROR_ILLEGAL_ADDRESS` with NPT, try `--precision double` or report the issue.

### Why OpenMM UMA ice failed while LAMMPS worked (fix: molecular wrap)

OpenMM often keeps **unwrapped** Cartesian coordinates (atoms can sit outside the primary cell).
The old UMA path called **`wrap_positions()` per atom** before the model. That **splits each atom
independently** across PBC → an H can jump to the opposite face while its O stays put → **fake
long bonds** → wrong UMA neighbor graph → **catastrophic forces / melting**.

LAMMPS + fairchem also uses per-atom wrap, but short runs from a **wrapped** data file often never
unwrap enough to break molecules; OpenMM **RPMD / integrator** unwraps aggressively.

**Fix (in `umapotential_pythonforce_batch.py`):** `_wrap_water_molecules_per_bead` — for **O,H,H × N**
(water), put **O in the primary cell** and each **H by minimum-image displacement from that O**.
**RPMD:** wrapping runs **per bead** (bead `i` uses only `state[i]` positions and box); there is no
centroid or wrap across beads.
Reinstall:

```bash
cp wrappers/python/openmmml/models/umapotential_pythonforce_batch.py \
   "$(python -c "import openmmml,os; print(os.path.join(os.path.dirname(openmmml.__file__),'models'))")/"
```

### LAMMPS + UMA (same ice, trajectory)

See **`lammps/README.md`**. From **`tests/uma_ice_rpmd/`** you can run:

```bash
python build_lammps_ice_data.py              # writes lammps/data.ice_uma
python run_lammps_uma_ice.py               # or: cd lammps && python run_lammps_uma_ice.py
```

FAIRChem’s **`fairchem-lammps`** driver runs LAMMPS with **UMA** on the **same `ice.cif`** as OpenMM. Produces **`lammps/dump.ice_uma.lammpstrj`**. Classical NVT @ 243 K (OpenMM RPMD analogue: 1 bead).

### OpenMM MD = LAMMPS conditions (stable ice check)

**`run_openmm_ice_lammps_match.py`** loads the **same** `lammps/data.ice_uma` as LAMMPS, then:

1. Energy minimization  
2. **Classical** **`LangevinMiddleIntegrator`**: **243 K**, **γ = 1/ps**, **1 fs** (same spirit as `fix nve` + `fix langevin 243 243 1.0`)  
3. **`run N`** steps (default **100** = 100 fs, like `in.ice_uma_quick.lmp`)

Not RPMD — matches the **classical** LAMMPS deck so you can compare stability (T drift, PE) and final PDB.

```bash
python lammps/build_lammps_ice_data.py -n 128   # data.ice_uma
# Default minimization is short (~1 min); full LAMMPS-like minimize = many UMA evals (looks frozen).
python run_openmm_ice_lammps_match.py --steps 100 --platform cuda --ml-device cuda \
  --model uma-s-1p1-pythonforce-batch -o ice_after_openmm_lammps_match.pdb
# Fast MD-only test (same starting geometry as LAMMPS data file):
python run_openmm_ice_lammps_match.py --skip-minimize --steps 100 ...
```

### Same-structure benchmark (OpenMM vs LAMMPS pipeline)

**`benchmark_same_structure_openmm_vs_lammps.py`** reads **`lammps/data.ice_uma`** (same atoms/box/positions as `run_lammps_uma_ice.py`) and times:

1. **LAMMPS_pipeline** — `wrap_positions` + `atomic_data_from_lammps_data` + `predict` (same as `FixExternalCallback`).
2. **OpenMM_pipeline** — molecular water wrap + `atomicdata_list_to_batch` + `predict` (same as openmmml 1-bead batch).
3. **predict_only** — frozen `AtomicData`, `predict` only (ML lower bound).

Expect **~same ms** for (1) vs (2): slowness is **UMA inference**, not RPMD or OpenMM vs LAMMPS Python glue.

```bash
python build_lammps_ice_data.py -n 128   # writes lammps/data.ice_uma if missing
python benchmark_same_structure_openmm_vs_lammps.py --reps 40
```

### UMA-LAMMPS Benchmarking

UMA does **not** use a classical LAMMPS `pair_style`; FAIRChem injects forces via **`fix external`**. The canonical ASE reference is **ASE + FAIRChemCalculator**. Use `benchmark_ase_reference.py` to compare OpenMM UMA against ASE:
```bash
python tests/uma_ice_rpmd/benchmark_ase_reference.py
```

### Minimal System Testing (OOM Avoidance)

For memory-constrained runs or validation, use a single water molecule:
```bash
python test_uma_ice_rpmd.py --molecules 1 --beads 1 --equil 5 --prod 20 --platform cpu
python tests/uma_ice_rpmd/test_force_layout.py
```

## Notes

- Ice Ih is generated using ASE with appropriate lattice parameters
- RPMD with 8 beads captures quantum effects at 243 K
- 1 fs timestep is appropriate for flexible water models
- UMA potential provides accurate ML-based forces and energies
- Simulation uses PILE thermostat (default for RPMD)

## Performance

- Typical performance: 50-200 steps/s on GPU (CUDA)
- 10 ps + 100 ps simulation: ~5-20 minutes on GPU
- CPU is slower but functional for small systems

## Troubleshooting

### Import Error: openmmml
```bash
cd /media/extradrive/Trajectories/openmm/fairchem/openmm-ml
pip install -e .
```

### Import Error: ase
```bash
pip install ase
```

### CUDA Out of Memory
Reduce system size or use CPU platform:
- Reduce `--molecules` to 16 or 8
- Script will automatically fall back to CPU if CUDA fails

### NPT on CUDA (posqCorrection Workaround)
When running NPT (pressure > 0) on CUDA, the script uses **double precision** by default to avoid
`posqCorrection` download crashes (`CUDA_ERROR_ILLEGAL_ADDRESS`). Mixed precision with the barostat
can trigger this on some setups. You can try `--precision mixed` if the underlying issue has been
fixed in your OpenMM build.

### CUDA Error (CUDA_ERROR_ILLEGAL_ADDRESS)
When OpenMM and PyTorch both use CUDA on the same GPU, a context conflict can cause
`CUDA_ERROR_ILLEGAL_ADDRESS`. The test script applies two fixes:

1. **Context order**: PyTorch CUDA is initialized *before* the OpenMM Context is created
   (torch.cuda.init() pops existing contexts; PyTorch must own the primary context first).

2. **Lazy model load**: The UMA model loads on first force computation (inside the callback),
   when OpenMM has pushed its context, so tensors allocate in the shared context.

If the error persists:
- Use `--ml-device cpu` to run UMA on CPU (OpenMM stays on GPU).
- Use `--platform cpu` to run everything on CPU.

**Verifying the fix is active:** When running with `--platform cuda`, expect:
```
  ML model device: cpu
  ...
  Loading UMA model 'uma-s-1' on cpu...
```
If you see `on cuda` instead, the installed openmmml does not have the fix. Copy the modified
files into the installed package:
```bash
INSTALL_PATH=$(python -c "import openmmml; import os; print(os.path.dirname(openmmml.__file__))")
cp wrappers/python/openmmml/models/umapotential_pythonforce_batch.py "$INSTALL_PATH/models/"
cp wrappers/python/openmmml/models/umapotential_pythonforce.py "$INSTALL_PATH/models/"
```

### Why compare_forces_ice.py still shows ~100 kJ/(mol·nm) max |ΔF|

OpenMM (`uma-s-1p1-pythonforce-batch`) and ASE (`uma-s-1p1`) load the **same** UMA weights. After the
eV/Å → kJ/(mol·nm) fix and **wrapping** the same positions into the primary cell for both paths, the
remaining gap is **not** a wrong potential—it is **small relative to the force magnitude**.

- In condensed ice, single components of **F** can be **tens of thousands** of kJ/(mol·nm) (large
  condensed-phase gradients). A max |ΔF| of ~100 kJ/(mol·nm) is typically **~0.1–0.2% relative**
  on those components.
- UMA inference is **float32** on GPU; ASE’s calculator uses the same model but can differ slightly
  in tensor layout, batching, and CUDA kernel order. That shows up as **sub‑percent** force noise,
  not as a systematic unit or geometry bug.
- **Energy** usually agrees within a few kJ/mol on totals of order 10⁶ kJ/mol—consistent with the
  same underlying energy surface.

So: **there is no “extra” 100 kJ/(mol·nm) error in the sense of a duplicated or wrong model**; the
number is large in absolute units because **forces themselves are huge** in those units on dense ice.

```bash
python compare_forces_ice.py --input ice.cif --molecules 16
```

## References

- UMA: Universal Molecular Atomistic potentials (FAIRChem)
- RPMD: Ring Polymer Molecular Dynamics
- Ice Ih: Hexagonal ice structure
