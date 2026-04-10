# UMA Ice RPMD Simulation

Test if ice remains frozen at 243K using UMA potential with RPMD.

## Layout (repo organization)

| Path | Role |
|------|------|
| [`scripts/`](scripts/) | Primary CLIs: OpenMM reference driver, i-PI/LAMMPS driver, plotting. |
| [`ice_rpmd_workflow/`](ice_rpmd_workflow/) | Importable helpers (`ice_order_parameters`, `rpmd_thermo_utils`, LAMMPS/OpenMM geometry). Flat shims at the parent directory keep legacy `from ice_order_parameters import …` working. |
| [`fixtures/golden/`](fixtures/golden/) | Small committed reference files for regression tests (optional). |
| [`ipi/`](ipi/) | i-PI XML and inputs; trajectory `*.xyz`, `RESTART`, and logs are **gitignored**—regenerate with `scripts/run_ipi_lammps_uma_rpmd.py`. |
| [`pipeline_out/`](pipeline_out/) | Default CSV/plot output directory (**gitignored**). Create it locally when running workflows. |

Thin **shims** at the old names (`run_openmm_rpmd_reference.py`, `run_ipi_lammps_uma_rpmd.py`, `plot_rpmd_comparison.py`, …) delegate to `scripts/` so existing commands keep working.

### Pytest (repo root)

```bash
# Fast subset (no FairChem download, no CUDA-marked tests, no slow tests)
python -m pytest tests/openmmml tests/uma_ice_rpmd -m "not cuda and not slow and not fairchem" -q

# OpenMMML static / source guards only
python -m pytest tests/openmmml -m "openmmml and not fairchem" -q

# Full FairChem parity (downloads UMA checkpoint; slow)
python -m pytest tests/openmmml/test_uma_fairchem_batch_parity.py -q
```

Requires **OpenMM with PythonForce** and matching plugins for GPU and many integration tests; without them, pytest skips or fails only those modules.

### Regenerating reference outputs

After changing order-parameter code or simulation settings, re-run your drivers and replace files under `fixtures/golden/` only when you intend to update baselines. Do not commit large trajectories; keep `pipeline_out/` and `ipi/*.xyz` local.

## Overview

This simulation tests the stability of ice at sub-freezing temperatures using:
- **UMA potential**: Small model **`uma-s-1p1-pythonforce-batch`** (HF removed `uma-s-1.pt`; `uma-s-1-pythonforce-batch` is remapped to p1 in current openmmml)
- **RPMD**: Ring Polymer Molecular Dynamics (default **8 beads**; use **`--beads 1`** for classical MD and **~LAMMPS-equivalent UMA cost** per step)
- **Temperature**: 243 K (30 degrees below freezing)
- **Ice Ih structure**: proton-ordered orthorhombic supercell from `generate_ice_ih.py`, chosen by **explicit replication** `--nx`, `--ny`, `--nz` (8 molecules per cell); otherwise GenIce CIF, GenIce2, or embedded CIF fallback

**LAMMPS force parity:** OpenMM UMA uses **per-atom** PBC wrap by default (same as Fairchem LAMMPS `fix external`). Override with **`--use-molecular-wrap-for-uma`** on `test_uma_ice_rpmd.py` if you need the old molecular wrap. See [docs/LAMMPS_OPENMM_UMA_PARITY.md](docs/LAMMPS_OPENMM_UMA_PARITY.md).

### Speed vs LAMMPS (why OpenMM can feel slow)

| Cause | What to do |
|-------|------------|
| **Multiple RPMD beads** | Each MD step runs UMA on **every bead** (batched, but still ~N× GNN work vs one classical config). LAMMPS ice uses **one config per step**. Use **`--beads 1`** to match that cost when benchmarking or developing. |
| **Same ML stack** | LAMMPS (`fix external`) and OpenMM both call Fairchem **`predict(AtomicData)`**—no faster LAMMPS-only API. |
| **CUDA sync (optional)** | openmmml avoids **`torch.cuda.synchronize()`** / **`empty_cache()`** in the hot path by default. Set **`OPENMMML_CUDA_SYNC_AFTER_FORCE=1`** only if debugging OOM or async errors (hurts throughput). |
| **Profiling** | **`OPENMMML_DEBUG_LOGS=1`** plus **`--debug-logs`** on `test_uma_ice_rpmd.py` prints per-step **prep / atomicdata / inference** ms so you can see where time goes. |
| **Turbo + batching (RPMD)** | Batched RPMD already calls UMA once per step with all beads. **`--inference-turbo-rpmd`** (or `createSystem(..., inference_settings='turbo_rpmd')`) turns on **TF32**, **no activation checkpointing**, and Fairchem **`merge_mole=True`** (constant composition across beads). Full Fairchem **turbo** also sets **`torch.compile`**; that is **off** here for UMA stability. **`--optimize-inference`** is TF32 + no checkpoint but **without** merge_mole (slightly lower peak VRAM in some cases). |

Example (fast classical comparison to LAMMPS):

```bash
OPENMMML_DEBUG_LOGS=1 python test_uma_ice_rpmd.py --molecules 128 --beads 1 --platform cuda \
  --ml-device cuda --temperature 243 --pressure 0 --dt 1.0 --equil 0 --prod 0.1 \
  --model uma-s-1p1-pythonforce-batch --debug-logs --minimal
```

### Ice order vs disorder (`run_openmm_ice_lammps_match.py`)

To quantify **melting / loss of crystal order** (oxygen framework):

| Quantity | Interpretation |
|----------|----------------|
| **Q6** (Steinhardt *l*=6, mean over O) | Higher in **ordered** ice-like shells; **drops** as structure becomes liquid-like / amorphous. |
| **q_tet** (Chau–Harding, mean over O) | **1** = perfect tetrahedron of 4 nearest O neighbors; **lower** when second shell or disorder breaks tetrahedral H-bond geometry. |

Requires **scipy** (uses `sph_harm_y`). First O–O shell is included with **r** in ~(2.0, 3.7) Å (see `ice_order_parameters.py`).

```bash
python run_openmm_ice_lammps_match.py --steps 5000 --report-every 100 --order-every 100 \
  --order-csv ice_order.csv --platform cuda --ml-device cuda --model uma-s-1p1-pythonforce-batch
# Or omit --order-every: --order-csv alone defaults order sampling to --report-every (CSV + Q6/q_tet on log lines).
```

Use **`--report-every`** ≤ **`--order-every`** (or equal) so order lines are printed when you expect. Post-process **`ice_order.csv`** (e.g. plot `q6_mean` vs `time_ps`).

### OpenMM vs LAMMPS same experiment (order parameters)

To compare **ice disorder (Q6, q_tet)** between OpenMM and LAMMPS on the same protocol:

1. **Build shared structure** (once):  
   `cd tests/uma_ice_rpmd && python lammps/build_lammps_ice_data.py -n 128`

2. **OpenMM** (10k steps, dump every 10, order CSV):  
   `python run_openmm_ice_lammps_match.py --steps 10000 --traj-every 10 --report-every 10 --order-csv ice_order.csv --platform cuda --ml-device cuda --model uma-s-1p1-pythonforce-batch`

3. **LAMMPS** (same length/interval via long deck or CLI):  
   `cd lammps && python run_lammps_uma_ice.py --infile in.ice_uma_openmm_match_long.lmp --device cuda`  
   Or override from short deck:  
   `python run_lammps_uma_ice.py --infile in.ice_uma_openmm_match.lmp --steps 10000 --dump-every 10 --thermo-every 10 --device cuda`

4. **Post-process LAMMPS trajectory** → same CSV schema as OpenMM:  
   `cd .. && python lammps_order_from_dump.py --data lammps/data.ice_uma --dump lammps/dump.openmm_match_long.lammpstrj -o ice_order_lammps.csv`

5. **Compare** `ice_order.csv` (OpenMM) vs `ice_order_lammps.csv` (LAMMPS): same columns; trajectories will differ (different Langevin RNG) but **trends** (q_tet / Q6 drift) should be comparable.

### RPMD benchmark: OpenMM vs i-PI + UMA

Compare **path-integral RPMD** (243 K, PILE-G / pile_g) between OpenMM and i-PI. The i-PI side follows the standard two-process model: an **i-PI server** (ring-polymer integrator + socket) and a **force client** that returns UMA energies/forces for each bead configuration.

Two force-client backends are supported via `--client`:

| `--client` | Stack | Prerequisites |
|------------|-------|---------------|
| `lammps` (default) | LAMMPS Python API + `fix ipi` + `fix external` (FairChem UMA) | `pip install fairchem-lammps fairchem-core lammps` (LAMMPS with MISC package) |
| `python` | `i-pi-py_driver` + ASE `FAIRChemCalculator` | `pip install ipi fairchem-core` |

The pipeline (`run_rpmd_32x32.sh`) uses **64 molecules** (ice Ih 2x2x2 supercell) and **32 beads** by default. Set `RPMD_IPI_CLIENT=python` to use the ASE driver instead of LAMMPS.

**Monitoring:** i-PI progress lives in `ipi/i-pi_run.log`. Bead trajectories (`ice__i-pi.traj_*.xyz`) show how many frames have been saved (see `stride` in `ipi/input.xml`).

**Why LAMMPS can look "stuck" at Step 0:** Under `fix ipi`, dynamics are owned by **i-PI** (see LAMMPS docs: initial LAMMPS coordinates are replaced by i-PI). The warning `No fixes with time integration` is expected. See **[docs/IPI_LAMMPS_OUTPUT_AND_PROGRESS.md](docs/IPI_LAMMPS_OUTPUT_AND_PROGRESS.md)**.

#### Quick start (orchestrated)

1. **Build LAMMPS data** (64 molecules = 2x2x2 cells):
   `python lammps/build_lammps_ice_data.py --nx 2 --ny 2 --nz 2 -o lammps/data.ice_uma_64`
2. **Run i-PI + UMA**:
   `python run_ipi_lammps_uma_rpmd.py --molecules 64 --beads 4 --steps 100`
   (uses LAMMPS client by default; add `--client python` for the ASE driver)
3. **Post-process i-PI trajectory** (path-integral order over beads):
   `python ipi_order_from_traj.py --traj ipi/ice__i-pi.traj_0.xyz --beads 4 -o pipeline_out/ice_order_ipi_rpmd.csv`
4. **Run OpenMM RPMD reference**:
   `python run_openmm_rpmd_reference.py --nx 2 --ny 2 --nz 2 --beads 4 --steps 100`
5. **Plot comparison (RPMD)**:
   `python plot_rpmd_comparison.py -o pipeline_out/rpmd_comparison.png`

#### Manual two-terminal workflow (no orchestration script)

If you prefer to run i-PI and the force client yourself:

```bash
# Terminal 1: prep + i-PI server
cd tests/uma_ice_rpmd
python lammps/build_lammps_ice_data.py --nx 2 --ny 2 --nz 2 -o lammps/data.ice_uma_64
python ipi/convert_lammps_to_ipi_xyz.py --data lammps/data.ice_uma_64 -o ipi/init.xyz
cd ipi
i-pi input.xml

# Terminal 2: LAMMPS force client (after i-PI prints "Created inet socket")
cd tests/uma_ice_rpmd
python -c "
from lammps import lammps
from fairchem.lammps.lammps_fc import FIX_EXT_ID, FIX_EXTERNAL_CMD, FixExternalCallback
from fairchem.core import pretrained_mlip
lmp = lammps(cmdargs=['-nocite', '-log', 'none'])
lmp._predictor = pretrained_mlip.get_predict_unit('uma-s-1p1', device='cuda')
lmp._task_name = 'omol'
lmp.commands_list([
    'units metal', 'atom_style atomic', 'boundary p p p',
    'read_data lammps/data.ice_uma_64',
    'comm_modify cutoff 12.0',
    'pair_style zero 1.0', 'pair_coeff * *',
])
lmp.command(FIX_EXTERNAL_CMD)
lmp.set_fix_external_callback(FIX_EXT_ID, FixExternalCallback(charge=0, spin=1), lmp)
lmp.command('fix ipi_client all ipi 127.0.0.1 65535')
lmp.command('run 1000')
del lmp._predictor; lmp.close()
"
```

Edit `ipi/input.xml` to set beads, steps, timestep, thermostat, and port before starting.

#### Troubleshooting: "i-PI did not become ready in time"

If i-PI crashes during startup (before opening its socket), the orchestrator times out. The most common cause is a **broken PyFFTW** installation:

```
ValueError: cannot reshape array of size 0 into shape (32,576)
```

Fix: `pip uninstall pyfftw` -- i-PI falls back to NumPy FFT, which works correctly.

The orchestrator also prepends a small `pyfftw` stub on `PYTHONPATH` for the **i-PI subprocess** when it detects broken aligned allocation, so NumPy FFT is used without uninstalling PyFFTW.

#### Troubleshooting: `UMA driver exited with code 1` (Python client)

`fairchem.core` imports `InferenceBatcher`, which pulls **Ray Serve** and **FastAPI**. A **corrupted or mismatched** `pydantic-core` install often surfaces as:

`ImportError: cannot import name 'validate_core_schema' from 'pydantic_core'`

Repair example (sync to versions your `pydantic` package expects; check with `pip show pydantic pydantic-core`):

```bash
pip install --ignore-installed "pydantic-core==2.41.5" "pydantic==2.12.5"
```

Use `--ignore-installed` if pip reports a broken `pydantic-core` uninstall (missing `RECORD`).

**Classical MD only** (LAMMPS UMA vs OpenMM UMA vs TIP4P/2005f — **no** RPMD curves): `python plot_md_comparison.py -o pipeline_out/md_order_comparison.png`. Defaults read `ice_order_uma_lammps.csv`, `ice_order_uma_openmm.csv`, `ice_order_classical.csv` from `pipeline_out/` (from `run_ice_pipeline.py` or equivalent). **Kinetic temperature / PE on the LAMMPS curve:** `run_lammps_uma_ice.py` writes **`lammps/lammps_uma_run.log`** (LAMMPS `-log`); `lammps_order_from_dump.py --thermo-log` fills **T_K** and **PE_kj_mol** by matching **Step** to each dump frame (`run_ice_pipeline.py` passes this automatically). Use `--no-log-file` on the LAMMPS runner only if you do not need those columns.

#### 32 beads × 64 molecules (full RPMD resolution)

Pipeline default supercell is **2×2×2** (64 H₂O). The script name `run_rpmd_32x32.sh` refers to **32 RPMD beads**.

Batched UMA on **all 32 beads at once** can exceed VRAM on ~12 GB GPUs. The OpenMM path supports **chunked inference**: set `OPENMMML_UMA_RPMD_CHUNK` (e.g. `4` or `8`) so each ML call evaluates at most that many beads at a time (see `openmmml` `umapotential_pythonforce_batch.py`). Also set `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` if you see fragmentation OOM.

RPMD order CSVs (OpenMM UMA, TIP4P, i-PI post-process) use **path-integral** means:
per-oxygen ⟨Q6⟩ = mean over beads of Q6 evaluated on each bead’s positions, then
`q6_mean` / `q6_std` / `q6_p10` / `q6_p50` / `q6_p90` / `q6_min` / `q6_max` over oxygens
(same pattern for `q_tet_*`), plus `n_q6_valid`, `n_q_tet_valid`, and `n_oxygen`.
Chau–Harding *q* is **not** bounded below by 0, so `q_tet_min` / `q_tet_p10` may be
negative for distorted shells. This avoids evaluating observables on **centroid**
coordinates.

Example (OpenMM UMA RPMD, default **1 ps** = 10000 steps @ 0.1 fs in `run_rpmd_32x32.sh`; override with `RPMD_STEPS`):

```bash
cd tests/uma_ice_rpmd
export OPENMMML_UMA_RPMD_CHUNK=4
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
python run_openmm_rpmd_reference.py --nx 2 --ny 2 --nz 2 --beads 32 --dt 0.1 --steps 10000 --platform cuda
```

`run_openmm_rpmd_reference.py` passes **`--equil 2.0`** (2 ps ring-polymer thermostat equilibration) by default; use **`--equil 0`** for quick smoke tests. Velocities are drawn **after** energy minimization in `test_uma_ice_rpmd.py` (consistent with TIP4P RPMD). Optional **`--nve-check PS`** on `test_uma_ice_rpmd.py` runs a short NVE segment first and prints energy-drift diagnostics.

End-to-end (**OpenMM UMA RPMD PILE-G → OpenMM UMA RPMD PILE-G → i-PI `pile_g` + LAMMPS UMA → order CSV → 3-way plot**):

```bash
bash run_rpmd_32x32.sh
```

The script **best-effort kills** prior `run_openmm_*`, `run_ipi_*`, `test_uma_ice_rpmd`, and `i-pi` processes, then runs both models (LAMMPS client by default; override with `RPMD_IPI_CLIENT=python`) with matched `RPMD_STEPS` / `RPMD_DT_FS` (default **10000** steps @ 0.1 fs → **1.0 ps**). **`RPMD_EQUIL_PS` defaults to 2** (OpenMM NVT equilibration before production order rows, aligned with `run_openmm_rpmd_reference.py`); use **`RPMD_EQUIL_PS=0`** only for quick smoke tests (stderr warning). Outputs include `pipeline_out/ice_order_tip4p_rpmd.csv`, `ice_order_openmm_rpmd.csv`, `ice_order_ipi_rpmd.csv`, and `pipeline_out/rpmd_comparison_32x32.png`.

**OpenMM UMA RPMD only** (skip TIP4P and i-PI): `bash run_openmm_uma_rpmd_only.sh` — writes `pipeline_out/ice_order_openmm_rpmd.csv` with the same defaults as step 2 of `run_rpmd_32x32.sh`. Optional: `RPMD_STEPS`, `RPMD_DT_FS`, `OPENMMML_UMA_RPMD_CHUNK`.

**UMA vs TIP4P only** (no i-PI CSV required): `plot_rpmd_comparison.py` accepts `--openmm` and `--tip4p` and omits `--ipi` if you only want the two OpenMM curves, e.g. `pipeline_out/uma_vs_tip4p_rpmd_order.png`. The TIP4P runner (`run_openmm_tip4p_rpmd.py`) fills `T_K` / `PE_kj_mol` on each order row so the temperature panel compares both models. Match `--beads` to your RPMD count if it is not 32.

i-PI 3 writes bead trajectories as `ipi/ice__i-pi.traj_00.xyz` … `traj_31.xyz`. Pass `--traj ipi/ice__i-pi.traj_00.xyz` to `ipi_order_from_traj.py` (zero-padded indices are detected automatically).

### End-to-end pipeline (1 ps, three models)

A **single script** runs the full workflow sequentially (no parallel simulations, to avoid GPU OOM):

1. **Build** proton-ordered ice **Ih** LAMMPS data from **`--nx/--ny/--nz`** (default **2×2×2** → 64 H₂O, same as the RPMD benchmark). Output: `lammps/data.ice_uma_<8·nx·ny·nz>`.
2. **Run** 1 ps NVT (Langevin τ = 0.1 ps, **0.1 fs** default for TIP4P/2005f, 243 K) with **TIP4P/2005f** (OpenMM).
3. **Run** 1 ps same conditions with **UMA** (OpenMM).
4. **Run** 1 ps same conditions with **UMA** (LAMMPS), then post-process dump → order CSV (with **`--thermo-log`** for **T_K** / **PE**).
5. **Plot** Q6 and q_tet vs time for all three runs → `order_vs_time.png`.

All three simulations use the **same initial structure**, **same box**, **same thermostat and timestep** (default **0.1 fs** for TIP4P/2005f stability). Run length is 1 ps (10000 steps @ 0.1 fs). Order parameters and trajectory frames are written every 10 fs.

**Run (foreground):**
```bash
cd tests/uma_ice_rpmd
python run_ice_pipeline.py --out-dir pipeline_out
# Explicit supercell (example: 1×2×2 → 32 H₂O):
# python run_ice_pipeline.py --out-dir pipeline_out --nx 1 --ny 2 --nz 2
```

**Run in background** (one command; log to `pipeline.log`):
```bash
cd tests/uma_ice_rpmd
nohup python run_ice_pipeline.py --out-dir pipeline_out >> pipeline.log 2>&1 &
# Monitor: tail -f pipeline.log
```

**Classical (TIP4P/2005f) stability:** The pipeline uses **0.1 fs** by default (paper-recommended). The TIP4P/2005f force field matches González & Abascal (J. Chem. Phys. 135, 224516, 2011), Table I and Eq. (4) (M-site d_OM^rel = 0.13194). The script applies **geometry relaxation** by default so minimization is stable. The pipeline uses **GPU (cuda)** for all runs by default. If the classical CSV still shows huge T_K, use **`--variable-step-classical`** so the classical step uses OpenMM’s **VariableLangevinIntegrator** (adaptive dt). The classical script asserts **PBC**: nonbonded cutoff (0.7 nm) &lt; half the smallest box length, and aborts if max force after minimization exceeds 1e5 kJ/(mol·nm).

```bash
nohup python run_ice_pipeline.py --out-dir pipeline_out --variable-step-classical >> pipeline.log 2>&1 &
```

Options: `--nx` / `--ny` / `--nz` (default **2** each), `--work-dir` (default: script directory), `--dt-fs` (default **0.1** fs for TIP4P/2005f stability), `--steps` (default: 10000 for 1 ps @ 0.1 fs), `--skip-lammps`, `--only-plot`, `--platform-classical` (default **cuda**), `--variable-step-classical`. Tests: `pytest test_ice_pipeline.py -v`.

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

### Per-atom vs molecular wrap (UMA + RPMD)

OpenMM often keeps **unwrapped** Cartesian coordinates. **Per-atom** `wrap_positions` can place O and H on
different images → long apparent O–H bonds → bad UMA graph. **`_wrap_water_molecules_per_bead`** (molecular)
avoids that by shifting each whole water into the cell with Hs by MIC from O.

**Current default:** `test_uma_ice_rpmd.py` uses **`use_atom_wrap_for_lammps_parity=True`** (per-atom, LAMMPS /
Fairchem parity). If you see runaway centroid `T_K` or obvious geometry artifacts, try
**`--use-molecular-wrap-for-uma`**. Order-parameter CSVs use **per-atom** `wrap_cartesian_orthorhombic` on
positions before Q6 (aligned with [`ipi_order_from_traj.py`](ipi_order_from_traj.py)).

Reinstall openmmml models after editing sources:

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
# Writes full trajectory by default: ice_after_openmm_lammps_match.dcd (every step). Ovito: open DCD + load topology from PDB.
# Optional multi-frame XYZ: --xyz traj.xyz   |  Speed (printed at end): ~0.2–1 ns/day typical for 128 H2O UMA on one GPU.
# Fast MD-only, no trajectory: --no-trajectory
python run_openmm_ice_lammps_match.py --skip-minimize --steps 100 --no-trajectory ...

**LAMMPS same settings:** `lammps/in.ice_uma_openmm_match.lmp` (+ optional `_notraj` for speed). See `lammps/README.md` § OpenMM-matched deck.
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

### Classical flexible water (TIP4P/2005f)

**`run_openmm_ice_classical_flex.py`** runs the same ICE setup (structure from `lammps/data.ice_uma`, 243 K, NVT Langevin) with the **TIP4P/2005f** flexible water model (González & Abascal, J. Chem. Phys. 135, 224516, 2011). Use it to compare classical vs UMA ice stability and order parameters (Q6, q_tet).

- **Force field**: TIP4P/2005f (Morse O–H, harmonic H–O–H, LJ + M-site). Parameters match the paper (SklogWiki / González & Abascal 2011), Table I and Eq. (4): M-site d_OM^rel = 0.13194 (average3 weights 0.73612, 0.13194, 0.13194). The XML is in `wrappers/python/openmm/app/data/tip4p2005f.xml` when building OpenMM from this repo; otherwise copy it or pass `--force-field /path/to/tip4p2005f.xml`.
- **Geometry relaxation**: By default the script sets each water to TIP4P/2005f equilibrium (r_OH = 0.09419 nm, H–O–H = 107.4°) before minimization, because the LAMMPS/CIF ice structure has rigid TIP4P/2005 geometry (r_OH ≈ 0.0957 nm, angle ≈ 104.5°). Use `--no-relax-geometry` only if your structure is already equilibrated with TIP4P/2005f.
- **Timestep**: 0.1 fs default (recommended for TIP4P/2005f; 0.2 fs may be unstable).
- **Output**: Same order-parameter CSV schema as `run_openmm_ice_lammps_match.py` for direct comparison.

```bash
python lammps/build_lammps_ice_data.py -n 128
python run_openmm_ice_classical_flex.py --steps 5000 --order-csv ice_order_classical.csv --platform cuda
# Compare with UMA: ice_order.csv (UMA) vs ice_order_classical.csv (TIP4P/2005f)
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
- OpenMM RPMD uses **PILE_G** by default (`--rpmd-thermostat pile-g`); use `--rpmd-thermostat pile` for PILE on all modes (older default)

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
