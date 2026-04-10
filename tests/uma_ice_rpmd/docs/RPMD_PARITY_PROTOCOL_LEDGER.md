# RPMD parity protocol ledger (OpenMM UMA vs i-PI + LAMMPS UMA)

This file records **CLI and conventions** used for `pipeline_out/rpmd_comparison.png` and related CSVs so dynamics/postprocessing can be reproduced and Phase-2 debugging stays grounded.

## PBC / postprocessing (Phase 1)

| Item | Convention |
|------|------------|
| i-PI order CSV | [`ipi_order_from_traj.py`](../ipi_order_from_traj.py) default **`--pbc-wrap` implicit = per-atom** (`wrap_cartesian_orthorhombic`), matching LAMMPS/Fairchem `wrap_positions`. Legacy: `--molecular-wrap`. |
| OpenMM order CSV | [`test_uma_ice_rpmd.py`](../test_uma_ice_rpmd.py) applies **`wrap_cartesian_orthorhombic`** to each bead’s positions **before** `ice_order_metrics_path_integral` (analysis-only; same as i-PI export). |
| UMA force path | Default **`use_atom_wrap_for_lammps_parity=True`** (omit `--use-molecular-wrap-for-uma`). |

## OpenMM UMA RPMD (typical `pipeline_out/ice_order_openmm_rpmd.csv`)

- **Driver:** `run_openmm_rpmd_reference.py` (invoked from `run_openmm_uma_rpmd_only.sh` or `test_uma_ice_rpmd.py`).
- **Shell defaults** ([`run_openmm_uma_rpmd_only.sh`](../run_openmm_uma_rpmd_only.sh), [`run_rpmd_32x32.sh`](../run_rpmd_32x32.sh)): `RPMD_STEPS=10000`, `RPMD_DT_FS=0.1`, **`RPMD_EQUIL_PS=2`** (ps NVT equil before production order rows), `OPENMMML_UMA_RPMD_CHUNK=4`, 32 beads, PILE-G, `nx ny nz = 2 2 2` (64 H₂O). Use **`RPMD_EQUIL_PS=0`** only for smoke tests (stderr warning in scripts).
- **Reference driver** ([`run_openmm_rpmd_reference.py`](../run_openmm_rpmd_reference.py)) also defaults **`--equil 2.0`** ps; **confirm** overrides when auditing a CSV.
- **Order CSV columns:** path-integral ⟨Q6⟩, ⟨q_tet⟩ over beads; **`T_K`** = centroid kinetic temperature (see [`KINETIC_TEMPERATURE_UMA_RPMD.md`](../KINETIC_TEMPERATURE_UMA_RPMD.md)).
- **Seed:** CLI default `284759` in `test_uma_ice_rpmd.py` / reference script (verify `--seed` for the run you plot).

## i-PI + LAMMPS UMA (`pipeline_out/ice_order_ipi_rpmd.csv`)

- **Orchestrator:** [`run_ipi_lammps_uma_rpmd.py`](../run_ipi_lammps_uma_rpmd.py).
- **i-PI input:** [`ipi/input.xml`](../ipi/input.xml) — thermostat `pile_g`, **`tau`** (fs), ring polymer **beads**, **timestep**, **seed** must match the intended comparison.
- **Trajectory stride:** `IPI_TRAJ_OUTPUT_STRIDE` in the orchestrator (must match `input.xml` trajectory/properties stride); `time_ps` in the order CSV follows **`Step:` × dt / 1000** when `Step:` is present in xyz comments.
- **`T_K` / `PE_kj_mol`:** merged from `ipi/ice__i-pi.md` via [`ipi_thermo_utils.py`](../ipi_thermo_utils.py) — **`temperature(nm=0)`** is the centroid-mode kinetic temperature line to compare against OpenMM **`T_K`** (see [`IPI_OPENMM_THERMOSTAT_PARITY.md`](IPI_OPENMM_THERMOSTAT_PARITY.md)).
- **Init:** LAMMPS `data` → `convert_lammps_to_ipi_xyz.py` → `ipi/init.xyz`.

## Phase 2 experiments (forces, thermostat, wrap A/B)

### Minimal OpenMM RPMD smoke (5 integration steps)

Use this to verify the stack (PythonForce + UMA + RPMD integrator) without a long run:

```bash
cd tests/uma_ice_rpmd
export OPENMM_PLUGIN_DIR=/path/to/your/openmm/build/lib/plugins   # required for uma-s-1p1-pythonforce-batch
export OPENMMML_UMA_RPMD_CHUNK=4
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
python run_openmm_rpmd_reference.py \
  --nx 2 --ny 2 --nz 2 --beads 4 --dt 0.1 --equil 0 --steps 5 \
  --order-every 1 --report-every-steps 1 \
  --optimize-inference-tf32-only \
  --order-csv pipeline_out/phase2_openmm_5step.csv \
  --platform cuda
```

Conda OpenMM without the Python ML plugin fails with `AttributeError: module 'openmm' has no attribute 'PythonForce'`.

### Longer checklist

1. **Static forces:**  
   `python compare_forces_openmm_vs_lammps.py --data lammps/data.ice_uma_64`  
   Use `--openmm-platform cpu` if CUDA/PTX mismatches. Requires an OpenMM build that exposes **`PythonForce`** (custom OpenMM from this tree + `OPENMM_PLUGIN_DIR`); stock conda OpenMM without the plugin will fail at `createSystem`.

2. **Forces along trajectories:**  
   `python compare_forces_during_dynamics.py --dump lammps/dump.pipeline_1ps_02fs.lammpstrj --data lammps/data.ice_uma_64 --frames 0 1 --device cpu`  
   Use a dump whose atom count matches `--data` (e.g. PIMD `lammps/pimd_output/dump.bead_000.lammpstrj` with `data.ice_uma_8`).  
   Fairchem on **PyTorch 2.4+** may raise `Could not infer dtype of numpy.bool_` inside UMA MoE heads (`escn_moe.py`); upgrade **`fairchem-core`** to a release that fixes numpy-boolean indexing, or use a slightly older PyTorch for these diagnostics.

3. **Thermostat mapping:** short jobs varying i-PI `--ipi-tau-fs` and OpenMM `--rpmd-friction` / `--rpmd-centroid-friction`; compare i-PI `temperature(nm=0)` vs OpenMM `T_K` on the same time base after postprocessing alignment.

4. **OpenMM wrap A/B (dynamics):** two runs with the same seed and `--equil`, with vs without `--use-molecular-wrap-for-uma`; compare `ice_order_*_rpmd.csv` drift (document bead count, steps, and equilibration).
