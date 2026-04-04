# LAMMPS vs OpenMM UMA parity (per-atom PBC wrap)

## Convention

OpenMM ice / benchmark scripts that target **LAMMPS `fix external` + Fairchem** use **`use_atom_wrap_for_lammps_parity=True`** on `MLPotential.createSystem`. That enables ASE `wrap_positions` on each atomic configuration (each RPMD bead independently in the batched path), matching `fairchem.lammps` preprocessing.

The default in [`umapotential_pythonforce_batch.py`](../../../wrappers/python/openmmml/models/umapotential_pythonforce_batch.py) remains `False` for backward compatibility outside this repo’s parity workflows.

## Main simulation entry point

[`test_uma_ice_rpmd.py`](../test_uma_ice_rpmd.py) defaults to **per-atom wrap** (`use_atom_wrap_for_lammps_parity=True`). Use **`--use-molecular-wrap-for-uma`** for molecular water wrap (`_wrap_water_molecules_per_bead`).

## i-PI / plotting parity

[`ipi_order_from_traj.py`](../ipi_order_from_traj.py) defaults to the same **per-atom** Cartesian wrap (`wrap_cartesian_orthorhombic`) before Q6 / q_tet. Legacy molecular grouping: **`--molecular-wrap`**. Full CLI and OpenMM vs i-PI ledger: [`RPMD_PARITY_PROTOCOL_LEDGER.md`](RPMD_PARITY_PROTOCOL_LEDGER.md).

## Verification commands

1. **Static forces** (same `data` file as LAMMPS):

   ```bash
   cd tests/uma_ice_rpmd
   python compare_forces_openmm_vs_lammps.py --data lammps/data.ice_uma_64
   ```

   If OpenMM raises `CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)`, your GPU driver is older than the CUDA/PTX used to build OpenMM. Either upgrade the NVIDIA driver, rebuild OpenMM against an older toolkit, or force CPU integration:

   `python compare_forces_openmm_vs_lammps.py --data ... --openmm-platform cpu --device cpu`

2. **Pytest regression** (small system, CPU LAMMPS-style path):

   ```bash
   pytest tests/uma_ice_rpmd/test_openmm_lammps_uma_force_parity.py -v
   ```

3. **NVE continuation from a LAMMPS dump** (OpenMM already uses atom wrap):

   ```bash
   python run_openmm_from_lammps_state.py --state-dump ... --data lammps/data.ice_uma ...
   ```

4. **NVT checklist** (matched seeds, dt, Langevin): see [`lammps/README.md`](../lammps/README.md) and [`run_openmm_ice_lammps_match.py`](../run_openmm_ice_lammps_match.py).

## NVE / thermostat

Parity of **trajectories** still requires identical **initial velocities**, **timestep**, **thermostat parameters**, and **minimize** protocol—not only matching UMA forces.
