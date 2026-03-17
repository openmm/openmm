#!/usr/bin/env python3
"""
Diagnose UMA energy discrepancy between OpenMM PythonForce and standalone computation.

Compares three energy computations at the same positions:
1. OpenMM PE via context.getState(getEnergy=True)  -- through PythonForce callback
2. Standalone direct AtomicData construction + predict_unit.predict()
3. Standalone LAMMPS-style atomic_data_from_lammps_data + predict_unit.predict()

Also checks: system force count, predict_unit.direct_forces flag, raw positions.
Reference: LAMMPS PE = -266277.22 eV at step 5000.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import torch

_SCRIPT_DIR = Path(__file__).resolve().parent
_LAMMPS_DIR = _SCRIPT_DIR / "lammps"

if os.getenv("OPENMM_PLUGIN_DIR") is None:
    for _base in [os.getenv("CONDA_PREFIX"), os.path.expanduser("~/miniconda3")]:
        if _base:
            _plugins = os.path.join(_base, "lib", "plugins")
            if os.path.isdir(_plugins):
                os.environ["OPENMM_PLUGIN_DIR"] = _plugins
                break

import openmm
from openmm import Context, CustomIntegrator, Platform, Vec3, unit
from openmm.app import Topology, element
from openmmml import MLPotential

sys.path.insert(0, str(_SCRIPT_DIR))
from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data
from run_openmm_from_lammps_state import read_lammps_state_dump, water_topology


LAMMPS_PE_EV = -266277.22
STATE_DUMP = _LAMMPS_DIR / "dump.nve_state_5000.lammpstrj"
DATA_FILE = _SCRIPT_DIR / "pipeline_out" / "shared_data.ice_uma"


def main():
    print("=" * 72)
    print("DIAGNOSTIC: UMA Energy Comparison (OpenMM vs Standalone vs LAMMPS)")
    print("=" * 72)

    _, cell, Z_data = parse_lammps_data(DATA_FILE)
    pos_ang, vel_ang_ps, types = read_lammps_state_dump(STATE_DUMP)
    n = pos_ang.shape[0]
    n_mol = n // 3
    Z = np.array([8 if t == 1 else 1 for t in types], dtype=np.int64)

    topology = water_topology(n_mol)
    positions_nm = pos_ang / 10.0
    pos_list = [
        Vec3(float(positions_nm[i, 0]), float(positions_nm[i, 1]), float(positions_nm[i, 2])) * unit.nanometer
        for i in range(n)
    ]
    vel_nm_ps = vel_ang_ps / 10.0
    vel_list = [
        Vec3(float(vel_nm_ps[i, 0]), float(vel_nm_ps[i, 1]), float(vel_nm_ps[i, 2])) * unit.nanometers / unit.picosecond
        for i in range(n)
    ]

    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    a = Vec3(lx / 10, 0, 0) * unit.nanometer
    b = Vec3(xy / 10, ly / 10, 0) * unit.nanometer
    c = Vec3(xz / 10, yz / 10, lz / 10) * unit.nanometer

    # --- Step 1: Create OpenMM system and get PE via PythonForce ---
    print("\n--- Step 1: OpenMM PythonForce PE ---")
    torch.cuda.init()
    potential = MLPotential("uma-s-1p1-pythonforce-batch")
    system = potential.createSystem(
        topology, task_name="omol", charge=0, spin=1,
        device="cuda", use_atom_wrap_for_lammps_parity=True,
        removeCMMotion=False,
    )

    print(f"  System forces: {system.getNumForces()}")
    for i in range(system.getNumForces()):
        f = system.getForce(i)
        print(f"    Force {i}: {type(f).__name__}")

    dt = 0.1 * unit.femtoseconds
    integrator = CustomIntegrator(dt)
    integrator.addComputePerDof("v", "v + 0.5*dt*f/m")
    integrator.addComputePerDof("x", "x + dt*v")
    integrator.addUpdateContextState()
    integrator.addComputePerDof("v", "v + 0.5*dt*f/m")

    platform = Platform.getPlatformByName("CUDA")
    props = {"Precision": "mixed", "DeviceIndex": "0", "DisablePmeStream": "true"}
    context = Context(system, integrator, platform, props)
    context.setPeriodicBoxVectors(a, b, c)
    context.setPositions(pos_list)
    context.setVelocities(vel_list)

    st = context.getState(getEnergy=True, getPositions=True)
    openmm_pe_kj = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    openmm_pe_ev = openmm_pe_kj / 96.4853

    openmm_pos_nm = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    openmm_pos_ang = openmm_pos_nm * 10.0

    print(f"  OpenMM PE: {openmm_pe_kj:.2f} kJ/mol = {openmm_pe_ev:.4f} eV")
    print(f"  LAMMPS PE: {LAMMPS_PE_EV:.4f} eV")
    print(f"  Difference: {openmm_pe_ev - LAMMPS_PE_EV:.4f} eV")

    pos_diff = np.max(np.abs(openmm_pos_ang - pos_ang))
    print(f"  Max position difference (OpenMM vs LAMMPS dump): {pos_diff:.10f} A")

    # --- Step 2: Standalone direct AtomicData + predict ---
    print("\n--- Step 2: Standalone predict (OpenMM-style AtomicData) ---")
    from fairchem.core.calculate import pretrained_mlip
    from fairchem.core.datasets.atomic_data import AtomicData
    from ase.geometry import wrap_positions
    from ase.data import atomic_numbers as ase_atomic_numbers

    predict_unit = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")

    print(f"  predict_unit.direct_forces = {predict_unit.direct_forces}")

    symbols = ["O", "H", "H"] * n_mol
    atomic_nums = np.array([ase_atomic_numbers[s] for s in symbols], dtype=np.int64)
    cell_ang = np.array([[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]], dtype=np.float64)

    pos_wrapped = wrap_positions(pos_ang.copy(), cell=cell_ang, pbc=True, eps=0)
    wrap_diff = np.max(np.abs(pos_wrapped - pos_ang))
    print(f"  Max wrap displacement: {wrap_diff:.10f} A")

    data_direct = AtomicData(
        pos=torch.tensor(pos_wrapped, dtype=torch.float32),
        atomic_numbers=torch.tensor(atomic_nums),
        cell=torch.tensor(cell_ang, dtype=torch.float32).unsqueeze(0),
        pbc=torch.tensor([True, True, True], dtype=torch.bool).unsqueeze(0),
        natoms=torch.tensor([n], dtype=torch.long),
        edge_index=torch.empty((2, 0), dtype=torch.long),
        cell_offsets=torch.empty((0, 3), dtype=torch.float32),
        nedges=torch.tensor([0], dtype=torch.long),
        charge=torch.LongTensor([0]),
        spin=torch.LongTensor([1]),
        fixed=torch.zeros(n, dtype=torch.long),
        tags=torch.zeros(n, dtype=torch.long),
        batch=torch.zeros(n, dtype=torch.long),
        dataset=["omol"],
    )

    # Test A: with torch.no_grad() (as OpenMM callback does)
    with torch.no_grad():
        pred_nograd = predict_unit.predict(data_direct)
        e_nograd = float(pred_nograd["energy"].cpu().item())
        f_nograd = pred_nograd["forces"].cpu().numpy().copy()

    print(f"  Energy (with torch.no_grad): {e_nograd:.4f} eV")
    print(f"  Diff from LAMMPS: {e_nograd - LAMMPS_PE_EV:.4f} eV")

    # Test B: WITHOUT torch.no_grad() (let predict handle its own context)
    data_direct2 = AtomicData(
        pos=torch.tensor(pos_wrapped, dtype=torch.float32),
        atomic_numbers=torch.tensor(atomic_nums),
        cell=torch.tensor(cell_ang, dtype=torch.float32).unsqueeze(0),
        pbc=torch.tensor([True, True, True], dtype=torch.bool).unsqueeze(0),
        natoms=torch.tensor([n], dtype=torch.long),
        edge_index=torch.empty((2, 0), dtype=torch.long),
        cell_offsets=torch.empty((0, 3), dtype=torch.float32),
        nedges=torch.tensor([0], dtype=torch.long),
        charge=torch.LongTensor([0]),
        spin=torch.LongTensor([1]),
        fixed=torch.zeros(n, dtype=torch.long),
        tags=torch.zeros(n, dtype=torch.long),
        batch=torch.zeros(n, dtype=torch.long),
        dataset=["omol"],
    )

    pred_grad = predict_unit.predict(data_direct2)
    e_grad = float(pred_grad["energy"].cpu().item())
    f_grad = pred_grad["forces"].cpu().numpy().copy()

    print(f"  Energy (no torch.no_grad): {e_grad:.4f} eV")
    print(f"  Diff from LAMMPS: {e_grad - LAMMPS_PE_EV:.4f} eV")
    print(f"  Energy diff (no_grad vs grad): {e_nograd - e_grad:.8f} eV")

    f_diff = np.abs(f_nograd - f_grad)
    print(f"  Force diff (no_grad vs grad): max={f_diff.max():.8f} mean={f_diff.mean():.8f} eV/A")

    # --- Step 3: LAMMPS-style AtomicData ---
    print("\n--- Step 3: Standalone predict (LAMMPS-style AtomicData) ---")
    from fairchem.lammps.lammps_fc import atomic_data_from_lammps_data

    cell_t = torch.tensor(cell_ang, dtype=torch.float32).unsqueeze(0)
    data_lammps = atomic_data_from_lammps_data(
        pos_wrapped, atomic_nums, n, cell_t,
        [True, True, True], "omol",
        charge=0, spin=1,
    )

    pred_lammps = predict_unit.predict(data_lammps)
    e_lammps = float(pred_lammps["energy"].cpu().item())
    f_lammps = pred_lammps["forces"].cpu().numpy().copy()

    print(f"  Energy (LAMMPS-style): {e_lammps:.4f} eV")
    print(f"  Diff from LAMMPS: {e_lammps - LAMMPS_PE_EV:.4f} eV")
    print(f"  Diff from direct (no_grad): {e_nograd - e_lammps:.8f} eV")

    f_dl = np.abs(f_nograd - f_lammps)
    print(f"  Force diff (direct vs LAMMPS-style): max={f_dl.max():.8f} mean={f_dl.mean():.8f} eV/A")

    # --- Summary ---
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"  LAMMPS PE:             {LAMMPS_PE_EV:>16.4f} eV")
    print(f"  OpenMM PE (callback):  {openmm_pe_ev:>16.4f} eV  (diff: {openmm_pe_ev - LAMMPS_PE_EV:+.4f})")
    print(f"  Standalone (no_grad):  {e_nograd:>16.4f} eV  (diff: {e_nograd - LAMMPS_PE_EV:+.4f})")
    print(f"  Standalone (grad):     {e_grad:>16.4f} eV  (diff: {e_grad - LAMMPS_PE_EV:+.4f})")
    print(f"  Standalone (LAMMPS):   {e_lammps:>16.4f} eV  (diff: {e_lammps - LAMMPS_PE_EV:+.4f})")
    print(f"  direct_forces:         {predict_unit.direct_forces}")
    print(f"  System forces count:   {system.getNumForces()}")
    print(f"  Max pos diff:          {pos_diff:.1e} A")
    print("=" * 72)

    if abs(e_nograd - e_grad) > 0.01:
        print("\n*** CRITICAL: torch.no_grad() CHANGES the energy! ***")
        print("    The outer torch.no_grad() in compute_uma_forces_single breaks")
        print("    gradient-based force computation. Remove it.")
    elif abs(openmm_pe_ev - e_nograd) > 0.01:
        print("\n*** CRITICAL: OpenMM PythonForce PE differs from standalone! ***")
        print("    The PythonForce callback path produces different energy.")
    elif abs(e_nograd - LAMMPS_PE_EV) > 0.01:
        print("\n*** NOTE: Standalone PE differs from LAMMPS PE ***")
        print("    Could be reference energy offset (constant, doesn't affect forces).")
    else:
        print("\n  All energies match. Issue is likely in the integrator.")


if __name__ == "__main__":
    main()
