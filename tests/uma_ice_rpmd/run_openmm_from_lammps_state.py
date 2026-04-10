#!/usr/bin/env python3
"""
Run OpenMM UMA from an exact LAMMPS state (positions + velocities) in pure NVE.

Reads a LAMMPS custom dump containing id type mass x y z vx vy vz,
converts units (A -> nm, A/ps -> nm/ps), creates a UMA system with a
Velocity-Verlet CustomIntegrator (no thermostat), and runs NVE dynamics
while logging Q6 and q_tet order parameters.

This is the definitive test: if Q6 tracks match between LAMMPS NVE and
OpenMM NVE from the same state, the codes are equivalent.

Usage:
  python run_openmm_from_lammps_state.py \
      --state-dump lammps/dump.nve_state_5000.lammpstrj \
      --data lammps/data.ice_uma \
      --steps 3000 --dt-fs 0.1 --order-every 100 \
      --order-csv pipeline_out/ice_order_openmm_nve.csv
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np

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

try:
    from ice_order_parameters import ice_order_metrics, positions_nm_to_oxygen_angstrom
except ImportError:
    ice_order_metrics = None  # type: ignore
    positions_nm_to_oxygen_angstrom = None  # type: ignore


def read_lammps_state_dump(dump_path: Path):
    """Parse a single-frame LAMMPS custom dump with id type mass x y z vx vy vz.

    Returns
    -------
    pos_ang : (n, 3) float64 -- positions in Angstrom, sorted by atom id
    vel_ang_ps : (n, 3) float64 -- velocities in A/ps, sorted by atom id
    types : (n,) int -- LAMMPS atom types, sorted by atom id
    """
    lines = dump_path.read_text().splitlines()
    i = 0
    n_lines = len(lines)
    while i < n_lines:
        if lines[i].startswith("ITEM: TIMESTEP"):
            step = int(lines[i + 1])
            i += 2
            break
        i += 1

    while i < n_lines and not lines[i].startswith("ITEM: NUMBER OF ATOMS"):
        i += 1
    nat = int(lines[i + 1])
    i += 2

    while i < n_lines and not lines[i].startswith("ITEM: ATOMS"):
        i += 1
    header = lines[i].split()[2:]
    idx_id = header.index("id")
    idx_type = header.index("type")
    idx_x = header.index("x")
    idx_vx = header.index("vx")
    i += 1

    rows = []
    for _ in range(nat):
        p = lines[i].split()
        aid = int(p[idx_id])
        typ = int(p[idx_type])
        x, y, z = float(p[idx_x]), float(p[idx_x + 1]), float(p[idx_x + 2])
        vx, vy, vz = float(p[idx_vx]), float(p[idx_vx + 1]), float(p[idx_vx + 2])
        rows.append((aid, typ, x, y, z, vx, vy, vz))
        i += 1

    rows.sort(key=lambda r: r[0])
    pos_ang = np.array([[r[2], r[3], r[4]] for r in rows], dtype=np.float64)
    vel_ang_ps = np.array([[r[5], r[6], r[7]] for r in rows], dtype=np.float64)
    types = np.array([r[1] for r in rows], dtype=np.int64)
    print(f"  Read state dump: step={step}, {nat} atoms")
    print(f"  Pos range: {pos_ang.min():.4f} .. {pos_ang.max():.4f} A")
    print(f"  Vel range: {vel_ang_ps.min():.4f} .. {vel_ang_ps.max():.4f} A/ps")
    return pos_ang, vel_ang_ps, types


def _make_velocity_verlet(dt_fs: float) -> CustomIntegrator:
    """Pure Velocity-Verlet NVE integrator (no thermostat)."""
    dt = dt_fs * unit.femtoseconds
    integrator = CustomIntegrator(dt)
    integrator.addComputePerDof("v", "v + 0.5*dt*f/m")
    integrator.addComputePerDof("x", "x + dt*v")
    integrator.addUpdateContextState()
    integrator.addComputePerDof("v", "v + 0.5*dt*f/m")
    return integrator


def water_topology(n_molecules: int) -> Topology:
    top = Topology()
    chain = top.addChain()
    for _ in range(n_molecules):
        r = top.addResidue("HOH", chain)
        top.addAtom("O", element.oxygen, r)
        top.addAtom("H1", element.hydrogen, r)
        top.addAtom("H2", element.hydrogen, r)
    return top


def main() -> None:
    ap = argparse.ArgumentParser(
        description="OpenMM UMA NVE from LAMMPS state — velocity-matched deterministic comparison"
    )
    ap.add_argument(
        "--state-dump", type=Path, required=True,
        help="LAMMPS dump with id type mass x y z vx vy vz (single frame at thermalization)",
    )
    ap.add_argument("--data", type=Path, default=_LAMMPS_DIR / "data.ice_uma",
                     help="LAMMPS data file (for box dimensions)")
    ap.add_argument("--model", default="uma-s-1p1-pythonforce-batch")
    ap.add_argument("--dt-fs", type=float, default=0.1, help="Timestep in fs (match LAMMPS)")
    ap.add_argument("--steps", type=int, default=3000, help="NVE steps")
    ap.add_argument("--platform", default="cuda", choices=["cuda", "cpu"])
    ap.add_argument("--ml-device", default=None)
    ap.add_argument("--report-every", type=int, default=100)
    ap.add_argument("--order-every", type=int, default=100, help="Q6/q_tet every N steps")
    ap.add_argument("--order-csv", type=Path, default=None, help="Output CSV for order params")
    args = ap.parse_args()

    if not args.state_dump.is_file():
        print(f"Missing state dump: {args.state_dump}")
        sys.exit(1)
    if not args.data.is_file():
        print(f"Missing data file: {args.data}")
        sys.exit(1)

    _, cell, Z_data = parse_lammps_data(args.data)
    pos_ang, vel_ang_ps, types = read_lammps_state_dump(args.state_dump)

    n = pos_ang.shape[0]
    assert n == Z_data.shape[0], f"Atom count mismatch: dump {n} vs data {Z_data.shape[0]}"
    n_mol = n // 3
    assert n % 3 == 0

    Z = np.array([8 if t == 1 else 1 for t in types], dtype=np.int64)
    for i in range(n_mol):
        assert Z[3 * i] == 8 and Z[3 * i + 1] == 1 and Z[3 * i + 2] == 1, \
            f"Unexpected atom order at molecule {i}: Z={Z[3*i:3*i+3]}"

    topology = water_topology(n_mol)

    # Unit conversion: A -> nm
    positions_nm = pos_ang / 10.0
    pos_list = [
        Vec3(float(positions_nm[i, 0]), float(positions_nm[i, 1]), float(positions_nm[i, 2])) * unit.nanometer
        for i in range(n)
    ]

    # Unit conversion: A/ps -> nm/ps
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

    print("=" * 72)
    print("OpenMM NVE from LAMMPS state — Velocity-Verlet (no thermostat)")
    print("=" * 72)
    print(f"  Atoms: {n}  molecules: {n_mol}")
    print(f"  Box (nm): a={lx/10:.4f}  b={ly/10:.4f}  c={lz/10:.4f}")
    print(f"  Model: {args.model}  platform: {args.platform}")
    print(f"  dt: {args.dt_fs} fs  steps: {args.steps}  total: {args.steps * args.dt_fs / 1000:.4f} ps")

    ml_dev = args.ml_device or ("cuda" if args.platform == "cuda" else None)
    if args.platform == "cuda" and ml_dev == "cuda":
        import torch
        if torch.cuda.is_available():
            torch.cuda.init()

    potential = MLPotential(args.model)
    system = potential.createSystem(
        topology,
        task_name="omol",
        charge=0,
        spin=1,
        device=ml_dev,
        use_atom_wrap_for_lammps_parity=True,
        removeCMMotion=False,
    )

    integrator = _make_velocity_verlet(args.dt_fs)
    print(f"  Integrator: Velocity-Verlet NVE (no thermostat)")
    print(f"  CMMotionRemover: DISABLED (matching LAMMPS NVE)")

    platform = Platform.getPlatformByName(args.platform.upper())
    props = {}
    if args.platform == "cuda":
        props = {"Precision": "mixed", "DeviceIndex": "0", "DisablePmeStream": "true"}

    context = Context(system, integrator, platform, props)
    context.setPeriodicBoxVectors(a, b, c)
    context.setPositions(pos_list)
    context.setVelocities(vel_list)

    st0 = context.getState(getEnergy=True, getPositions=True, getVelocities=True)
    pe0 = st0.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    ke0 = st0.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
    v0 = st0.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
    v_max = np.max(np.abs(v0))
    print(f"\n  Initial PE: {pe0:.2f} kJ/mol  KE: {ke0:.2f} kJ/mol  E_tot: {pe0+ke0:.2f} kJ/mol")
    print(f"  Max |v|: {v_max:.6f} nm/ps")

    box_ang = np.array([float(lx), float(ly), float(lz)], dtype=np.float64)
    dof = 3 * n
    R = unit.MOLAR_GAS_CONSTANT_R

    def _instantaneous_T_K(st_) -> float:
        ke = st_.getKineticEnergy()
        return (2 * ke / (dof * R)).value_in_unit(unit.kelvin)

    order_csv_f = None
    if args.order_csv:
        order_csv_f = open(args.order_csv, "w")
        order_csv_f.write("step,time_ps,T_K,PE_kj_mol,KE_kj_mol,E_tot_kj_mol,q6_mean,q6_std,q_tet_mean,q_tet_std\n")

    if ice_order_metrics is not None:
        pos_nm_init = st0.getPositions()
        pos_ang_init = positions_nm_to_oxygen_angstrom(pos_nm_init)
        om0 = ice_order_metrics(pos_ang_init, box_ang)
        print(f"  Initial Q6={om0.q6_mean:.4f}  q_tet={om0.q_tet_mean:.4f}")
        if order_csv_f is not None:
            order_csv_f.write(
                f"0,0.000000,{_instantaneous_T_K(st0):.4f},{pe0:.6f},{ke0:.6f},{pe0+ke0:.6f},"
                f"{om0.q6_mean:.6f},{om0.q6_std:.6f},{om0.q_tet_mean:.6f},{om0.q_tet_std:.6f}\n"
            )
            order_csv_f.flush()

    print(f"\n--- NVE MD {args.steps} steps @ {args.dt_fs} fs ---")
    dt_ps = args.dt_fs / 1000.0
    t_md = time.perf_counter()

    for step in range(1, args.steps + 1):
        integrator.step(1)
        if step % args.report_every == 0 or step == args.steps:
            st = context.getState(getEnergy=True)
            pe = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            ke = st.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
            temp = _instantaneous_T_K(st)
            e_tot = pe + ke
            line = (
                f"  step {step:5d}  t={step * dt_ps:.4f} ps  T={temp:7.2f} K  "
                f"PE={pe:14.2f}  KE={ke:14.2f}  E_tot={e_tot:14.2f} kJ/mol"
            )
            if (
                ice_order_metrics is not None
                and args.order_every > 0
                and (step % args.order_every == 0 or step == args.steps)
            ):
                pos_nm = context.getState(getPositions=True).getPositions()
                pos_a = positions_nm_to_oxygen_angstrom(pos_nm)
                om = ice_order_metrics(pos_a, box_ang)
                line += f"  Q6={om.q6_mean:.4f}  q_tet={om.q_tet_mean:.4f}"
                if order_csv_f is not None:
                    order_csv_f.write(
                        f"{step},{step * dt_ps:.6f},{temp:.4f},{pe:.6f},{ke:.6f},{e_tot:.6f},"
                        f"{om.q6_mean:.6f},{om.q6_std:.6f},{om.q_tet_mean:.6f},{om.q_tet_std:.6f}\n"
                    )
                    order_csv_f.flush()
            print(line, flush=True)

    elapsed = time.perf_counter() - t_md
    ps_done = args.steps * dt_ps
    steps_per_s = args.steps / elapsed if elapsed > 0 else 0.0
    print(f"\n  MD wall time: {elapsed:.1f} s")
    print(f"  {steps_per_s:.2f} steps/s  |  {ps_done:.4f} ps integrated")

    st_final = context.getState(getEnergy=True)
    pe_f = st_final.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    ke_f = st_final.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Final PE: {pe_f:.2f}  KE: {ke_f:.2f}  E_tot: {pe_f+ke_f:.2f} kJ/mol")
    print(f"  Energy drift: {(pe_f+ke_f) - (pe0+ke0):.6f} kJ/mol")

    if order_csv_f is not None:
        order_csv_f.close()
        print(f"  Order CSV: {args.order_csv}")
    print("=" * 72)


if __name__ == "__main__":
    main()
