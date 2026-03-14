#!/usr/bin/env python3
"""
OpenMM classical MD matching LAMMPS in.ice_uma_quick.lmp conditions:

  Same initial structure: lammps/data.ice_uma (read_data)
  minimize (energy relaxation)
  T = 243 K, timestep 1 fs, Langevin friction 1 ps^-1 (same as fix langevin 243 243 1.0)
  Then N steps of MD (default 100 = 100 fs like LAMMPS quick deck)

Uses standard LangevinMiddleIntegrator (classical NVT), not RPMD, so the dynamical
ensemble matches LAMMPS NVE + Langevin as closely as OpenMM allows.

Prerequisites:
  pip install openmm openmmml fairchem-core
  Build data:  python lammps/build_lammps_ice_data.py -n 128

Usage:
  cd tests/uma_ice_rpmd
  python run_openmm_ice_lammps_match.py --steps 100 --platform cuda --ml-device cuda
  # Full trajectory (DCD every step + optional XYZ for Ovito):
  python run_openmm_ice_lammps_match.py --steps 100 --trajectory run.dcd --xyz run.xyz
  # Track ice disorder (Q6 ↓, q_tet ↓ as structure melts):
  python run_openmm_ice_lammps_match.py --steps 5000 --order-every 100 --order-csv order.csv
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

from openmm import (
    Context,
    LangevinMiddleIntegrator,
    LocalEnergyMinimizer,
    Platform,
    Vec3,
    unit,
)
from openmm.app import Topology, element
from openmmml import MLPotential

# Reuse parser from benchmark (same file layout)
sys.path.insert(0, str(_SCRIPT_DIR))
from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data

try:
    from ice_order_parameters import ice_order_metrics, positions_nm_to_oxygen_angstrom
except ImportError:
    ice_order_metrics = None  # type: ignore
    positions_nm_to_oxygen_angstrom = None  # type: ignore


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
        description="OpenMM UMA ice MD = LAMMPS in.ice_uma_quick conditions (classical Langevin)"
    )
    ap.add_argument("--data", type=Path, default=_LAMMPS_DIR / "data.ice_uma")
    ap.add_argument("--model", default="uma-s-1p1-pythonforce-batch")
    ap.add_argument("--temperature", type=float, default=243.0)
    ap.add_argument("--friction", type=float, default=1.0, help="1/ps, LAMMPS damp 1.0")
    ap.add_argument("--dt-fs", type=float, default=1.0)
    ap.add_argument("--steps", type=int, default=100, help="MD steps (100 @ 1fs = 100 fs)")
    ap.add_argument(
        "--minimize-iters",
        type=int,
        default=150,
        help="Max minimizer steps (each step = full UMA eval ~0.3–1 s for 128 H2O). "
        "2000+ can look 'stuck' for tens of minutes. LAMMPS allow 10k — match with --minimize-iters 10000 if patient.",
    )
    ap.add_argument("--minimize-tol", type=float, default=10.0, help="OpenMM kJ/mol/nm (looser = fewer steps)")
    ap.add_argument(
        "--skip-minimize",
        action="store_true",
        help="Skip minimization (fast path; use if geometry already relaxed or debugging MD only)",
    )
    ap.add_argument("--platform", default="cuda", choices=["cuda", "cpu"])
    ap.add_argument("--ml-device", default=None)
    ap.add_argument("--seed", type=int, default=284759, help="LAMMPS velocity seed analogue")
    ap.add_argument("-o", "--output", type=Path, default=None, help="Final PDB")
    ap.add_argument(
        "--trajectory",
        type=Path,
        default=None,
        metavar="PATH.dcd",
        help="Write full trajectory (DCD). Default: <output_stem>.dcd if -o set, else ice_openmm_lammps_match.dcd",
    )
    ap.add_argument(
        "--traj-every",
        type=int,
        default=1,
        help="Write DCD/XYZ every N steps (1 = every step, like LAMMPS dump 1)",
    )
    ap.add_argument(
        "--xyz",
        type=Path,
        default=None,
        metavar="PATH.xyz",
        help="Optional extended trajectory as multi-frame XYZ (Ovito-friendly; O/H labels)",
    )
    ap.add_argument(
        "--no-trajectory",
        action="store_true",
        help="No DCD (faster MD; only final PDB). Default is to write <output>.dcd every step.",
    )
    ap.add_argument("--report-every", type=int, default=10)
    ap.add_argument(
        "--order-every",
        type=int,
        default=0,
        help="If >0, print mean Q6 + q_tet (ice order) every N steps. Requires scipy. "
        "Q6 = Steinhardt l=6 (drops as ice disorders); q_tet = tetrahedral order 0–1.",
    )
    ap.add_argument(
        "--order-csv",
        type=Path,
        default=None,
        metavar="PATH.csv",
        help="Write order metrics CSV (also sets --order-every to --report-every if you omit --order-every)",
    )
    args = ap.parse_args()

    if args.order_csv is not None and args.order_every <= 0:
        args.order_every = max(1, args.report_every)
        print(
            f"  NOTE: --order-csv without --order-every → using --order-every={args.order_every} "
            f"(same as report interval). CSV updates every {args.order_every} steps.",
            flush=True,
        )

    if not args.data.is_file():
        print(f"Missing {args.data} — run: python {_LAMMPS_DIR / 'build_lammps_ice_data.py'} -n 128")
        sys.exit(1)

    pos_ang, cell, Z = parse_lammps_data(args.data)
    n = Z.shape[0]
    assert n % 3 == 0
    n_mol = n // 3
    # O,H,H order must match data file (type 1 O, 2 H — build script order)
    for i in range(n_mol):
        assert Z[3 * i] == 8 and Z[3 * i + 1] == 1 and Z[3 * i + 2] == 1

    topology = water_topology(n_mol)
    positions_nm = pos_ang / 10.0
    pos_list = [
        Vec3(float(positions_nm[i, 0]), float(positions_nm[i, 1]), float(positions_nm[i, 2])) * unit.nanometer
        for i in range(n)
    ]
    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    a = Vec3(lx / 10, 0, 0) * unit.nanometer
    b = Vec3(xy / 10, ly / 10, 0) * unit.nanometer
    c = Vec3(xz / 10, yz / 10, lz / 10) * unit.nanometer
    box = (a, b, c)

    print("=" * 72)
    print("OpenMM MD — LAMMPS-matched (data.ice_uma + Langevin 243 K, 1 fs)")
    print("=" * 72)
    print(f"  Atoms: {n}  molecules: {n_mol}")
    print(f"  Box (nm): a={a[0].value_in_unit(unit.nanometer):.4f} …")
    print(f"  Model: {args.model}  platform: {args.platform}")

    ml_dev = args.ml_device or ("cuda" if args.platform == "cuda" else None)
    if args.platform == "cuda" and ml_dev == "cuda":
        import torch
        if torch.cuda.is_available():
            torch.cuda.init()

    potential = MLPotential(args.model)
    system = potential.createSystem(
        topology, task_name="omol", charge=0, spin=1, device=ml_dev
    )

    integrator = LangevinMiddleIntegrator(
        args.temperature * unit.kelvin,
        args.friction / unit.picosecond,
        args.dt_fs * unit.femtoseconds,
    )
    integrator.setRandomNumberSeed(args.seed)

    platform = Platform.getPlatformByName(args.platform.upper())
    props = {}
    if args.platform == "cuda":
        props = {"Precision": "mixed", "DeviceIndex": "0", "DisablePmeStream": "true"}

    context = Context(system, integrator, platform, props)
    context.setPeriodicBoxVectors(a, b, c)
    context.setPositions(pos_list)

    print("\n--- Minimize (cf. LAMMPS minimize) ---")
    if args.skip_minimize:
        print("  --skip-minimize: no minimization (instant).")
    else:
        # Each iteration calls UMA once — no console output until done → looks "stuck".
        est_s = args.minimize_iters * 0.4 * (n / 384.0)
        print(
            f"  WARNING: Minimizer runs up to {args.minimize_iters} UMA force evals (~{est_s:.0f} s wall time "
            f"if all used). No progress lines until finished. Use --skip-minimize or lower --minimize-iters to test MD.",
            flush=True,
        )
        t0 = time.perf_counter()
        LocalEnergyMinimizer.minimize(
            context, maxIterations=args.minimize_iters, tolerance=args.minimize_tol
        )
        st = context.getState(getEnergy=True)
        print(f"  PE after minimize: {st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.2f} kJ/mol")
        print(f"  Minimize wall time: {time.perf_counter() - t0:.1f} s", flush=True)

    context.setVelocitiesToTemperature(args.temperature * unit.kelvin, args.seed)

    topology.setPeriodicBoxVectors((a, b, c))
    out = args.output or (
        _SCRIPT_DIR / f"ice_openmm_lammps_match_T{args.temperature:.0f}_steps{args.steps}.pdb"
    )
    if args.no_trajectory:
        dcd_path = None
    elif args.trajectory is not None:
        dcd_path = args.trajectory
    else:
        dcd_path = out.with_suffix(".dcd")

    dcd = None
    xyz_f = None
    symbols = ["O", "H", "H"] * n_mol
    if dcd_path:
        from openmm.app import DCDFile

        dcd_f = open(dcd_path, "wb")
        dcd = DCDFile(
            dcd_f,
            topology,
            args.dt_fs * unit.femtoseconds,
            interval=args.traj_every,
        )
        print(f"  Trajectory DCD: {dcd_path} (every {args.traj_every} step(s))")
    if args.xyz:
        xyz_f = open(args.xyz, "w")
        print(f"  Trajectory XYZ: {args.xyz}")

    # Orthorhombic cell (Å) → nm; DCD triclinic pack can fail if angles aren’t plain floats
    _uc_nm = Vec3(float(lx) * 0.1, float(ly) * 0.1, float(lz) * 0.1) * unit.nanometer

    def _write_traj_frame(step_index: int, pos_nm) -> None:
        if dcd is not None and (step_index == 0 or step_index % args.traj_every == 0):
            dcd.writeModel(pos_nm, unitCellDimensions=_uc_nm)
        if xyz_f is not None and (step_index == 0 or step_index % args.traj_every == 0):
            pos_ang = np.array(
                [[p[i].value_in_unit(unit.nanometer) * 10.0 for i in range(3)] for p in pos_nm]
            )
            xyz_f.write(f"{n}\n")
            xyz_f.write(f"Step {step_index} OpenMM UMA\n")
            for i in range(n):
                xyz_f.write(f"{symbols[i]:2s} {pos_ang[i,0]:.8f} {pos_ang[i,1]:.8f} {pos_ang[i,2]:.8f}\n")
            xyz_f.flush()

    st0 = context.getState(getPositions=True)
    _write_traj_frame(0, st0.getPositions())

    box_ang = np.array([float(lx), float(ly), float(lz)], dtype=np.float64)
    order_csv_f = None
    if args.order_every > 0:
        if ice_order_metrics is None:
            print("  WARNING: --order-every ignored (ice_order_parameters import failed)", flush=True)
        else:
            print(
                f"  Order params every {args.order_every} steps: Q6 (Steinhardt), q_tet (1=ideal tetrahedron)",
                flush=True,
            )
            if args.order_csv:
                order_csv_f = open(args.order_csv, "w")
                order_csv_f.write("step,time_ps,T_K,PE_kj_mol,q6_mean,q6_std,q_tet_mean,q_tet_std\n")

    print(f"\n--- MD {args.steps} steps @ {args.dt_fs} fs (cf. LAMMPS run 100) ---")
    dt_ps = args.dt_fs / 1000.0
    t_md = time.perf_counter()
    # Same T as OpenMM StateDataReporter: T = 2*KE/(dof*R), KE = getKineticEnergy() [kJ/mol]
    dof = 3 * n
    R = unit.MOLAR_GAS_CONSTANT_R

    def _instantaneous_T_K(st_) -> float:
        if hasattr(integrator, "computeSystemTemperature"):
            return integrator.computeSystemTemperature().value_in_unit(unit.kelvin)
        ke = st_.getKineticEnergy()
        return (2 * ke / (dof * R)).value_in_unit(unit.kelvin)

    write_traj = dcd is not None or xyz_f is not None
    for step in range(1, args.steps + 1):
        integrator.step(1)
        if write_traj and step % args.traj_every == 0:
            stp = context.getState(getPositions=True)
            _write_traj_frame(step, stp.getPositions())
        if step % args.report_every == 0 or step == args.steps:
            st = context.getState(getEnergy=True, getVelocities=True)
            pe = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            temp = _instantaneous_T_K(st)
            line = f"  step {step:5d}  t={step * dt_ps:.4f} ps  T={temp:7.2f} K  PE={pe:14.2f} kJ/mol"
            if (
                args.order_every > 0
                and ice_order_metrics is not None
                and (step % args.order_every == 0 or step == args.steps)
            ):
                pos_nm = context.getState(getPositions=True).getPositions()
                pos_ang = positions_nm_to_oxygen_angstrom(pos_nm)
                om = ice_order_metrics(pos_ang, box_ang)
                line += f"  Q6={om.q6_mean:.3f}  q_tet={om.q_tet_mean:.3f}"
                if order_csv_f is not None:
                    order_csv_f.write(
                        f"{step},{step * dt_ps:.6f},{temp:.4f},{pe:.6f},"
                        f"{om.q6_mean:.6f},{om.q6_std:.6f},{om.q_tet_mean:.6f},{om.q_tet_std:.6f}\n"
                    )
                    order_csv_f.flush()
            print(line)

    elapsed = time.perf_counter() - t_md
    ps_done = args.steps * dt_ps
    ns = ps_done * 1e-3
    ns_per_day = ns / elapsed * 86400.0 if elapsed > 0 else 0.0
    steps_per_s = args.steps / elapsed if elapsed > 0 else 0.0
    print(f"  MD wall time: {elapsed:.1f} s")
    print(f"  Simulation speed: {ns_per_day:.3f} ns/day  |  {steps_per_s:.2f} steps/s  |  {ps_done:.4f} ps integrated")

    if dcd is not None:
        dcd_f.close()
        print(f"  Trajectory (DCD): {dcd_path}  ({1 + args.steps // args.traj_every} frames @ every {args.traj_every} step)")
    if xyz_f is not None:
        xyz_f.close()
        print(f"  Trajectory (XYZ): {args.xyz}")
    if order_csv_f is not None:
        order_csv_f.close()
        print(f"  Order CSV: {args.order_csv}")

    from openmm.app import PDBFile

    # PDB CRYST1 uses % float formatting; orthogonal cell avoids Quantity angles from triclinic.
    topology.setUnitCellDimensions(
        Vec3(float(lx) * 0.1, float(ly) * 0.1, float(lz) * 0.1) * unit.nanometer
    )
    PDBFile.writeFile(
        topology,
        context.getState(getPositions=True).getPositions(),
        open(out, "w"),
        keepIds=True,
    )
    print(f"  Final PDB: {out}")
    print("=" * 72)


if __name__ == "__main__":
    main()
