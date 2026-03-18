#!/usr/bin/env python3
"""
OpenMM RPMD with TIP4P/2005f flexible water on ice.

Runs path-integral RPMD (PILE thermostat) using TIP4P/2005f forces instead of UMA.
Same structure source as UMA/LAMMPS runs (LAMMPS data file). Order parameters
use centroid bead positions (Q6, q_tet) for comparison with OpenMM UMA and
i-PI+LAMMPS UMA RPMD.

Prerequisites:
  pip install openmm scipy
  Build data:  python lammps/build_lammps_ice_data.py -n 128

Usage:
  cd tests/uma_ice_rpmd
  python run_openmm_tip4p_rpmd.py --data lammps/data.ice_uma_128 --beads 8 --dt 0.1 --prod 100
  python run_openmm_tip4p_rpmd.py --molecules 32 --beads 4 --order-csv pipeline_out/ice_order_tip4p_rpmd.csv --seed 284760
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
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
    LocalEnergyMinimizer,
    Platform,
    RPMDIntegrator,
    Vec3,
    unit,
)
from openmm.app import ForceField, Modeller, PDBFile, PME

sys.path.insert(0, str(_SCRIPT_DIR))
from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data

try:
    from ice_order_parameters import ice_order_metrics
except ImportError:
    ice_order_metrics = None

_TIP4P2005F_R_OH_NM = 0.09419
_TIP4P2005F_ANGLE_RAD = 1.87448  # 107.4°


def _tip4p2005f_forcefield_path() -> Path:
    import openmm.app as _app
    pkg = Path(_app.__file__).parent / "data" / "tip4p2005f.xml"
    if pkg.is_file():
        return pkg
    repo = _SCRIPT_DIR.parent.parent / "wrappers/python/openmm/app/data/tip4p2005f.xml"
    if repo.is_file():
        return repo
    raise FileNotFoundError("tip4p2005f.xml not found. Pass --force-field path.")


def _water_pdb_from_lammps(pos_ang: np.ndarray, cell: np.ndarray, n_mol: int) -> str:
    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    lines = [f"CRYST1{lx:9.3f}{ly:9.3f}{lz:9.3f}  90.00  90.00  90.00 P 1           1"]
    for i in range(n_mol):
        o, h1, h2 = pos_ang[3 * i], pos_ang[3 * i + 1], pos_ang[3 * i + 2]
        res = i + 1
        lines.append(f"HETATM{i*3+1:5d}  O   HOH  {res:4d}    {o[0]:8.3f}{o[1]:8.3f}{o[2]:8.3f}  1.00  0.00           O")
        lines.append(f"HETATM{i*3+2:5d}  H1  HOH  {res:4d}    {h1[0]:8.3f}{h1[1]:8.3f}{h1[2]:8.3f}  1.00  0.00           H")
        lines.append(f"HETATM{i*3+3:5d}  H2  HOH  {res:4d}    {h2[0]:8.3f}{h2[1]:8.3f}{h2[2]:8.3f}  1.00  0.00           H")
    lines.append("END")
    return "\n".join(lines)


def _relax_water_geometry_to_tip4p2005f(pos_nm: np.ndarray, n_mol: int, box_nm: np.ndarray) -> np.ndarray:
    out = np.array(pos_nm, dtype=np.float64, copy=True)
    r_eq = _TIP4P2005F_R_OH_NM
    half_angle = _TIP4P2005F_ANGLE_RAD / 2.0
    w1, w2 = 0.73612, 0.13194
    for m in range(n_mol):
        oi, h1i, h2i, mi = 4 * m, 4 * m + 1, 4 * m + 2, 4 * m + 3
        o, h1, h2 = out[oi], out[h1i], out[h2i]

        def wrap(d):
            d = np.array(d, dtype=np.float64)
            for k in range(3):
                L = box_nm[k]
                d[k] -= np.round(d[k] / L) * L
            return d

        d1 = wrap(h1 - o)
        d2 = wrap(h2 - o)
        bisector = d1 + d2
        b_norm = np.linalg.norm(bisector)
        if b_norm < 1e-10:
            continue
        b_hat = bisector / b_norm
        u = d1 - np.dot(d1, b_hat) * b_hat
        u_norm = np.linalg.norm(u)
        if u_norm < 1e-10:
            u = d2 - np.dot(d2, b_hat) * b_hat
            u_norm = np.linalg.norm(u)
        u_hat = u / u_norm if u_norm >= 1e-10 else np.array([1.0, 0.0, 0.0])
        h1_new = o + r_eq * (np.cos(half_angle) * b_hat + np.sin(half_angle) * u_hat)
        h2_new = o + r_eq * (np.cos(half_angle) * b_hat - np.sin(half_angle) * u_hat)
        out[h1i], out[h2i] = h1_new, h2_new
        out[mi] = w1 * o + w2 * h1_new + w2 * h2_new
    return out


def _mic_centroid(positions_per_bead: list[np.ndarray], box_matrix: np.ndarray) -> np.ndarray:
    ref = positions_per_bead[0]
    inv_box = np.linalg.inv(box_matrix)
    accum = np.zeros_like(ref)
    for b in range(len(positions_per_bead)):
        delta = positions_per_bead[b] - ref
        frac = delta @ inv_box.T
        frac -= np.round(frac)
        accum += frac @ box_matrix
    return ref + accum / len(positions_per_bead)


def main() -> None:
    ap = argparse.ArgumentParser(description="OpenMM TIP4P/2005f RPMD ice simulation")
    ap.add_argument("--data", type=Path, default=None)
    ap.add_argument("--molecules", type=int, default=128)
    ap.add_argument("--beads", type=int, default=8)
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--prod", type=float, default=100.0)
    ap.add_argument("--temperature", type=float, default=243.0)
    ap.add_argument("--order-csv", type=Path, default=None)
    ap.add_argument("--order-every", type=int, default=10)
    ap.add_argument("--seed", type=int, default=284759)
    ap.add_argument("--platform", default="cuda")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--force-field", type=Path, default=None)
    ap.add_argument("--skip-minimize", action="store_true")
    args = ap.parse_args()

    data_path = args.data or (_LAMMPS_DIR / f"data.ice_uma_{args.molecules}")
    if not data_path.is_file():
        print(f"Missing {data_path}. Run: python lammps/build_lammps_ice_data.py -n {args.molecules} -o {data_path}")
        sys.exit(1)

    ff_path = args.force_field or _tip4p2005f_forcefield_path()
    if not ff_path.is_file():
        print(f"Force field not found: {ff_path}")
        sys.exit(1)

    pos_ang, cell, Z = parse_lammps_data(data_path)
    n_real = Z.shape[0]
    assert n_real % 3 == 0
    n_mol = n_real // 3
    if args.molecules is not None and n_mol > args.molecules:
        n_mol = args.molecules
        n_real = n_mol * 3
        pos_ang = pos_ang[:n_real]
        Z = Z[:n_real]

    pdb_str = _water_pdb_from_lammps(pos_ang, cell, n_mol)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
        f.write(pdb_str)
        pdb_path = f.name
    try:
        pdb = PDBFile(pdb_path)
    finally:
        Path(pdb_path).unlink(missing_ok=True)

    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    a = Vec3(lx / 10, 0, 0) * unit.nanometer
    b = Vec3(xy / 10, ly / 10, 0) * unit.nanometer
    c = Vec3(xz / 10, yz / 10, lz / 10) * unit.nanometer
    box_matrix = np.array([
        [lx / 10, 0, 0],
        [xy / 10, ly / 10, 0],
        [xz / 10, yz / 10, lz / 10],
    ])

    ff = ForceField(str(ff_path))
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addExtraParticles(ff)
    n_total = modeller.topology.getNumAtoms()
    assert n_total == 4 * n_mol

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=0.7 * unit.nanometer,
        constraints=None,
        rigidWater=False,
        ewaldErrorTolerance=0.0005,
    )

    dt_fs = args.dt
    friction = 1.0
    integrator = RPMDIntegrator(
        args.beads,
        args.temperature * unit.kelvin,
        friction / unit.picosecond,
        dt_fs * unit.femtoseconds,
    )
    integrator.setThermostatType(RPMDIntegrator.Pile)
    integrator.setRandomNumberSeed(args.seed)

    platform = Platform.getPlatformByName(args.platform.upper())
    props = {"Precision": "mixed", "DeviceIndex": "0"} if args.platform == "cuda" else {}
    context = Context(system, integrator, platform, props)
    context.setPeriodicBoxVectors(a, b, c)
    context.setPositions(modeller.positions)

    pos_nm = np.array([[modeller.positions[i][j].value_in_unit(unit.nanometer) for j in range(3)] for i in range(n_total)])
    box_nm = np.array([lx, ly, lz]) / 10.0
    pos_nm = _relax_water_geometry_to_tip4p2005f(pos_nm, n_mol, box_nm)
    context.setPositions([Vec3(float(p[0]), float(p[1]), float(p[2])) * unit.nanometer for p in pos_nm])

    if not args.skip_minimize:
        LocalEnergyMinimizer.minimize(context, tolerance=1.0, maxIterations=2000)
    min_pos = context.getState(getPositions=True).getPositions()

    np.random.seed(args.seed)
    context.setPositions(min_pos)
    context.setVelocitiesToTemperature(args.temperature * unit.kelvin)
    base_vel = context.getState(getVelocities=True).getVelocities()
    for i in range(args.beads):
        vel_list = []
        for v in base_vel:
            vx = v[0].value_in_unit(unit.nanometer / unit.picosecond) + 0.01 * np.random.randn()
            vy = v[1].value_in_unit(unit.nanometer / unit.picosecond) + 0.01 * np.random.randn()
            vz = v[2].value_in_unit(unit.nanometer / unit.picosecond) + 0.01 * np.random.randn()
            vel_list.append(Vec3(vx, vy, vz) * unit.nanometer / unit.picosecond)
        integrator.setVelocities(i, vel_list)
        integrator.setPositions(i, min_pos)

    dt_ps = dt_fs / 1000.0
    total_steps = int(args.prod / dt_ps)
    order_every = args.order_every or 10
    report_every = max(1, int(1.0 / dt_ps))

    box_orth_ang = np.linalg.norm(box_matrix, axis=1) * 10.0
    z_arr = np.array([8, 1, 1, 0] * n_mol, dtype=np.int64)

    order_csv_f = None
    if args.order_csv:
        args.order_csv.parent.mkdir(parents=True, exist_ok=True)
        order_csv_f = open(args.order_csv, "w")
        order_csv_f.write("step,time_ps,T_K,PE_kj_mol,q6_mean,q6_std,q_tet_mean,q_tet_std\n")

    if ice_order_metrics and order_csv_f:
        bead0 = integrator.getState(0, getPositions=True)
        init_pos = np.array(bead0.getPositions(asNumpy=True).value_in_unit(unit.nanometer)) * 10.0
        om = ice_order_metrics(init_pos, box_orth_ang, z=z_arr)
        order_csv_f.write(f"0,0.000000,,,{om.q6_mean:.6f},{om.q6_std:.6f},{om.q_tet_mean:.6f},{om.q_tet_std:.6f}\n")

    print(f"TIP4P/2005f RPMD: {n_mol} mol, {args.beads} beads, {dt_fs} fs, {args.prod} ps, {total_steps} steps")
    start = time.time()
    for step in range(1, total_steps + 1):
        integrator.step(1)
        if step % report_every == 0:
            t_ps = step * dt_ps
            state = integrator.getState(0, getEnergy=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            ns_day = (step / (time.time() - start)) * dt_ps * 86400 / 1000.0
            print(f"  step {step:6d}  t={t_ps:.3f} ps  PE={pe:.1f}  {ns_day:.2f} ns/day")
        if order_csv_f and step % order_every == 0 and ice_order_metrics:
            bead_pos_list = []
            for b in range(args.beads):
                s = integrator.getState(b, getPositions=True)
                bead_pos_list.append(s.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
            centroid = _mic_centroid(bead_pos_list, box_matrix)
            pos_ang = centroid * 10.0
            om = ice_order_metrics(pos_ang, box_orth_ang, z=z_arr)
            t_ps = step * dt_ps
            order_csv_f.write(f"{step},{t_ps:.6f},,,{om.q6_mean:.6f},{om.q6_std:.6f},{om.q_tet_mean:.6f},{om.q_tet_std:.6f}\n")
            order_csv_f.flush()

    if order_csv_f:
        order_csv_f.close()
        print(f"Order CSV: {args.order_csv}")
    print("Done.")


if __name__ == "__main__":
    main()
