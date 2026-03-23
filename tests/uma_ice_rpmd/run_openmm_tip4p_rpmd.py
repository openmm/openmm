#!/usr/bin/env python3
"""
OpenMM RPMD with TIP4P/2005f flexible water on ice.

Runs path-integral RPMD using TIP4P/2005f forces instead of UMA. Default bath:
**PILE_G** (Bussi/SVR on centroid normal mode + PILE on internal modes), matching
OpenMM ``RPMDIntegrator.PileG`` and the UMA ice script defaults. Same structure
source as UMA/LAMMPS runs (LAMMPS data file). Order parameters use centroid bead
positions (Q6, q_tet) for comparison with OpenMM UMA RPMD.

Prerequisites:
  pip install openmm scipy
  Build data:  python lammps/build_lammps_ice_data.py --nx 2 --ny 2 --nz 4 -o lammps/data.ice_uma_128

Usage:
  cd tests/uma_ice_rpmd
  python run_openmm_tip4p_rpmd.py --data lammps/data.ice_uma_128 --beads 8 --dt 0.1 --prod 100
  python run_openmm_tip4p_rpmd.py --molecules 64 --beads 4 --order-csv pipeline_out/ice_order_tip4p_rpmd.csv --seed 284760
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
    from ice_order_parameters import (
        ICE_ORDER_PI_CSV_HEADER,
        format_ice_order_pi_csv_row,
        ice_order_metrics_path_integral,
    )
except ImportError:
    ICE_ORDER_PI_CSV_HEADER = ""
    format_ice_order_pi_csv_row = None  # type: ignore
    ice_order_metrics_path_integral = None  # type: ignore

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


def _auto_pme_cutoff_nm(cell: np.ndarray, n_mol: int) -> float:
    """Choose PME cutoff (nm) strictly below half the shortest box side.

    Orthorhombic LAMMPS cell: ``cell[i,i]`` in Ångström → nm is ``/10``; half-box is ``/20``.
    Previously we used 0.7 nm for ``n_mol > 32``, which is invalid for small ice supercells
    (e.g. 2×2×2 → ~0.45 nm half width) and triggers OpenMM's box-size check.
    """
    lx, ly, lz = float(cell[0, 0]), float(cell[1, 1]), float(cell[2, 2])
    half_min_nm = min(abs(lx), abs(ly), abs(lz)) / 20.0
    preferred = 0.4 if n_mol <= 64 else 0.7
    return min(preferred, 0.99 * half_min_nm)


def _centroid_kinetic_temperature_and_pe_kj_mol(
    integrator: RPMDIntegrator,
    n_beads: int,
    masses: np.ndarray,
    vel_buf: np.ndarray,
) -> tuple[float, float]:
    """Kinetic temperature (K) and PE (kJ/mol) consistent with UMA RPMD order CSV.

    ``T`` is from the centroid velocity (mean over beads), with equipartition using
    all particles with positive mass (same ``k_B`` convention as ``test_uma_ice_rpmd``).
    PE is from bead 0.
    """
    n_total = masses.shape[0]
    for b in range(n_beads):
        st = integrator.getState(b, getVelocities=True)
        vel_buf[b, :, :] = st.getVelocities(asNumpy=True).value_in_unit(
            unit.nanometer / unit.picosecond
        )
    centroid_vel = np.mean(vel_buf, axis=0)
    ke = 0.0
    for j in range(n_total):
        m = float(masses[j])
        if m <= 0.0:
            continue
        ke += 0.5 * m * float(np.sum(centroid_vel[j] ** 2))
    ndof = 3 * int(np.count_nonzero(masses > 0.0))
    if ndof <= 0:
        t_k = float("nan")
    else:
        gas_k = 8.314e-3  # kJ/(mol*K), m in dalton, v in nm/ps → kJ/mol for KE
        t_k = (2.0 * ke) / (ndof * gas_k)
    st0 = integrator.getState(0, getEnergy=True)
    pe = st0.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    return t_k, pe


def main() -> None:
    ap = argparse.ArgumentParser(description="OpenMM TIP4P/2005f RPMD ice simulation")
    ap.add_argument("--data", type=Path, default=None)
    ap.add_argument("--molecules", type=int, default=128)
    ap.add_argument("--beads", type=int, default=8)
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument(
        "--prod",
        type=float,
        default=1.0,
        help="Production time in ps (default: 1.0). Ignored if --steps is set.",
    )
    ap.add_argument(
        "--steps",
        type=int,
        default=None,
        help="Production steps (overrides --prod). prod_ps = steps * dt / 1000",
    )
    ap.add_argument("--temperature", type=float, default=243.0)
    ap.add_argument("--order-csv", type=Path, default=None)
    ap.add_argument("--order-every", type=int, default=10)
    ap.add_argument("--seed", type=int, default=284759)
    ap.add_argument("--platform", default="cuda")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--force-field", type=Path, default=None)
    ap.add_argument(
        "--cutoff-nm",
        type=float,
        default=None,
        help=(
            "Nonbonded cutoff (nm). Default: min(0.4 if molecules<=64 else 0.7, 0.99×half shortest box) "
            "so PME cutoff stays below half the periodic box."
        ),
    )
    ap.add_argument("--skip-minimize", action="store_true")
    ap.add_argument(
        "--rpmd-thermostat",
        type=str,
        default="pile-g",
        choices=["pile-g", "pile", "none"],
        help="pile-g = PILE_G (default); pile = PILE on all modes; none = NVE RPMD",
    )
    ap.add_argument(
        "--rpmd-friction",
        type=float,
        default=1.0,
        help="PILE internal friction 1/ps (default: 1.0)",
    )
    ap.add_argument(
        "--rpmd-centroid-friction",
        type=float,
        default=0.5,
        help="PILE_G Bussi centroid coupling 1/ps (default: 0.5)",
    )
    args = ap.parse_args()

    dt_fs = args.dt
    dt_ps = dt_fs / 1000.0
    if args.steps is not None:
        prod_ps = args.steps * dt_ps
    else:
        prod_ps = args.prod

    data_path = args.data or (_LAMMPS_DIR / f"data.ice_uma_{args.molecules}")
    if not data_path.is_file():
        m = args.molecules
        if m == 32:
            _nx, _ny, _nz = 1, 2, 2
        elif m == 64:
            _nx, _ny, _nz = 2, 2, 2
        elif m == 128:
            _nx, _ny, _nz = 2, 2, 4
        else:
            _nx = _ny = _nz = None
        if _nx is not None:
            hint = (
                f"python lammps/build_lammps_ice_data.py --nx {_nx} --ny {_ny} --nz {_nz} -o {data_path}"
            )
        else:
            hint = (
                f"python lammps/build_lammps_ice_data.py --nx NX --ny NY --nz NZ -o {data_path}  "
                f"(8×NX×NY×NZ = molecule count)"
            )
        print(f"Missing {data_path}. Run: {hint}", file=sys.stderr)
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

    cutoff_nm = (
        float(args.cutoff_nm)
        if args.cutoff_nm is not None
        else _auto_pme_cutoff_nm(cell, n_mol)
    )
    _half_min_nm = min(abs(cell[0, 0]), abs(cell[1, 1]), abs(cell[2, 2])) / 20.0
    print(
        f"  PME nonbonded cutoff: {cutoff_nm:.4f} nm (orthorhombic half min box: {_half_min_nm:.4f} nm)",
        flush=True,
    )
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=None,
        rigidWater=False,
        ewaldErrorTolerance=0.0005,
    )

    _ts = args.rpmd_thermostat.replace("-", "_")
    integrator = RPMDIntegrator(
        args.beads,
        args.temperature * unit.kelvin,
        args.rpmd_friction / unit.picosecond,
        dt_fs * unit.femtoseconds,
    )
    integrator.setRandomNumberSeed(args.seed)
    if _ts == "pile_g":
        integrator.setThermostatType(RPMDIntegrator.PileG)
        integrator.setCentroidFriction(args.rpmd_centroid_friction / unit.picosecond)
    elif _ts == "pile":
        integrator.setThermostatType(RPMDIntegrator.Pile)
    elif _ts == "none":
        integrator.setThermostatType(RPMDIntegrator.NoneThermo)
    else:
        raise ValueError(f"Invalid --rpmd-thermostat: {args.rpmd_thermostat!r}")

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

    total_steps = int(prod_ps / dt_ps)
    order_every = args.order_every or 10
    report_every = max(1, int(1.0 / dt_ps))

    box_orth_ang = np.linalg.norm(box_matrix, axis=1) * 10.0
    z_arr = np.array([8, 1, 1, 0] * n_mol, dtype=np.int64)

    masses = np.array(
        [system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(n_total)],
        dtype=np.float64,
    )
    vel_buf = np.zeros((args.beads, n_total, 3), dtype=np.float64)

    order_csv_f = None
    if args.order_csv:
        args.order_csv.parent.mkdir(parents=True, exist_ok=True)
        # Explicit flush after each row so monitoring with `tail -f` works during long runs.
        order_csv_f = open(args.order_csv, "w", encoding="utf-8", newline="\n")
        order_csv_f.write(ICE_ORDER_PI_CSV_HEADER)
        order_csv_f.flush()

    if ice_order_metrics_path_integral and format_ice_order_pi_csv_row and order_csv_f:
        bead_ang = []
        for b in range(args.beads):
            st = integrator.getState(b, getPositions=True)
            bead_ang.append(
                np.array(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)) * 10.0
            )
        om0 = ice_order_metrics_path_integral(bead_ang, box_orth_ang, z=z_arr)
        order_csv_f.write(format_ice_order_pi_csv_row(0, 0.0, "", "", om0))
        order_csv_f.flush()

    print(
        f"TIP4P/2005f RPMD: {n_mol} mol, {args.beads} beads, {dt_fs} fs, {prod_ps} ps, {total_steps} steps | "
        f"thermostat={args.rpmd_thermostat}"
    )
    start = time.time()
    for step in range(1, total_steps + 1):
        integrator.step(1)
        if step % report_every == 0:
            t_ps = step * dt_ps
            state = integrator.getState(0, getEnergy=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            ns_day = (step / (time.time() - start)) * dt_ps * 86400 / 1000.0
            print(f"  step {step:6d}  t={t_ps:.3f} ps  PE={pe:.1f}  {ns_day:.2f} ns/day")
        if order_csv_f and step % order_every == 0 and ice_order_metrics_path_integral:
            bead_ang = []
            for b in range(args.beads):
                s = integrator.getState(b, getPositions=True)
                bead_ang.append(
                    np.array(s.getPositions(asNumpy=True).value_in_unit(unit.nanometer)) * 10.0
                )
            om = ice_order_metrics_path_integral(bead_ang, box_orth_ang, z=z_arr)
            t_ps = step * dt_ps
            t_k_str, pe_str = "", ""
            t_k_val, pe_val = _centroid_kinetic_temperature_and_pe_kj_mol(
                integrator, args.beads, masses, vel_buf
            )
            if np.isfinite(t_k_val):
                t_k_str = f"{t_k_val:.4f}"
            pe_str = f"{pe_val:.6f}"
            order_csv_f.write(format_ice_order_pi_csv_row(step, t_ps, t_k_str, pe_str, om))
            order_csv_f.flush()

    if order_csv_f:
        order_csv_f.close()
        print(f"Order CSV: {args.order_csv}")
    print("Done.")


if __name__ == "__main__":
    main()
