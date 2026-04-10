#!/usr/bin/env python3
"""
OpenMM classical MD with TIP4P/2005f flexible water on the same ICE setup as UMA.

  Same initial structure: lammps/data.ice_uma (read_data)
  Force field: TIP4P/2005f (González & Abascal, J. Chem. Phys. 135, 224516, 2011)
  Minimize, then NVT Langevin 243 K. Default 0.1 fs (paper-recommended for TIP4P/2005f); 0.2 fs may be unstable.
  Outputs same order-parameter CSV schema (Q6, q_tet) as run_openmm_ice_lammps_match.py
  for comparison with UMA and LAMMPS.

Prerequisites:
  pip install openmm
  Build data:  python lammps/build_lammps_ice_data.py --nx 2 --ny 2 --nz 4 -o lammps/data.ice_uma

Usage:
  cd tests/uma_ice_rpmd
  python run_openmm_ice_classical_flex.py --steps 1000 --platform cuda
  python run_openmm_ice_classical_flex.py --steps 5000 --order-csv ice_order_classical.csv
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
    VariableLangevinIntegrator,
    unit,
)
from openmm.app import ForceField, Modeller, PDBFile, PME, Topology, element

sys.path.insert(0, str(_SCRIPT_DIR))
from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data

try:
    from ice_order_parameters import ice_order_metrics
except ImportError:
    ice_order_metrics = None  # type: ignore


def _max_force_norm(state) -> float:
    """Max |force| over all particles, in kJ/(mol·nm)."""
    forces = state.getForces()
    max_f = 0.0
    for f in forces:
        fx, fy, fz = f[0].value_in_unit(unit.kilojoules_per_mole / unit.nanometer), \
            f[1].value_in_unit(unit.kilojoules_per_mole / unit.nanometer), \
            f[2].value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        max_f = max(max_f, (fx*fx + fy*fy + fz*fz)**0.5)
    return max_f


def _tip4p2005f_forcefield_path() -> Path:
    """Resolve path to tip4p2005f.xml (repo or installed package)."""
    import openmm.app as _app
    pkg_data = Path(_app.__file__).parent / "data" / "tip4p2005f.xml"
    if pkg_data.is_file():
        return pkg_data
    repo_xml = _SCRIPT_DIR.parent.parent / "wrappers/python/openmm/app/data/tip4p2005f.xml"
    if repo_xml.is_file():
        return repo_xml
    raise FileNotFoundError(
        "tip4p2005f.xml not found. Use OpenMM built from this repo or pass --force-field path."
    )


# TIP4P/2005f equilibrium geometry (González & Abascal 2011): r_OH nm, H-O-H rad
_TIP4P2005F_R_OH_NM = 0.09419
_TIP4P2005F_ANGLE_RAD = 1.87448  # 107.4°


def _relax_water_geometry_to_tip4p2005f(pos_nm: np.ndarray, n_mol: int, box_nm: np.ndarray) -> np.ndarray:
    """
    Set each water's O-H and H-O-H to TIP4P/2005f equilibrium, keeping O and bisector direction.
    pos_nm: (n_total, 3) with n_total = 4*n_mol (O, H1, H2, M per molecule). In nm.
    Returns new positions (n_total, 3). Applies PBC when computing vectors.
    """
    out = np.array(pos_nm, dtype=np.float64, copy=True)
    r_eq = _TIP4P2005F_R_OH_NM
    half_angle = _TIP4P2005F_ANGLE_RAD / 2.0
    # Paper Table I: d_OM^rel = 0.13194 => average3 weights w2=w3=0.13194, w1=0.73612
    w1, w2 = 0.73612, 0.13194

    for m in range(n_mol):
        oi, h1i, h2i, mi = 4 * m, 4 * m + 1, 4 * m + 2, 4 * m + 3
        o = out[oi]
        h1 = out[h1i]
        h2 = out[h2i]
        # PBC wrap for vectors
        def wrap(d):
            d = np.array(d, dtype=np.float64)
            for axis in range(3):
                L = box_nm[axis]
                d[axis] -= np.round(d[axis] / L) * L
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
        if u_norm >= 1e-10:
            u_hat = u / u_norm
        else:
            u_hat = np.array([1.0, 0.0, 0.0]) if abs(b_hat[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            u_hat = u_hat - np.dot(u_hat, b_hat) * b_hat
            u_hat = u_hat / np.linalg.norm(u_hat)
        h1_new = o + r_eq * (np.cos(half_angle) * b_hat + np.sin(half_angle) * u_hat)
        h2_new = o + r_eq * (np.cos(half_angle) * b_hat - np.sin(half_angle) * u_hat)
        out[h1i] = h1_new
        out[h2i] = h2_new
        out[mi] = w1 * o + w2 * h1_new + w2 * h2_new
    return out


def _water_pdb_from_lammps(pos_ang: np.ndarray, cell: np.ndarray, n_mol: int) -> str:
    """Build PDB content (O,H1,H2 per residue) so Modeller.addExtraParticles works."""
    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    lines = [
        f"CRYST1{lx:9.3f}{ly:9.3f}{lz:9.3f}  90.00  90.00  90.00 P 1           1"
    ]
    for i in range(n_mol):
        o = pos_ang[3 * i]
        h1 = pos_ang[3 * i + 1]
        h2 = pos_ang[3 * i + 2]
        res = i + 1
        lines.append(f"HETATM{i*3+1:5d}  O   HOH  {res:4d}    {o[0]:8.3f}{o[1]:8.3f}{o[2]:8.3f}  1.00  0.00           O")
        lines.append(f"HETATM{i*3+2:5d}  H1  HOH  {res:4d}    {h1[0]:8.3f}{h1[1]:8.3f}{h1[2]:8.3f}  1.00  0.00           H")
        lines.append(f"HETATM{i*3+3:5d}  H2  HOH  {res:4d}    {h2[0]:8.3f}{h2[1]:8.3f}{h2[2]:8.3f}  1.00  0.00           H")
    lines.append("END")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="OpenMM TIP4P/2005f flexible water ice MD (same structure as UMA/LAMMPS)"
    )
    ap.add_argument("--data", type=Path, default=_LAMMPS_DIR / "data.ice_uma")
    ap.add_argument(
        "--force-field",
        type=Path,
        default=None,
        help="Path to tip4p2005f.xml (default: auto-detect from repo or package)",
    )
    ap.add_argument("--temperature", type=float, default=243.0)
    ap.add_argument("--friction", type=float, default=1.0, help="1/ps")
    ap.add_argument("--dt-fs", type=float, default=0.1, help="Timestep in fs (0.1 recommended for TIP4P/2005f; 0.2 may be unstable)")
    ap.add_argument("--steps", type=int, default=1000)
    ap.add_argument("--minimize-iters", type=int, default=2000, help="Max minimization iterations (default 2000 for flexible water)")
    ap.add_argument("--minimize-tol", type=float, default=10.0)
    ap.add_argument("--skip-minimize", action="store_true")
    ap.add_argument(
        "--relax-geometry",
        action="store_true",
        default=True,
        help="Set each water to TIP4P/2005f equilibrium (r_OH=0.09419 nm, H-O-H=107.4°) before minimize (default: True). Required when initial structure has rigid TIP4P/2005 geometry (e.g. ice CIF).",
    )
    ap.add_argument("--no-relax-geometry", action="store_false", dest="relax_geometry", help="Skip geometry relaxation (use if structure is already TIP4P/2005f-equilibrated).")
    ap.add_argument("--platform", default="cuda", choices=["cuda", "cpu"])
    ap.add_argument(
        "--variable-step",
        action="store_true",
        help="Use VariableLangevinIntegrator for stability (adapts dt; run until simulated time = steps*dt)",
    )
    ap.add_argument(
        "--variable-step-tol",
        type=float,
        default=0.01,
        help="Error tolerance for variable-step integrator (default 0.01)",
    )
    ap.add_argument(
        "--variable-step-chunk-fs",
        type=float,
        default=None,
        help="Max time chunk per stepTo() in fs for variable-step (default: same as --dt-fs). Use 0.1 for stricter stepping.",
    )
    ap.add_argument("--seed", type=int, default=284759)
    ap.add_argument("--cutoff", type=float, default=0.7, help="Nonbonded cutoff in nm (default 0.7; reduce for small boxes, e.g. 0.4 for 32 molecules)")
    ap.add_argument("-o", "--output", type=Path, default=None)
    ap.add_argument("--trajectory", type=Path, default=None)
    ap.add_argument("--traj-every", type=int, default=1)
    ap.add_argument("--xyz", type=Path, default=None)
    ap.add_argument("--no-trajectory", action="store_true")
    ap.add_argument("--report-every", type=int, default=10)
    ap.add_argument("--order-every", type=int, default=0)
    ap.add_argument("--order-csv", type=Path, default=None)
    args = ap.parse_args()

    if args.order_csv is not None and args.order_every <= 0:
        args.order_every = max(1, args.report_every)
        print(
            f"  NOTE: --order-csv without --order-every → using --order-every={args.order_every}",
            flush=True,
        )

    if not args.data.is_file():
        print(
            f"Missing {args.data} — run: python {_LAMMPS_DIR / 'build_lammps_ice_data.py'} "
            f"--nx 2 --ny 2 --nz 4 -o {args.data}"
        )
        sys.exit(1)

    if args.force_field is not None and not args.force_field.is_file():
        print(f"Force field file not found: {args.force_field}")
        sys.exit(1)
    ff_path = args.force_field if args.force_field is not None else _tip4p2005f_forcefield_path()

    pos_ang, cell, Z = parse_lammps_data(args.data)
    n_real = Z.shape[0]
    assert n_real % 3 == 0
    n_mol = n_real // 3
    for i in range(n_mol):
        assert Z[3 * i] == 8 and Z[3 * i + 1] == 1 and Z[3 * i + 2] == 1

    pdb_string = _water_pdb_from_lammps(pos_ang, cell, n_mol)
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
        f.write(pdb_string)
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

    forcefield = ForceField(str(ff_path))
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addExtraParticles(forcefield)
    n_total = modeller.topology.getNumAtoms()
    assert n_total == 4 * n_mol

    # Sanity-check structure: O,H,H,M order per residue; positions finite
    for res in modeller.topology.residues():
        names = [a.name for a in res.atoms()]
        assert names == ["O", "H1", "H2", "M"], f"Expected O,H1,H2,M per residue, got {names}"
    pos_check = modeller.positions
    pos_nm_arr = np.array(
        [[pos_check[i][j].value_in_unit(unit.nanometer) for j in range(3)] for i in range(n_total)]
    )
    assert not np.any(np.isnan(pos_nm_arr)) and not np.any(np.isinf(pos_nm_arr))
    print("  Structure check: O,H1,H2,M per residue OK; positions finite", flush=True)

    cutoff_nm = args.cutoff
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=None,
        rigidWater=False,
        ewaldErrorTolerance=0.0005,
    )
    # PBC sanity: cutoff must be < half the smallest box length to avoid spurious interactions.
    half_box_nm = min(float(lx), float(ly), float(lz)) / 20.0  # Angstrom -> nm, half
    assert cutoff_nm < half_box_nm, (
        f"PBC: nonbondedCutoff {cutoff_nm} nm must be < half smallest box {half_box_nm:.3f} nm. "
        "Increase box or reduce cutoff."
    )
    print(f"  PBC check: cutoff {cutoff_nm} nm < half_box {half_box_nm:.3f} nm OK", flush=True)

    if args.variable_step:
        # Variable-step Langevin for stability (adapts dt to keep error below tolerance).
        integrator = VariableLangevinIntegrator(
            args.temperature * unit.kelvin,
            args.friction / unit.picosecond,
            args.variable_step_tol,
        )
        # OpenMM expects max step size in picoseconds; 0.2 fs = 0.0002 ps
        integrator.setMaximumStepSize((args.dt_fs / 1000.0) * unit.picoseconds)
        target_time_ps = args.steps * (args.dt_fs / 1000.0)
        print(f"  Integrator: VariableLangevin (errorTol={args.variable_step_tol}, maxStep={args.dt_fs} fs)")
        max_step = integrator.getMaximumStepSize()
        max_step_ps = max_step.value_in_unit(unit.picoseconds) if hasattr(max_step, "value_in_unit") else float(max_step)
        print(f"  getMaximumStepSize() = {max_step} (0 => not applied)", flush=True)
        if max_step_ps <= 0.0 or max_step_ps > 0.001:
            print(
                "  WARNING: getMaximumStepSize() is 0 or >1 fs; variable-step may use large internal steps. "
                "Try --variable-step-chunk-fs 0.1 to force 0.1 fs advances.",
                flush=True,
            )
        print(f"  Run until simulated time >= {target_time_ps:.4f} ps", flush=True)
    else:
        integrator = LangevinMiddleIntegrator(
            args.temperature * unit.kelvin,
            args.friction / unit.picosecond,
            args.dt_fs * unit.femtoseconds,
        )
        target_time_ps = None
    integrator.setRandomNumberSeed(args.seed)

    platform = Platform.getPlatformByName(args.platform.upper())
    props = {}
    if args.platform == "cuda":
        props = {"Precision": "mixed", "DeviceIndex": "0", "DisablePmeStream": "true"}

    context = Context(system, integrator, platform, props)
    context.setPeriodicBoxVectors(a, b, c)
    context.setPositions(modeller.positions)

    # Unit check: positions and box in expected range (nm)
    box_nm = context.getState(getPositions=False).getPeriodicBoxVectors()
    def _vec_norm_nm(v):
        return (v[0]**2 + v[1]**2 + v[2]**2)**0.5
    box_lens_nm = np.array([
        _vec_norm_nm(box_nm[0]).value_in_unit(unit.nanometer),
        _vec_norm_nm(box_nm[1]).value_in_unit(unit.nanometer),
        _vec_norm_nm(box_nm[2]).value_in_unit(unit.nanometer),
    ])
    for i, L_val in enumerate(box_lens_nm):
        assert 0.5 < L_val < 100.0, f"Box vector {i} length {L_val} nm unexpected (expect ~1.5 nm)"
    pos0 = context.getState(getPositions=True).getPositions()
    # First two O positions (indices 0 and 4); expect separation ~0.28 nm
    o0 = np.array([pos0[0][j].value_in_unit(unit.nanometer) for j in range(3)])
    o1 = np.array([pos0[4][j].value_in_unit(unit.nanometer) for j in range(3)])
    dr = o1 - o0
    dr = dr - np.round(dr / box_lens_nm) * box_lens_nm
    oo_dist = float((dr**2).sum()**0.5)
    assert 0.2 < oo_dist < 1.0, f"O-O distance {oo_dist:.4f} nm unexpected (expect ~0.28 nm)"
    print(f"  Unit check: box ~{box_lens_nm[0]:.3f} nm, O-O(0,1)={oo_dist:.3f} nm OK", flush=True)

    if getattr(args, "relax_geometry", False):
        pos0_arr = np.array(
            [[pos0[i][j].value_in_unit(unit.nanometer) for j in range(3)] for i in range(n_total)]
        )
        pos_relaxed = _relax_water_geometry_to_tip4p2005f(pos0_arr, n_mol, box_lens_nm)
        relaxed_positions = [Vec3(pos_relaxed[i][0], pos_relaxed[i][1], pos_relaxed[i][2]) * unit.nanometer for i in range(n_total)]
        context.setPositions(relaxed_positions)
        print("  Applied --relax-geometry: each water set to r_OH=0.09419 nm, H-O-H=107.4°", flush=True)

    print("=" * 72)
    print("OpenMM TIP4P/2005f flexible water — ICE benchmark (same structure as UMA)")
    print("=" * 72)
    print(f"  Molecules: {n_mol}  atoms: {n_total} (3 O,H,H + 1 M per water)")
    print(f"  T = {args.temperature} K  dt = {args.dt_fs} fs  platform: {args.platform}")

    # Pre-minimization geometry (TIP4P/2005f equilibrium: r_OH=0.09419 nm, H-O-H=107.4°)
    def _pos_nm(pos_list, i, j=3):
        return np.array([pos_list[i][k].value_in_unit(unit.nanometer) for k in range(j)])
    def _dist_nm(pos_list, a, b, box_nm=None):
        d = _pos_nm(pos_list, b) - _pos_nm(pos_list, a)
        if box_nm is not None:
            for axis in range(3):
                L = box_nm[axis]
                d[axis] -= np.round(d[axis] / L) * L
        return float(np.sqrt((d**2).sum()))
    def _angle_rad(pos_list, o, h1, h2):
        u = _pos_nm(pos_list, h1) - _pos_nm(pos_list, o)
        v = _pos_nm(pos_list, h2) - _pos_nm(pos_list, o)
        u = u - np.round(u / box_lens_nm) * box_lens_nm
        v = v - np.round(v / box_lens_nm) * box_lens_nm
        nu, nv = np.linalg.norm(u), np.linalg.norm(v)
        if nu < 1e-10 or nv < 1e-10:
            return float("nan")
        return np.arccos(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    pos_pre = context.getState(getPositions=True).getPositions()
    oh_vals, angle_vals = [], []
    for m in range(min(10, n_mol)):
        oi, h1i, h2i = 4 * m, 4 * m + 1, 4 * m + 2
        oh_vals.append(_dist_nm(pos_pre, oi, h1i, box_lens_nm))
        oh_vals.append(_dist_nm(pos_pre, oi, h2i, box_lens_nm))
        angle_vals.append(_angle_rad(pos_pre, oi, h1i, h2i) * 180.0 / np.pi)
    min_oo_pre = 1e9
    for i in range(min(32, n_mol)):
        for j in range(i + 1, min(32, n_mol)):
            min_oo_pre = min(min_oo_pre, _dist_nm(pos_pre, 4 * i, 4 * j, box_lens_nm))
    st_pre = context.getState(getEnergy=True, getForces=True)
    pe_pre = st_pre.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    max_f_pre = _max_force_norm(st_pre)
    print(f"  [pre-min] PE={pe_pre:.1f} kJ/mol  max|F|={max_f_pre:.1f} kJ/(mol·nm)", flush=True)
    print(f"  [pre-min] O-H (first 10 mol): min={min(oh_vals):.4f} max={max(oh_vals):.4f} nm (TIP4P/2005f eq 0.09419)", flush=True)
    print(f"  [pre-min] H-O-H (first 10 mol): {np.nanmin(angle_vals):.1f}-{np.nanmax(angle_vals):.1f} deg (eq 107.4)", flush=True)
    print(f"  [pre-min] min O-O (first 32): {min_oo_pre:.4f} nm", flush=True)

    if not args.skip_minimize:
        t0 = time.perf_counter()
        LocalEnergyMinimizer.minimize(
            context, maxIterations=args.minimize_iters, tolerance=args.minimize_tol
        )
        st = context.getState(getEnergy=True, getForces=True, getPositions=True)
        pe = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"  PE after minimize: {pe:.2f} kJ/mol")
        max_f = _max_force_norm(st)
        print(f"  Max |force| after minimize: {max_f:.2f} kJ/(mol·nm)", flush=True)
        max_force_threshold = 1e5  # kJ/(mol·nm); sane runs are typically < 1e4
        if max_f > max_force_threshold:
            print(
                f"  ERROR: max |force| {max_f:.2f} exceeds threshold {max_force_threshold:.0f}. "
                "Force field matches González & Abascal J. Chem. Phys. 135, 224516 (2011); "
                "see check_tip4p2005f_single_water.py. If pre-min O-H / H-O-H are far from 0.09419 nm and 107.4°, "
                "the initial structure (e.g. ice CIF) may be incompatible with TIP4P/2005f equilibrium; "
                "equilibrate with this model first or use a structure already relaxed with TIP4P/2005f.",
                flush=True,
            )
            sys.exit(1)
        pos_nm = st.getPositions()
        # O-H distances first molecule (O=0, H1=1, H2=2)
        def _dist(a, b):
            d = np.array([pos_nm[b][j].value_in_unit(unit.nanometer) - pos_nm[a][j].value_in_unit(unit.nanometer) for j in range(3)])
            return float((d**2).sum()**0.5)
        oh1, oh2 = _dist(0, 1), _dist(0, 2)
        print(f"  O-H distances (first mol): {oh1:.4f}, {oh2:.4f} nm (expect ~0.096)", flush=True)
        # Min O-O over first 20 O's (indices 0,4,...,76)
        o_inds = [4 * i for i in range(min(20, n_mol))]
        min_oo = 1e9
        for i, ii in enumerate(o_inds):
            for jj in o_inds[i+1:]:
                min_oo = min(min_oo, _dist(ii, jj))
        print(f"  Min O-O (first 20): {min_oo:.4f} nm (expect > 0.25)", flush=True)
        print(f"  Minimize wall time: {time.perf_counter() - t0:.1f} s", flush=True)
    else:
        print("  --skip-minimize: no minimization.")
        # Diagnostic: initial max force (unminimized structure)
        st0 = context.getState(getEnergy=True, getForces=True)
        print(f"  [diagnostic] Initial (unminimized) PE={st0.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.2f} kJ/mol  max|F|={_max_force_norm(st0):.2f} kJ/(mol·nm)", flush=True)

    context.setVelocitiesToTemperature(args.temperature * unit.kelvin, args.seed)

    # Diagnostic: KE and T after velocity init (use real-atom DOF for T)
    dof_real = 3 * n_real
    R = unit.MOLAR_GAS_CONSTANT_R
    st_init = context.getState(getEnergy=True)
    ke_init = st_init.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
    T_init = (2 * st_init.getKineticEnergy() / (dof_real * R)).value_in_unit(unit.kelvin)
    print(f"  [diagnostic] After setVelocitiesToTemperature: KE={ke_init:.2f} kJ/mol  T={T_init:.2f} K", flush=True)

    modeller.topology.setPeriodicBoxVectors((a, b, c))
    out = args.output or (
        _SCRIPT_DIR / f"ice_classical_flex_T{args.temperature:.0f}_steps{args.steps}.pdb"
    )
    if args.no_trajectory:
        dcd_path = None
    elif args.trajectory is not None:
        dcd_path = args.trajectory
    else:
        dcd_path = out.with_suffix(".dcd")

    dcd = None
    xyz_f = None
    symbols_4 = ["O", "H", "H", "M"] * n_mol
    if dcd_path:
        from openmm.app import DCDFile

        dcd_f = open(dcd_path, "wb")
        dcd = DCDFile(
            dcd_f,
            modeller.topology,
            args.dt_fs * unit.femtoseconds,
            interval=args.traj_every,
        )
        print(f"  Trajectory DCD: {dcd_path} (every {args.traj_every} step(s))")
    if args.xyz:
        xyz_f = open(args.xyz, "w")
        print(f"  Trajectory XYZ: {args.xyz}")

    _uc_nm = Vec3(float(lx) * 0.1, float(ly) * 0.1, float(lz) * 0.1) * unit.nanometer
    box_ang = np.array([float(lx), float(ly), float(lz)], dtype=np.float64)
    z_4 = np.array([8, 1, 1, 0] * n_mol, dtype=np.int64)

    def _write_traj_frame(step_index: int, pos_nm) -> None:
        if dcd is not None and (step_index == 0 or step_index % args.traj_every == 0):
            dcd.writeModel(pos_nm, unitCellDimensions=_uc_nm)
        if xyz_f is not None and (step_index == 0 or step_index % args.traj_every == 0):
            pos_a = np.array(
                [[p[i].value_in_unit(unit.nanometer) * 10.0 for i in range(3)] for p in pos_nm],
                dtype=np.float64,
            )
            xyz_f.write(f"{n_total}\n")
            xyz_f.write(f"Step {step_index} OpenMM TIP4P/2005f\n")
            for i in range(n_total):
                xyz_f.write(
                    f"{symbols_4[i]:2s} {pos_a[i,0]:.8f} {pos_a[i,1]:.8f} {pos_a[i,2]:.8f}\n"
                )
            xyz_f.flush()

    st0 = context.getState(getPositions=True)
    _write_traj_frame(0, st0.getPositions())

    order_csv_f = None
    if args.order_every > 0:
        if ice_order_metrics is None:
            print("  WARNING: --order-every ignored (ice_order_parameters import failed)", flush=True)
        else:
            print(
                f"  Order params every {args.order_every} steps: Q6, q_tet",
                flush=True,
            )
            if args.order_csv:
                order_csv_f = open(args.order_csv, "w")
                order_csv_f.write("step,time_ps,T_K,PE_kj_mol,q6_mean,q6_std,q_tet_mean,q_tet_std\n")

    dt_ps = args.dt_fs / 1000.0
    # Virtual sites (M) have zero mass; getKineticEnergy() includes only real atoms → DOF = 3*n_real
    dof = 3 * n_real
    R = unit.MOLAR_GAS_CONSTANT_R

    def _T_K(st) -> float:
        ke = st.getKineticEnergy()
        return (2 * ke / (dof * R)).value_in_unit(unit.kelvin)

    use_variable_step = args.variable_step and target_time_ps is not None
    if use_variable_step:
        print(f"\n--- MD until time >= {target_time_ps:.4f} ps (variable step) ---")
        max_steps_guard = max(args.steps * 3, 100_000)
    else:
        print(f"\n--- MD {args.steps} steps @ {args.dt_fs} fs ---")
    t_md = time.perf_counter()
    write_traj = dcd is not None or xyz_f is not None

    step = 0
    last_t_ps = 0.0

    if use_variable_step:
        # Drive by time chunks so kernel gets finite maxTime. Optionally use smaller chunk via --variable-step-chunk-fs.
        chunk_fs = args.variable_step_chunk_fs if args.variable_step_chunk_fs is not None else args.dt_fs
        chunk_ps = chunk_fs / 1000.0
        n_samples = max(1, int(round(target_time_ps / chunk_ps)))
        sample_times_ps = np.linspace(chunk_ps, target_time_ps, num=n_samples)
        if chunk_fs != args.dt_fs:
            print(f"  Using variable-step chunk: {chunk_fs} fs ({chunk_ps:.6f} ps) -> {n_samples} samples", flush=True)
        for step_index, t_ps in enumerate(sample_times_ps, start=1):
            integrator.stepTo(t_ps * unit.picoseconds)
            t_ps_now = context.getTime().value_in_unit(unit.picosecond)
            last_t_ps = t_ps_now
            if write_traj and (step_index == 1 or step_index % args.traj_every == 0):
                stp = context.getState(getPositions=True)
                _write_traj_frame(step_index, stp.getPositions())
            if step_index % args.report_every == 0 or step_index == n_samples:
                st = context.getState(getEnergy=True, getVelocities=True)
                pe = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                temp = _T_K(st)
                if step_index == 1:
                    ke1 = st.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
                    print(f"  [diagnostic] After 1st step: PE={pe:.2f} kJ/mol  KE={ke1:.2f} kJ/mol  T={temp:.2f} K", flush=True)
                line = f"  step {step_index:5d}  t={t_ps_now:.4f} ps  T={temp:7.2f} K  PE={pe:14.2f} kJ/mol"
                if args.order_every > 0 and ice_order_metrics is not None and (step_index % args.order_every == 0 or step_index == n_samples):
                    pos_nm = context.getState(getPositions=True).getPositions()
                    pos_ang_all = np.array(
                        [[p[i].value_in_unit(unit.nanometer) * 10.0 for i in range(3)] for p in pos_nm],
                        dtype=np.float64,
                    )
                    om = ice_order_metrics(pos_ang_all, box_ang, z=z_4)
                    line += f"  Q6={om.q6_mean:.3f}  q_tet={om.q_tet_mean:.3f}"
                    if order_csv_f is not None:
                        order_csv_f.write(
                            f"{step_index},{t_ps_now:.6f},{temp:.4f},{pe:.6f},"
                            f"{om.q6_mean:.6f},{om.q6_std:.6f},{om.q_tet_mean:.6f},{om.q_tet_std:.6f}\n"
                        )
                        order_csv_f.flush()
                print(line)
    else:
        while True:
            step += 1
            integrator.step(1)
            t_ps = step * dt_ps
            last_t_ps = t_ps
            if step >= args.steps:
                break
            if write_traj and step % args.traj_every == 0:
                stp = context.getState(getPositions=True)
                _write_traj_frame(step, stp.getPositions())
            if step % args.report_every == 0 or step == args.steps:
                st = context.getState(getEnergy=True, getVelocities=True)
                pe = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                temp = _T_K(st)
                if step == 1:
                    ke1 = st.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
                    print(f"  [diagnostic] After 1st step: PE={pe:.2f} kJ/mol  KE={ke1:.2f} kJ/mol  T={temp:.2f} K", flush=True)
                line = f"  step {step:5d}  t={t_ps:.4f} ps  T={temp:7.2f} K  PE={pe:14.2f} kJ/mol"
                if (
                    args.order_every > 0
                    and ice_order_metrics is not None
                    and (step % args.order_every == 0 or step == args.steps)
                ):
                    pos_nm = context.getState(getPositions=True).getPositions()
                    pos_ang_all = np.array(
                        [[p[i].value_in_unit(unit.nanometer) * 10.0 for i in range(3)] for p in pos_nm],
                        dtype=np.float64,
                    )
                    om = ice_order_metrics(pos_ang_all, box_ang, z=z_4)
                    line += f"  Q6={om.q6_mean:.3f}  q_tet={om.q_tet_mean:.3f}"
                    if order_csv_f is not None:
                        order_csv_f.write(
                            f"{step},{t_ps:.6f},{temp:.4f},{pe:.6f},"
                            f"{om.q6_mean:.6f},{om.q6_std:.6f},{om.q_tet_mean:.6f},{om.q_tet_std:.6f}\n"
                        )
                        order_csv_f.flush()
                print(line)

    elapsed = time.perf_counter() - t_md
    ps_done = last_t_ps
    ns = ps_done * 1e-3
    ns_per_day = ns / elapsed * 86400.0 if elapsed > 0 else 0.0
    n_steps_done = (n_samples if use_variable_step else step)
    steps_per_s = n_steps_done / elapsed if elapsed > 0 else 0.0
    print(f"  MD wall time: {elapsed:.1f} s")
    print(f"  Speed: {ns_per_day:.3f} ns/day  |  {steps_per_s:.2f} steps/s")

    if dcd is not None:
        dcd_f.close()
        print(f"  Trajectory (DCD): {dcd_path}")
    if xyz_f is not None:
        xyz_f.close()
        print(f"  Trajectory (XYZ): {args.xyz}")
    if order_csv_f is not None:
        order_csv_f.close()
        print(f"  Order CSV: {args.order_csv}")

    modeller.topology.setUnitCellDimensions(
        Vec3(float(lx) * 0.1, float(ly) * 0.1, float(lz) * 0.1) * unit.nanometer
    )
    PDBFile.writeFile(
        modeller.topology,
        context.getState(getPositions=True).getPositions(),
        open(out, "w"),
        keepIds=True,
    )
    print(f"  Final PDB: {out}")
    print("=" * 72)


if __name__ == "__main__":
    main()
