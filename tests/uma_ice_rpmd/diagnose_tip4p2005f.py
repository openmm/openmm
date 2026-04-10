#!/usr/bin/env python3
"""
Systematic diagnostics for TIP4P/2005f multi-molecule energy/force anomalies.

Steps (per plan):
  1. Energy decomposition by force group (Morse / angle / nonbonded)
  2. Nonbonded exception audit for intramolecular pairs (incl. M-site)
  3. NoCutoff vs PME total PE
  4. Reference vs CUDA (same geometry)
  5. Scaling PE/N vs molecule count (rebuild LAMMPS data per N)
  6. Optional: TIP4P-Ew rigid sanity check on same box

Usage:
  cd tests/uma_ice_rpmd
  python diagnose_tip4p2005f.py
  python diagnose_tip4p2005f.py --data lammps/data.ice_uma_64 --skip-scaling
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent
_LAMMPS_DIR = _SCRIPT_DIR / "lammps"
_ROOT = _SCRIPT_DIR.parent.parent

if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data
from openmm import Context, LangevinMiddleIntegrator, NonbondedForce, Platform, Vec3, unit
from openmm.app import ForceField, Modeller, PDBFile, PME, NoCutoff

# Reuse classical runner helpers
import run_openmm_ice_classical_flex as flex


def _tip4p2005f_path() -> Path:
    return flex._tip4p2005f_forcefield_path()


def _build_modeller_from_data(data_path: Path):
    pos_ang, cell, Z = parse_lammps_data(data_path)
    n_real = Z.shape[0]
    assert n_real % 3 == 0
    n_mol = n_real // 3
    pdb_string = flex._water_pdb_from_lammps(pos_ang, cell, n_mol)
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
    ff_path = _tip4p2005f_path()
    forcefield = ForceField(str(ff_path))
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addExtraParticles(forcefield)
    n_total = modeller.topology.getNumAtoms()
    assert n_total == 4 * n_mol
    return forcefield, modeller, (a, b, c), n_mol


def _vec_len_nm(v: Vec3) -> float:
    comps = v.value_in_unit(unit.nanometer)
    return float((comps[0] ** 2 + comps[1] ** 2 + comps[2] ** 2) ** 0.5)


def _positions_with_relax(modeller, n_mol: int, box_vecs) -> list:
    """Apply default geometry relaxation (same as production)."""
    a, b, c = box_vecs
    box_lens_nm = np.array([_vec_len_nm(a), _vec_len_nm(b), _vec_len_nm(c)])
    pos_check = modeller.positions
    pos_nm_arr = np.array(
        [[pos_check[i][j].value_in_unit(unit.nanometer) for j in range(3)] for i in range(4 * n_mol)]
    )
    pos_relaxed = flex._relax_water_geometry_to_tip4p2005f(pos_nm_arr, n_mol, box_lens_nm)
    return [
        Vec3(pos_relaxed[i][0], pos_relaxed[i][1], pos_relaxed[i][2]) * unit.nanometer
        for i in range(4 * n_mol)
    ]


def _create_system_typed(
    forcefield: ForceField,
    topology,
    *,
    nonbonded_method,
    cutoff_nm: float | None,
    ewald_tol: float = 0.0005,
):
    kwargs = dict(constraints=None, rigidWater=False, ewaldErrorTolerance=ewald_tol)
    if nonbonded_method is NoCutoff:
        system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff, **kwargs)
    else:
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=PME,
            nonbondedCutoff=cutoff_nm * unit.nanometer,
            **kwargs,
        )
    return system


def _assign_force_groups(system) -> dict[str, int]:
    """Assign every force a distinct group so summed group energies match total PE."""
    names: dict[str, int] = {}
    for fi in range(system.getNumForces()):
        f = system.getForce(fi)
        f.setForceGroup(fi)
        names[f"{type(f).__name__}[{fi}]"] = fi
    return names


def _energy_by_group(context: Context, n_groups: int) -> dict[int, float]:
    out = {}
    for g in range(n_groups):
        st = context.getState(getEnergy=True, groups={1 << g})
        out[g] = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    return out


def _sum_group_energies(context: Context, n_forces: int) -> float:
    return sum(
        context.getState(getEnergy=True, groups={1 << g}).getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole
        )
        for g in range(n_forces)
    )


def _max_force_norm(state) -> float:
    forces = state.getForces()
    max_f = 0.0
    for f in forces:
        fx = f[0].value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        fy = f[1].value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        fz = f[2].value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        max_f = max(max_f, (fx * fx + fy * fy + fz * fz) ** 0.5)
    return max_f


def _audit_nb_exceptions(system, n_mol: int) -> tuple[list[tuple[int, int]], list[str]]:
    nb = None
    for fi in range(system.getNumForces()):
        f = system.getForce(fi)
        if isinstance(f, NonbondedForce):
            nb = f
            break
    assert nb is not None
    excluded_pairs: set[tuple[int, int]] = set()
    for i in range(nb.getNumExceptions()):
        p1, p2, q, sig, eps = nb.getExceptionParameters(i)
        i1, i2 = min(p1, p2), max(p1, p2)
        excluded_pairs.add((i1, i2))
    missing: list[str] = []
    for m in range(n_mol):
        o, h1, h2, ms = 4 * m, 4 * m + 1, 4 * m + 2, 4 * m + 3
        for a, b, label in [
            (o, h1, "O-H1"),
            (o, h2, "O-H2"),
            (h1, h2, "H1-H2"),
            (ms, o, "M-O"),
            (ms, h1, "M-H1"),
            (ms, h2, "M-H2"),
        ]:
            pair = (min(a, b), max(a, b))
            if pair not in excluded_pairs:
                missing.append(f"mol{m} {label} pair ({pair[0]},{pair[1]})")
    return list(excluded_pairs), missing


def _run_context(
    system,
    positions,
    box_vecs,
    platform_name: str,
    precision: str | None = None,
):
    a, b, c = box_vecs
    integrator = LangevinMiddleIntegrator(
        243 * unit.kelvin, 10.0 / unit.picosecond, 0.0001 * unit.picoseconds
    )
    platform = Platform.getPlatformByName(platform_name)
    props = {}
    if platform_name.upper() == "CUDA":
        props["Precision"] = precision or "mixed"
        props["DeviceIndex"] = "0"
    context = Context(system, integrator, platform, props)
    context.setPeriodicBoxVectors(a, b, c)
    context.setPositions(positions)
    return context


def _pe_total(context: Context) -> float:
    st = context.getState(getEnergy=True)
    return st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)


def _min_oo_angstrom(pos_ang: np.ndarray, cell: np.ndarray, n_mol: int) -> float:
    """Minimum O–O distance (Å) with orthorhombic MIC; O at indices 3*i."""
    box = np.array([[cell[0, 0], 0.0, 0.0], [0.0, cell[1, 1], 0.0], [0.0, 0.0, cell[2, 2]]])

    def mic(dr: np.ndarray) -> np.ndarray:
        out = dr.copy()
        for k in range(3):
            L = box[k, k]
            out[k] -= np.round(out[k] / L) * L
        return out

    o_idx = [3 * i for i in range(n_mol)]
    min_d = 1e9
    for ii, i in enumerate(o_idx):
        for j in o_idx[ii + 1 :]:
            d = np.linalg.norm(mic(pos_ang[j] - pos_ang[i]))
            min_d = min(min_d, d)
    return float(min_d)


def _half_box_nm(box_vecs) -> float:
    a, b, c = box_vecs
    return 0.5 * min(_vec_len_nm(a), _vec_len_nm(b), _vec_len_nm(c))


def _safe_cutoff(box_vecs: tuple, default: float = 0.7) -> float:
    half = _half_box_nm(box_vecs)
    return min(default, half * 0.95)


def diagnose_data_file(data_path: Path, *, relax: bool = True) -> None:
    print("=" * 72)
    print(f"Data: {data_path}")
    print("=" * 72)
    pos_ang_pre, cell_pre, Z_pre = parse_lammps_data(data_path)
    n_mol_pre = Z_pre.shape[0] // 3
    d_oo = _min_oo_angstrom(pos_ang_pre, cell_pre, n_mol_pre)
    print(
        f"  Structure: min O–O (MIC) = {d_oo:.3f} Å  "
        f"(~1.7 Å ⇒ bad CIF truncation; ~2.2+ Å often OK for small ice Ic cells)"
    )
    forcefield, modeller, box_vecs, n_mol = _build_modeller_from_data(data_path)
    if relax:
        positions = _positions_with_relax(modeller, n_mol, box_vecs)
    else:
        positions = modeller.positions
    cutoff = _safe_cutoff(box_vecs)

    # --- Step 2: exclusion audit on PME system ---
    sys_pme = _create_system_typed(
        forcefield, modeller.topology, nonbonded_method=PME, cutoff_nm=cutoff
    )
    pairs, missing = _audit_nb_exceptions(sys_pme, n_mol)
    print(f"\n[Exceptions] NonbondedForce has {len(pairs)} exception pairs.")
    if missing:
        print(f"  MISSING intramolecular exclusions ({len(missing)}):")
        for line in missing[:20]:
            print(f"    {line}")
        if len(missing) > 20:
            print(f"    ... and {len(missing) - 20} more")
    else:
        print("  All O-H / H-H / M-* intramolecular pairs present as exceptions.")

    # --- Step 1: energy decomposition (PME) ---
    group_map = _assign_force_groups(sys_pme)
    n_grp = sys_pme.getNumForces()
    ctx = _run_context(sys_pme, positions, box_vecs, "Reference")
    eg = _energy_by_group(ctx, n_groups=n_grp)
    print("\n[Energy decomposition] PME, Reference platform, groups:")
    print(
        "  (Note: NonbondedForce+PME may report 0 on a single group on some builds; "
        "trust full PE.)"
    )
    for name, gid in sorted(group_map.items(), key=lambda x: x[1]):
        print(f"  {name}: {eg[gid]:.2f} kJ/mol")
    total_ref = _sum_group_energies(ctx, n_grp)
    total_full = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoules_per_mole
    )
    print(f"  Sum(groups): {total_ref:.2f} kJ/mol  full PE: {total_full:.2f} kJ/mol")
    print(f"  PE/mol ≈ {total_full / n_mol:.2f} kJ/mol")

    st = ctx.getState(getEnergy=True, getForces=True)
    print(f"  max|F| = {_max_force_norm(st):.2f} kJ/(mol·nm)")

    # --- Step 3: NoCutoff ---
    sys_nc = _create_system_typed(
        forcefield, modeller.topology, nonbonded_method=NoCutoff, cutoff_nm=None
    )
    ctx_nc = _run_context(sys_nc, positions, box_vecs, "Reference")
    pe_nc = _pe_total(ctx_nc)
    print(f"\n[NoCutoff vs PME] NoCutoff total PE = {pe_nc:.2f} kJ/mol")
    print(f"  PME total (full)             = {total_full:.2f} kJ/mol")

    # --- Step 4: CUDA if available ---
    try:
        sys_pme2 = _create_system_typed(
            forcefield, modeller.topology, nonbonded_method=PME, cutoff_nm=cutoff
        )
        ctx_cuda = _run_context(sys_pme2, positions, box_vecs, "CUDA", "mixed")
        pe_cuda_mixed = _pe_total(ctx_cuda)
        sys_pme3 = _create_system_typed(
            forcefield, modeller.topology, nonbonded_method=PME, cutoff_nm=cutoff
        )
        ctx_cuda_d = _run_context(sys_pme3, positions, box_vecs, "CUDA", "double")
        pe_cuda_double = _pe_total(ctx_cuda_d)
        print("\n[Platform] PME total PE:")
        print(f"  Reference:  {total_full:.2f} kJ/mol")
        print(f"  CUDA mixed: {pe_cuda_mixed:.2f} kJ/mol")
        print(f"  CUDA double:{pe_cuda_double:.2f} kJ/mol")
    except Exception as e:
        print(f"\n[Platform] CUDA skip: {e}")


def _build_data_for_n(n: int, out_path: Path) -> None:
    build = _LAMMPS_DIR / "build_lammps_ice_data.py"
    subprocess.run(
        [sys.executable, str(build), "-n", str(n), "-o", str(out_path)],
        cwd=str(_SCRIPT_DIR),
        check=True,
    )


def scaling_test(ns: list[int]) -> None:
    print("\n" + "=" * 72)
    print("Scaling: PE / N vs molecule count")
    print("=" * 72)
    ff_path = _tip4p2005f_path()
    forcefield = ForceField(str(ff_path))
    with tempfile.TemporaryDirectory() as td:
        tdir = Path(td)
        for n in ns:
            dp = tdir / f"data_{n}"
            try:
                _build_data_for_n(n, dp)
            except subprocess.CalledProcessError as e:
                print(f"  N={n}: build failed: {e}")
                continue
            _, modeller, box_vecs, n_mol = _build_modeller_from_data(dp)
            assert n_mol == n
            positions = _positions_with_relax(modeller, n_mol, box_vecs)
            cutoff = _safe_cutoff(box_vecs)
            sys_pme = _create_system_typed(
                forcefield, modeller.topology, nonbonded_method=PME, cutoff_nm=cutoff
            )
            ctx = _run_context(sys_pme, positions, box_vecs, "Reference")
            pe = _pe_total(ctx)
            print(f"  N={n:3d}: PE={pe:12.2f} kJ/mol  PE/N={pe / n:10.2f} kJ/mol")


def tip4pew_sanity(data_path: Path) -> None:
    """Rigid TIP4P-Ew on same PDB topology (3 atoms/residue, no M in PDB)."""
    print("\n" + "=" * 72)
    print("Reference: tip4pew.xml rigid water, PME (sanity)")
    print("=" * 72)
    pos_ang, cell, Z = parse_lammps_data(data_path)
    n_mol = Z.shape[0] // 3
    pdb_string = flex._water_pdb_from_lammps(pos_ang, cell, n_mol)
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
    ff = ForceField("tip4pew.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addExtraParticles(ff)
    cutoff = _safe_cutoff((a, b, c))
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=cutoff * unit.nanometer,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    ctx = _run_context(system, modeller.positions, (a, b, c), "Reference")
    pe = _pe_total(ctx)
    if pe != pe:  # NaN
        print("  TIP4P-Ew rigid PE = NaN (skip: check box/cutoff or platform)")
    else:
        print(f"  TIP4P-Ew rigid PE = {pe:.2f} kJ/mol  PE/N = {pe / n_mol:.2f} kJ/mol")


def main() -> int:
    if os.getenv("OPENMM_PLUGIN_DIR") is None:
        for _base in [os.getenv("CONDA_PREFIX"), os.path.expanduser("~/miniconda3")]:
            if _base:
                _plugins = os.path.join(_base, "lib", "plugins")
                if os.path.isdir(_plugins):
                    os.environ["OPENMM_PLUGIN_DIR"] = _plugins
                    break

    ap = argparse.ArgumentParser(description="Diagnose TIP4P/2005f multi-molecule energies")
    ap.add_argument(
        "--data",
        type=Path,
        nargs="*",
        default=[
            _LAMMPS_DIR / "data.ice_uma",
            _LAMMPS_DIR / "data.ice_uma_64",
            _LAMMPS_DIR / "data.ice_uma_32",
        ],
        help="LAMMPS data files to diagnose",
    )
    ap.add_argument("--skip-scaling", action="store_true")
    ap.add_argument("--skip-tip4pew", action="store_true")
    ap.add_argument("--no-relax", action="store_true", help="Skip geometry relaxation")
    args = ap.parse_args()

    for dp in args.data:
        if not dp.is_file():
            print(f"Skip missing file: {dp}")
            continue
        diagnose_data_file(dp, relax=not args.no_relax)
        if not args.skip_tip4pew:
            tip4pew_sanity(dp)

    if not args.skip_scaling:
        scaling_test([2, 4, 8, 16, 32, 64, 128])

    return 0


if __name__ == "__main__":
    sys.exit(main())
