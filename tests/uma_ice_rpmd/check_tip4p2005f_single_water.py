#!/usr/bin/env python3
"""
Single-water sanity check for TIP4P/2005f force field.

Builds one water with tip4p2005f XML at equilibrium geometry (r_eq, angle from XML),
computes PE and max |force|, and asserts they are in a sane range. Used to verify
force field units before running full ice simulations.

Usage:
  cd tests/uma_ice_rpmd
  python check_tip4p2005f_single_water.py
  python check_tip4p2005f_single_water.py --force-field path/to/tip4p2005f.xml
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent

if os.getenv("OPENMM_PLUGIN_DIR") is None:
    for _base in [os.getenv("CONDA_PREFIX"), os.path.expanduser("~/miniconda3")]:
        if _base:
            _plugins = os.path.join(_base, "lib", "plugins")
            if os.path.isdir(_plugins):
                os.environ["OPENMM_PLUGIN_DIR"] = _plugins
                break

import openmm
from openmm import Context, Platform, unit
from openmm.app import ForceField, Modeller, PDBFile, NoCutoff


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


def _make_single_water_pdb(req_nm: float, angle_rad: float) -> str:
    """PDB content for one HOH at equilibrium: O at origin, H1 along x, H2 in xy plane."""
    # O at origin; H1 at (req, 0, 0); H2 at req*(cos(angle), sin(angle), 0) so H1-O-H2 = angle
    o = np.array([0.0, 0.0, 0.0])
    h1 = np.array([req_nm, 0.0, 0.0])
    h2 = np.array([req_nm * np.cos(angle_rad), req_nm * np.sin(angle_rad), 0.0])
    # PDB: cols 31-38, 39-46, 47-54 for x,y,z (8.3f each); 77-78 element (O or H)
    o_ang, h1_ang, h2_ang = o * 10, h1 * 10, h2 * 10
    def atom_line(serial, name_4, x, y, z, elem):
        return f"HETATM{serial:5d}  {name_4:<4}HOH     1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {elem:>2}"
    lines = [
        atom_line(1, "O", o_ang[0], o_ang[1], o_ang[2], "O"),
        atom_line(2, "H1", h1_ang[0], h1_ang[1], h1_ang[2], "H"),
        atom_line(3, "H2", h2_ang[0], h2_ang[1], h2_ang[2], "H"),
        "END",
    ]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="Single-water TIP4P/2005f energy sanity check")
    ap.add_argument("--force-field", type=Path, default=None, help="Path to tip4p2005f.xml")
    ap.add_argument("--max-pe-magnitude", type=float, default=1e5,
                    help="Max allowed |PE| per water in kJ/mol (default 1e5)")
    ap.add_argument("--max-force", type=float, default=1e5,
                    help="Max allowed max |force| in kJ/(mol·nm) (default 1e5)")
    args = ap.parse_args()

    ff_path = args.force_field if args.force_field is not None else _tip4p2005f_forcefield_path()
    if not ff_path.is_file():
        print(f"Force field not found: {ff_path}")
        return 1

    # Equilibrium geometry from XML: req=0.09419 nm, angle=1.87448 rad
    req_nm = 0.09419
    angle_rad = 1.87448
    pdb_string = _make_single_water_pdb(req_nm, angle_rad)
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
        f.write(pdb_string)
        pdb_path = f.name
    try:
        pdb = PDBFile(pdb_path)
    finally:
        Path(pdb_path).unlink(missing_ok=True)

    forcefield = ForceField(str(ff_path))
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addExtraParticles(forcefield)
    n_total = modeller.topology.getNumAtoms()
    assert n_total == 4, "Expected 4 atoms (O, H1, H2, M)"

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
        rigidWater=False,
    )
    integrator = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.001 * unit.picoseconds)
    platform = Platform.getPlatformByName("Reference")
    context = Context(system, integrator, platform)
    context.setPositions(modeller.positions)

    state = context.getState(getEnergy=True, getForces=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    forces = state.getForces()
    max_f = 0.0
    for f in forces:
        fx = f[0].value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        fy = f[1].value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        fz = f[2].value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        max_f = max(max_f, (fx*fx + fy*fy + fz*fz)**0.5)

    print(f"Single water (TIP4P/2005f equilibrium): PE = {pe:.4f} kJ/mol  max|F| = {max_f:.4f} kJ/(mol·nm)")

    if abs(pe) > args.max_pe_magnitude:
        print(f"FAIL: |PE| = {abs(pe):.2f} exceeds threshold {args.max_pe_magnitude}. Check force field units.")
        return 1
    if max_f > args.max_force:
        print(f"FAIL: max|force| = {max_f:.2f} exceeds threshold {args.max_force}. Check force field units.")
        return 1
    print("PASS: PE and forces in sane range.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
