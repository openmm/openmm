#!/usr/bin/env python3
"""
Export OpenMM-UMA minimized ice structure as a LAMMPS data file for use as shared
initial structure in the pipeline (--shared-initial-structure).

Both OpenMM and LAMMPS can then start MD from this same minimized configuration
with --skip-minimize to isolate integrator/force differences from minimization differences.

Usage:
  cd tests/uma_ice_rpmd
  python export_minimized_ice.py --data lammps/data.ice_uma -o lammps/data.ice_uma.min
  python run_ice_pipeline.py --shared-initial-structure lammps/data.ice_uma.min ...
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

from openmm import Context, LocalEnergyMinimizer, Platform, Vec3, VerletIntegrator, unit
from openmm.app import Topology, element
from openmmml import MLPotential

sys.path.insert(0, str(_SCRIPT_DIR))
from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data
from run_openmm_ice_lammps_match import water_topology


def write_lammps_data(path: Path, header_lines: list[str], pos_ang: np.ndarray) -> None:
    """Write LAMMPS data file with same header as original but new positions.
    header_lines must include up to and including 'Atoms' and the blank line after it.
    Atoms section format: id type x y z (no mol id).
    """
    lines = list(header_lines)
    n = pos_ang.shape[0]
    for i in range(n):
        # type: 1 = O, 2 = H (same order as data.ice_uma: O,H,H per molecule)
        typ = 1 if (i % 3 == 0) else 2
        lines.append(f"{i + 1} {typ} {pos_ang[i, 0]:.8f} {pos_ang[i, 1]:.8f} {pos_ang[i, 2]:.8f}")
    path.write_text("\n".join(lines) + "\n")


def extract_header_and_atoms_start(data_path: Path) -> tuple[list[str], int]:
    """Return (lines before and including 'Atoms' line, 1-based line index of first atom row)."""
    lines = data_path.read_text().splitlines()
    header: list[str] = []
    first_atom_line = 0
    for i, line in enumerate(lines):
        if line.strip() == "Atoms":
            header = lines[: i + 2]  # include "Atoms" and blank
            first_atom_line = i + 2
            break
    return header, first_atom_line


def main() -> None:
    ap = argparse.ArgumentParser(description="OpenMM UMA minimize ice, export as LAMMPS data file")
    ap.add_argument("--data", type=Path, default=_LAMMPS_DIR / "data.ice_uma")
    ap.add_argument("-o", "--output", type=Path, required=True, help="Output LAMMPS data file (e.g. data.ice_uma.min)")
    ap.add_argument("--model", default="uma-s-1p1-pythonforce-batch")
    ap.add_argument("--minimize-iters", type=int, default=150)
    ap.add_argument("--minimize-tol", type=float, default=10.0)
    ap.add_argument("--platform", default="cuda", choices=["cuda", "cpu"])
    ap.add_argument("--ml-device", default=None)
    args = ap.parse_args()

    if not args.data.is_file():
        print(f"Missing {args.data}", file=sys.stderr)
        sys.exit(1)

    pos_ang, cell, Z = parse_lammps_data(args.data)
    n = Z.shape[0]
    assert n % 3 == 0
    n_mol = n // 3

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

    ml_dev = args.ml_device or ("cuda" if args.platform == "cuda" else None)
    potential = MLPotential(args.model)
    system = potential.createSystem(topology, task_name="omol", charge=0, spin=1, device=ml_dev)
    integrator = VerletIntegrator(1.0 * unit.femtoseconds)

    platform = Platform.getPlatformByName(args.platform.upper())
    props = {"Precision": "mixed", "DeviceIndex": "0"} if args.platform == "cuda" else {}
    context = Context(system, integrator, platform, props)
    context.setPeriodicBoxVectors(a, b, c)
    context.setPositions(pos_list)

    print("Minimizing with OpenMM UMA...")
    t0 = time.perf_counter()
    LocalEnergyMinimizer.minimize(
        context, maxIterations=args.minimize_iters, tolerance=args.minimize_tol
    )
    st = context.getState(getPositions=True, getEnergy=True)
    pos_nm = np.array([[p[0].value_in_unit(unit.nanometer), p[1].value_in_unit(unit.nanometer), p[2].value_in_unit(unit.nanometer)] for p in st.getPositions()])
    pos_ang_out = pos_nm * 10.0
    print(f"  PE after minimize: {st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.2f} kJ/mol")
    print(f"  Wall time: {time.perf_counter() - t0:.1f} s")

    header, _ = extract_header_and_atoms_start(args.data)
    write_lammps_data(args.output, header, pos_ang_out)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
