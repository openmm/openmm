#!/usr/bin/env python3
"""Standalone UMA energy evaluation on init.xyz and init_openmm_rpmd_64.xyz.

Bypasses i-PI/LAMMPS entirely to isolate whether the 323 eV PE gap
is from the structure itself or the evaluation pathway.

Expected results:
  - LAMMPS minimized energy: -133,146 eV
  - i-PI reported step-0 PE: -132,823 eV  (323 eV higher — the mystery)
"""

import sys
import re
import numpy as np
import torch
from ase import Atoms
from ase.io import read
from ase.geometry import wrap_positions
from fairchem.core import pretrained_mlip
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator


def parse_ipi_xyz(path: str) -> Atoms:
    """Parse i-PI extended XYZ with CELL(abcABC) comment."""
    with open(path) as f:
        n = int(f.readline().strip())
        comment = f.readline().strip()
        symbols, positions = [], []
        for _ in range(n):
            parts = f.readline().split()
            symbols.append(parts[0])
            positions.append([float(x) for x in parts[1:4]])

    cell_match = re.search(
        r"CELL\(abcABC\):\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
        comment,
    )
    if cell_match:
        a, b, c = float(cell_match.group(1)), float(cell_match.group(2)), float(cell_match.group(3))
        cell = [a, b, c]
    else:
        raise ValueError(f"Could not parse cell from: {comment}")

    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    return atoms


def evaluate_structure(atoms: Atoms, label: str, calc) -> dict:
    """Evaluate UMA energy on an Atoms object."""
    print(f"\n--- {label} ---")
    print(f"  N atoms: {len(atoms)}")
    print(f"  Cell:    {atoms.cell.lengths()}")
    print(f"  PBC:     {atoms.pbc}")
    pos = atoms.get_positions()
    print(f"  Pos X:   [{pos[:,0].min():.4f}, {pos[:,0].max():.4f}]")
    print(f"  Pos Y:   [{pos[:,1].min():.4f}, {pos[:,1].max():.4f}]")
    print(f"  Pos Z:   [{pos[:,2].min():.4f}, {pos[:,2].max():.4f}]")

    atoms.info["charge"] = 0
    atoms.info["spin"] = 1
    atoms.calc = calc
    energy_eV = float(atoms.get_potential_energy())
    forces = np.asarray(atoms.get_forces(), dtype=np.float64)
    fmax = np.abs(forces).max()
    fnorm = np.linalg.norm(forces)

    print(f"  UMA energy:     {energy_eV:.4f} eV")
    print(f"  UMA energy:     {energy_eV * 96.4853:.2f} kJ/mol")
    print(f"  Force max comp: {fmax:.6f} eV/Å")
    print(f"  Force two-norm: {fnorm:.6f} eV/Å")
    return {"energy_eV": energy_eV, "fmax": fmax, "fnorm": fnorm}


def main():
    device = sys.argv[1] if len(sys.argv) > 1 else "cuda"

    predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device=device)
    calc = FAIRChemCalculator(predictor, task_name="omol")

    # 1) i-PI init.xyz (coordinates in [-L/2, L/2])
    ipi_atoms = parse_ipi_xyz("ipi/init.xyz")
    r_ipi_raw = evaluate_structure(ipi_atoms.copy(), "i-PI init.xyz (raw, [-L/2,L/2])", calc)

    # 2) i-PI init.xyz with wrapped coords (to [0, L))
    ipi_wrapped = ipi_atoms.copy()
    ipi_wrapped.positions = wrap_positions(
        ipi_wrapped.positions, ipi_wrapped.cell, pbc=[True, True, True]
    )
    r_ipi_wrap = evaluate_structure(ipi_wrapped, "i-PI init.xyz (wrapped to [0,L))", calc)

    # 3) OpenMM init (coordinates in [0, L))
    openmm_atoms = read("pipeline_out/init_openmm_rpmd_64.xyz")
    r_openmm = evaluate_structure(openmm_atoms.copy(), "OpenMM init (raw)", calc)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  i-PI raw energy:     {r_ipi_raw['energy_eV']:.4f} eV")
    print(f"  i-PI wrapped energy: {r_ipi_wrap['energy_eV']:.4f} eV")
    print(f"  OpenMM energy:       {r_openmm['energy_eV']:.4f} eV")
    print()
    print(f"  Gap (i-PI raw - OpenMM):     {r_ipi_raw['energy_eV'] - r_openmm['energy_eV']:.4f} eV")
    print(f"  Gap (i-PI wrapped - OpenMM): {r_ipi_wrap['energy_eV'] - r_openmm['energy_eV']:.4f} eV")
    print(f"  Gap (raw - wrapped):         {r_ipi_raw['energy_eV'] - r_ipi_wrap['energy_eV']:.4f} eV")
    print()
    print("Reference values:")
    print("  LAMMPS minimize final:  -133146.25 eV")
    print("  i-PI step-0 PE:        -132823 eV (from .md file)")
    print("  OpenMM step-0 PE:      -133146 eV (from CSV)")


if __name__ == "__main__":
    main()
