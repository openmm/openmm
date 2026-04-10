#!/usr/bin/env python3
"""Diagnose atom ordering: compare LAMMPS read_data positions with init.xyz atom-by-atom.

If positions are identical but energies differ after scatter, the issue is
in how LAMMPS's internal state changes between the two evaluations.
"""

import re
import sys
import numpy as np
from pathlib import Path


def parse_ipi_xyz(path: str):
    with open(path) as f:
        n = int(f.readline().strip())
        comment = f.readline().strip()
        symbols, positions = [], []
        for _ in range(n):
            parts = f.readline().split()
            symbols.append(parts[0])
            positions.append([float(x) for x in parts[1:4]])
    return np.array(positions), symbols


def parse_lammps_data_atoms(path: str):
    """Parse atom positions from LAMMPS data file, in file order."""
    with open(path) as f:
        text = f.read()
    lines = text.splitlines()
    type_to_sym = {1: "O", 2: "H"}
    positions, symbols, tags = [], [], []
    in_atoms = False
    for line in lines:
        if line.strip().startswith("Atoms"):
            in_atoms = True
            continue
        if in_atoms and line.strip() == "":
            if positions:
                break
            continue
        if in_atoms:
            parts = line.split()
            if len(parts) >= 5:
                tag = int(parts[0])
                atype = int(parts[1])
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                tags.append(tag)
                symbols.append(type_to_sym.get(atype, "X"))
                positions.append([x, y, z])
    return np.array(positions), symbols, tags


def main():
    data_path = "ipi/lammps_data_ipi_client.data"
    xyz_path = "ipi/init.xyz"

    data_pos, data_sym, data_tags = parse_lammps_data_atoms(data_path)
    xyz_pos, xyz_sym = parse_ipi_xyz(xyz_path)

    n = len(data_pos)
    print(f"LAMMPS data: {n} atoms")
    print(f"init.xyz:    {len(xyz_pos)} atoms")
    assert len(xyz_pos) == n

    print(f"\n--- First 10 atoms ---")
    print(f"{'Idx':>4} {'Tag':>5} {'Data Sym':>8} {'XYZ Sym':>8} "
          f"{'Data X':>12} {'XYZ X':>12} {'ΔX':>12}")
    for i in range(min(10, n)):
        dx = data_pos[i, 0] - xyz_pos[i, 0]
        print(f"{i:4d} {data_tags[i]:5d} {data_sym[i]:>8} {xyz_sym[i]:>8} "
              f"{data_pos[i, 0]:12.6f} {xyz_pos[i, 0]:12.6f} {dx:12.2e}")

    diffs = data_pos - xyz_pos
    max_diff = np.abs(diffs).max()
    rms_diff = np.sqrt(np.mean(diffs ** 2))
    print(f"\nMax position difference: {max_diff:.2e} Å")
    print(f"RMS position difference: {rms_diff:.2e} Å")

    sym_match = all(d == x for d, x in zip(data_sym, xyz_sym))
    print(f"Symbol order match: {sym_match}")
    if not sym_match:
        mismatches = [(i, data_sym[i], xyz_sym[i])
                      for i in range(n) if data_sym[i] != xyz_sym[i]]
        print(f"  {len(mismatches)} symbol mismatches!")
        for i, ds, xs in mismatches[:10]:
            print(f"    Atom {i}: data={ds}, xyz={xs}")

    tag_order = np.array(data_tags)
    print(f"\nTag order: min={tag_order.min()}, max={tag_order.max()}, "
          f"sorted={np.all(np.diff(tag_order) > 0)}")

    print(f"\nData symbols summary: {data_sym[:6]}...")
    print(f"XYZ  symbols summary: {xyz_sym[:6]}...")

    # Also check if LAMMPS internal order differs from tag order
    from lammps import lammps
    lmp = lammps(cmdargs=["-nocite", "-log", "none", "-echo", "none"])
    lmp.commands_list([
        "units metal",
        "atom_style atomic",
        "boundary p p p",
        "pair_style zero 1.0",
        f"read_data {Path(data_path).resolve()}",
        "pair_coeff * *",
    ])

    x_internal = np.array(lmp.numpy.extract_atom("x")[:n], dtype=np.float64)
    types_internal = np.array(lmp.numpy.extract_atom("type")[:n], dtype=int)
    tags_internal = np.array(lmp.numpy.extract_atom("id")[:n], dtype=int)

    print(f"\n--- LAMMPS internal order ---")
    print(f"Tags[0:10]: {tags_internal[:10].tolist()}")
    print(f"Types[0:10]: {types_internal[:10].tolist()}")
    print(f"Tags sorted: {np.all(np.diff(tags_internal) > 0)}")

    data_vs_internal = np.abs(data_pos - x_internal).max()
    xyz_vs_internal = np.abs(xyz_pos - x_internal).max()
    print(f"\nData pos vs LAMMPS internal pos max diff: {data_vs_internal:.2e} Å")
    print(f"XYZ pos vs LAMMPS internal pos max diff: {xyz_vs_internal:.2e} Å")

    # Check type mapping
    internal_syms = ["O" if t == 1 else "H" for t in types_internal]
    type_match_data = all(d == i for d, i in zip(data_sym, internal_syms))
    type_match_xyz = all(x == i for x, i in zip(xyz_sym, internal_syms))
    print(f"\nData symbols match LAMMPS internal types: {type_match_data}")
    print(f"XYZ symbols match LAMMPS internal types: {type_match_xyz}")

    if not type_match_xyz:
        mismatches = [(i, xyz_sym[i], internal_syms[i])
                      for i in range(n) if xyz_sym[i] != internal_syms[i]]
        print(f"  {len(mismatches)} mismatches between XYZ symbols and LAMMPS types!")
        for idx, xs, ls in mismatches[:10]:
            print(f"    Atom {idx}: xyz={xs}, lammps_type={ls} "
                  f"(tag={tags_internal[idx]}, xyz_pos={xyz_pos[idx]}, lammps_pos={x_internal[idx]})")

    lmp.close()


if __name__ == "__main__":
    main()
