#!/usr/bin/env python3
"""Convert PDB trajectory to XYZ format for Ovito visualization.

Uses extended XYZ format with explicit Properties= and Lattice= so Ovito
can unambiguously map columns (species, positions) and load the simulation cell.

Usage:
  python pdb_to_xyz.py ice_rpmd_uma_T243_b8.pdb -o ice_rpmd.xyz
  python pdb_to_xyz.py trajectory.pdb                  # writes trajectory.xyz
"""
from pathlib import Path
import argparse
from ase.io import read, write


def main():
    parser = argparse.ArgumentParser(description="Convert PDB trajectory to XYZ for Ovito")
    parser.add_argument("input", type=Path, help="Input PDB file")
    parser.add_argument("-o", "--output", type=Path, default=None, help="Output XYZ file")
    args = parser.parse_args()

    out = args.output or args.input.with_suffix(".xyz")
    if out.suffix.lower() != ".xyz":
        out = out.with_suffix(".xyz")

    atoms_list = read(args.input, index=":")
    if not atoms_list:
        raise SystemExit(f"No frames found in {args.input}")

    # Extended XYZ: explicit Properties=species:S:1:pos:R:3 and Lattice= for Ovito
    write(out, atoms_list, format="extxyz")
    print(f"Wrote {len(atoms_list)} frames to {out} (extended XYZ)")


if __name__ == "__main__":
    main()
