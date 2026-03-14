#!/usr/bin/env python3
"""
Convert LAMMPS custom .lammpstrj (id type mass x y z ...) → extended XYZ for Ovito.

Usage:
  python lammpstrj_to_extxyz.py dump.ice_uma_quick.lammpstrj data.ice_uma -o ice_ovito.extxyz
"""
from __future__ import annotations

import argparse
from pathlib import Path


def _read_box_angstrom(data_path: Path) -> tuple[float, float, float, float, float, float]:
    lx = ly = lz = xy = xz = yz = 0.0
    for line in data_path.read_text().splitlines():
        if "xlo xhi" in line:
            a, b = map(float, line.split()[:2])
            lx = b - a
        elif "ylo yhi" in line:
            a, b = map(float, line.split()[:2])
            ly = b - a
        elif "zlo zhi" in line:
            a, b = map(float, line.split()[:2])
            lz = b - a
        elif "xy xz yz" in line:
            p = line.split()
            xy, xz, yz = float(p[0]), float(p[1]), float(p[2])
    return lx, ly, lz, xy, xz, yz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lammpstrj", type=Path)
    ap.add_argument("data", type=Path, help="LAMMPS data file (for lattice)")
    ap.add_argument("-o", "--output", type=Path, default=None)
    args = ap.parse_args()
    lx, ly, lz, xy, xz, yz = _read_box_angstrom(args.data)
    lattice = f'"{lx:.10f} 0 0  {xy:.10f} {ly:.10f} 0  {xz:.10f} {yz:.10f} {lz:.10f}"'
    props = "species:S:1:pos:R:3"
    out = args.output or args.lammpstrj.with_suffix(".extxyz")
    type_to_el = {1: "O", 2: "H"}

    lines = args.lammpstrj.read_text().splitlines()
    blocks: list[tuple[int, list[tuple[str, float, float, float]]]] = []
    i, n = 0, len(lines)
    while i < n:
        if not lines[i].startswith("ITEM: TIMESTEP"):
            i += 1
            continue
        step = int(lines[i + 1])
        i += 2
        if i >= n or not lines[i].startswith("ITEM: NUMBER OF ATOMS"):
            continue
        nat = int(lines[i + 1])
        i += 2
        while i < n and not lines[i].startswith("ITEM: ATOMS"):
            i += 1
        if i >= n:
            break
        header = lines[i].split()[2:]
        itype = header.index("type")
        ix = header.index("x")
        i += 1
        atoms = []
        for _ in range(nat):
            p = lines[i].split()
            t = int(float(p[itype]))
            atoms.append((type_to_el.get(t, "X"), float(p[ix]), float(p[ix + 1]), float(p[ix + 2])))
            i += 1
        blocks.append((step, atoms))

    with out.open("w") as f:
        for step, atoms in blocks:
            f.write(f"{len(atoms)}\n")
            f.write(f'Lattice={lattice} Properties={props} step={step}\n')
            for el, x, y, z in atoms:
                f.write(f"{el}  {x:.8f}  {y:.8f}  {z:.8f}\n")
    print(f"Wrote {out} ({len(blocks)} frames)")


if __name__ == "__main__":
    main()
