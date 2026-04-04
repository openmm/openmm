#!/usr/bin/env python3
"""
Convert LAMMPS data file to i-PI xyz format for initialization.

i-PI xyz format:
  Line 1: natoms
  Line 2: # CELL(abcABC): a b c alpha beta gamma
  Lines 3+: symbol x y z (one per atom)

Uses Angstrom for positions and cell (i-PI will convert via units).
Atom order must match LAMMPS data (O,H,H per molecule) for fix ipi parity.

By default, **orthorhombic** cells are shifted so atomic coordinates live in the
symmetric box ``[-Lx/2, Lx/2] × …`` (same convention as LAMMPS ``fix ipi`` after
it re-centers the cell). This keeps i-PI's initial Cartesian positions aligned
with the client's ``read_data`` box and avoids unnecessary PBC wraps on step 0.
Use ``--no-center-for-fix-ipi`` for legacy ``[0, L)`` coordinates.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent
_ROOT = _SCRIPT_DIR.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def shift_orthorhombic_positions_to_symmetric_cell(
    positions: np.ndarray,
    *,
    xlo: float,
    xhi: float,
    ylo: float,
    yhi: float,
    zlo: float,
    zhi: float,
    xy: float,
    xz: float,
    yz: float,
) -> np.ndarray:
    """Subtract cell midpoints so coords sit in ``[-L/2, L/2]`` (orthorhombic only)."""
    if abs(xy) + abs(xz) + abs(yz) > 1e-12:
        return positions
    cx = 0.5 * (xlo + xhi)
    cy = 0.5 * (ylo + yhi)
    cz = 0.5 * (zlo + zhi)
    return positions - np.array([cx, cy, cz], dtype=np.float64)


def _lammps_cell_to_abc_deg(lx: float, ly: float, lz: float, xy: float, xz: float, yz: float):
    """Convert LAMMPS restricted triclinic to a,b,c,alpha,beta,gamma (degrees)."""
    a_vec = np.array([lx, 0.0, 0.0])
    b_vec = np.array([xy, ly, 0.0])
    c_vec = np.array([xz, yz, lz])
    a = np.linalg.norm(a_vec)
    b = np.linalg.norm(b_vec)
    c = np.linalg.norm(c_vec)
    alpha = np.degrees(np.arccos(np.clip(np.dot(b_vec, c_vec) / (b * c), -1, 1)))
    beta = np.degrees(np.arccos(np.clip(np.dot(a_vec, c_vec) / (a * c), -1, 1)))
    gamma = np.degrees(np.arccos(np.clip(np.dot(a_vec, b_vec) / (a * b), -1, 1)))
    return a, b, c, alpha, beta, gamma


def main():
    p = argparse.ArgumentParser(description="LAMMPS data -> i-PI xyz")
    p.add_argument("--data", type=Path, default=_ROOT / "lammps/data.ice_uma")
    p.add_argument("-o", "--output", type=Path, default=_SCRIPT_DIR / "init.xyz")
    p.add_argument(
        "--no-center-for-fix-ipi",
        action="store_true",
        help="Keep LAMMPS Cartesian coords in [xlo,xhi) (legacy); default centers to [-L/2,L/2].",
    )
    args = p.parse_args()

    text = args.data.read_text()
    lines = text.splitlines()
    # Parse LAMMPS data
    i = 0
    while i < len(lines):
        line = lines[i]
        i += 1
        if "atoms" in line and line.split()[0].isdigit():
            natoms = int(line.split()[0])
            break
    # Skip to xlo xhi etc
    while i < len(lines):
        line = lines[i]
        i += 1
        parts = line.split()
        if len(parts) >= 4 and "xlo" in line and "xhi" in line:
            xlo, xhi = float(parts[0]), float(parts[1])
            break
    ly_line = lines[i]
    i += 1
    lz_line = lines[i]
    i += 1
    xy_line = lines[i]
    i += 1
    ylo, yhi = float(ly_line.split()[0]), float(ly_line.split()[1])
    zlo, zhi = float(lz_line.split()[0]), float(lz_line.split()[1])
    xy, xz, yz = 0.0, 0.0, 0.0
    if "xy" in xy_line and "xz" in xy_line and "yz" in xy_line:
        xy, xz, yz = [float(x) for x in xy_line.split()[:3]]
    lx, ly, lz = xhi - xlo, yhi - ylo, zhi - zlo

    # Skip to Masses then Atoms
    while i < len(lines):
        if "Atoms" in lines[i]:
            i += 1
            break
        i += 1

    type_to_symbol = {1: "O", 2: "H"}
    positions = []
    symbols = []
    while len(positions) < natoms and i < len(lines):
        parts = lines[i].split()
        i += 1
        if len(parts) >= 5:
            atype = int(parts[1])
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            positions.append([x, y, z])
            symbols.append(type_to_symbol.get(atype, "X"))

    pos_arr = np.asarray(positions, dtype=np.float64)
    if not args.no_center_for_fix_ipi:
        pos_arr = shift_orthorhombic_positions_to_symmetric_cell(
            pos_arr,
            xlo=xlo,
            xhi=xhi,
            ylo=ylo,
            yhi=yhi,
            zlo=zlo,
            zhi=zhi,
            xy=xy,
            xz=xz,
            yz=yz,
        )

    # i-PI xyz format: comment line must contain CELL(abcABC): a b c alpha beta gamma
    a, b, c, alpha, beta, gamma = _lammps_cell_to_abc_deg(lx, ly, lz, xy, xz, yz)
    out_lines = [
        str(natoms),
        f"# CELL(abcABC): {a:.6f} {b:.6f} {c:.6f} {alpha:.6f} {beta:.6f} {gamma:.6f}  cell{{angstrom}}",
    ]
    for sym, (px, py, pz) in zip(symbols, pos_arr):
        out_lines.append(f"{sym:>4} {px:16.8e} {py:16.8e} {pz:16.8e}")

    args.output.write_text("\n".join(out_lines) + "\n")
    print(f"Wrote {args.output} ({natoms} atoms)")


if __name__ == "__main__":
    main()
