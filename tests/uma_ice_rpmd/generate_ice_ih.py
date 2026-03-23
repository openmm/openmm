#!/usr/bin/env python3
"""
Generate a periodic ice Ih (proton-ordered) crystal structure as a PDB file.

The orthorhombic unit cell of ice Ih contains 8 water molecules. The supercell
is built by replicating this unit cell nx * ny * nz times, giving a total of
8 * nx * ny * nz molecules. The box dimensions are exact integer multiples of
the crystallographic cell, so the structure tiles seamlessly under periodic
boundary conditions with no gaps or extra space.

Proton ordering satisfies the Bernal-Fowler ice rules (each oxygen donates
exactly 2 and accepts exactly 2 hydrogen bonds), corresponding to the ice XI
polymorph.

Usage:
    python generate_ice_ih.py 5 3 3                    # writes ice_Ih_5x3x3.pdb
    python generate_ice_ih.py 3 2 2 -o my_ice.pdb      # custom output name
"""

import argparse
import sys
import numpy as np
from itertools import product


# ---------------------------------------------------------------------------
# Crystallographic constants
# ---------------------------------------------------------------------------
# Orthorhombic unit cell parameters for ice Ih (Angstroms)
A, B, C = 4.497, 7.789, 7.322
CELL = np.array([A, B, C])

# Internal z-puckering parameter
Z1 = 0.0622

# O-H covalent bond length (Angstroms)
R_OH = 0.9572

# 8 oxygen fractional coordinates in the orthorhombic unit cell, derived from
# the hexagonal P6_3/mmc ice Ih structure via cell transformation.
O_FRAC = np.array([
    [0.0, 1.0 / 3, Z1],
    [0.5, 5.0 / 6, Z1],
    [0.5, 1.0 / 6, 1.0 - Z1],
    [0.0, 2.0 / 3, 1.0 - Z1],
    [0.5, 1.0 / 6, 0.5 + Z1],
    [0.0, 2.0 / 3, 0.5 + Z1],
    [0.0, 1.0 / 3, 0.5 - Z1],
    [0.5, 5.0 / 6, 0.5 - Z1],
])


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def _build_neighbor_list(O_cart: np.ndarray, cell: np.ndarray, cutoff: float = 3.0):
    """Find all O-O neighbors within `cutoff`, including periodic images."""
    n = len(O_cart)
    neighbors: dict[int, list] = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(n):
            for sx, sy, sz in product([-1, 0, 1], repeat=3):
                shift = np.array([sx, sy, sz], dtype=float)
                rj = O_cart[j] + shift * cell
                dr = rj - O_cart[i]
                d = np.linalg.norm(dr)
                if 0.1 < d < cutoff:
                    neighbors[i].append((j, shift.copy(), dr.copy(), d))
    return neighbors


def _enumerate_bonds(neighbors: dict):
    """Return a list of unique undirected O-O bonds."""
    bonds = []
    seen: set[tuple] = set()
    for i, nbrs in neighbors.items():
        for j, shift, dr, d in nbrs:
            k1 = (i, j, tuple(shift.astype(int)))
            k2 = (j, i, tuple((-shift).astype(int)))
            if k1 not in seen and k2 not in seen:
                seen.add(k1)
                bonds.append((i, j, shift.copy(), dr.copy()))
    return bonds


def _solve_proton_ordering(bonds: list, n_atoms: int = 8) -> list[int]:
    """
    Backtracking solver for ice-rule proton ordering.

    Returns a list of orientations (0 or 1) for each bond.  orientation[k] == 0
    means atom i is the donor for bond k; orientation[k] == 1 means atom j is
    the donor.
    """
    n_bonds = len(bonds)
    orientation = [0] * n_bonds

    def _count(atom: int, depth: int) -> int:
        c = 0
        for k in range(depth):
            bi, bj, _, _ = bonds[k]
            if orientation[k] == 0 and bi == atom:
                c += 1
            elif orientation[k] == 1 and bj == atom:
                c += 1
        return c

    def _bt(depth: int) -> bool:
        if depth == n_bonds:
            return all(_count(at, depth) == 2 for at in range(n_atoms))
        bi, bj, _, _ = bonds[depth]
        for orient in (0, 1):
            orientation[depth] = orient
            donor = bi if orient == 0 else bj
            if _count(donor, depth + 1) > 2:
                continue
            ok = True
            for at in range(n_atoms):
                cur = _count(at, depth + 1)
                rem = sum(
                    1
                    for k in range(depth + 1, n_bonds)
                    if bonds[k][0] == at or bonds[k][1] == at
                )
                if cur > 2 or cur + rem < 2:
                    ok = False
                    break
            if ok and _bt(depth + 1):
                return True
        return False

    if not _bt(0):
        raise RuntimeError("Failed to find a valid proton ordering")
    return orientation


def _place_hydrogens(O_cart, bonds, orientation):
    """Place two H atoms per O along its donor bond directions."""
    donor_vecs: dict[int, list] = {i: [] for i in range(len(O_cart))}
    for k, (i, j, _shift, dr) in enumerate(bonds):
        if orientation[k] == 0:
            donor_vecs[i].append(dr.copy())
        else:
            donor_vecs[j].append(-dr.copy())

    waters = []
    for i in range(len(O_cart)):
        O = O_cart[i]
        d1 = donor_vecs[i][0]
        d2 = donor_vecs[i][1]
        H1 = O + R_OH * d1 / np.linalg.norm(d1)
        H2 = O + R_OH * d2 / np.linalg.norm(d2)
        waters.append((O.copy(), H1.copy(), H2.copy()))
    return waters


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def build_unit_cell():
    """Build the 8-molecule orthorhombic ice Ih unit cell (O + 2H each)."""
    O_cart = O_FRAC * CELL
    neighbors = _build_neighbor_list(O_cart, CELL)
    bonds = _enumerate_bonds(neighbors)
    orientation = _solve_proton_ordering(bonds)
    return _place_hydrogens(O_cart, bonds, orientation)


def replicate(unit_cell_waters, nx: int, ny: int, nz: int):
    """
    Replicate the unit cell and wrap atoms into the supercell box.

    Returns (waters, supercell) where waters is a list of (O, H1, H2) tuples
    and supercell is a 3-element array of box dimensions.
    """
    supercell = CELL * np.array([nx, ny, nz])
    waters = []
    for ix, iy, iz in product(range(nx), range(ny), range(nz)):
        shift = np.array([ix, iy, iz], dtype=float) * CELL
        for O, H1, H2 in unit_cell_waters:
            O_new = O + shift
            H1_new = (H1 + shift).copy()
            H2_new = (H2 + shift).copy()
            # Wrap H atoms that sit outside the supercell (bonds crossing
            # the unit-cell boundary in the original cell).
            for H in (H1_new, H2_new):
                for dim in range(3):
                    while H[dim] < 0:
                        H[dim] += supercell[dim]
                    while H[dim] >= supercell[dim]:
                        H[dim] -= supercell[dim]
            waters.append((O_new, H1_new, H2_new))
    return waters, supercell


def write_pdb(waters, supercell, nx, ny, nz, path: str):
    """Write the ice crystal to a PDB file with CRYST1 record."""
    n_mol = len(waters)
    vol = supercell[0] * supercell[1] * supercell[2]
    density = (n_mol * 18.015 / 6.022e23) / (vol * 1e-24)

    with open(path, "w") as f:
        f.write(
            f"CRYST1{supercell[0]:9.3f}{supercell[1]:9.3f}{supercell[2]:9.3f}"
            f"  90.00  90.00  90.00 P 1           1\n"
        )
        f.write(f"REMARK   Ice Ih crystal, orthorhombic supercell {nx}x{ny}x{nz}\n")
        f.write(f"REMARK   {n_mol} water molecules, proton-ordered (ice XI)\n")
        f.write(
            f"REMARK   Unit cell: a={A:.3f} b={B:.3f} c={C:.3f} Angstroms\n"
        )
        f.write(f"REMARK   Density: {density:.4f} g/cm3\n")

        serial = 0
        for mol_idx, (O, H1, H2) in enumerate(waters):
            res = (mol_idx % 9999) + 1
            for name, pos, elem in [
                ("OW ", O, "O"),
                ("HW1", H1, "H"),
                ("HW2", H2, "H"),
            ]:
                serial += 1
                f.write(
                    f"ATOM  {serial:5d}  {name} HOH {res:4d}    "
                    f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                    f"  1.00  0.00           {elem}\n"
                )
        f.write("END\n")

    return density


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate a periodic ice Ih crystal PDB file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "The orthorhombic unit cell contains 8 water molecules.\n"
            "Total molecules = 8 * nx * ny * nz.\n\n"
            "Examples:\n"
            "  %(prog)s 3 2 2                  # 96 molecules,  ~13.5 x 15.6 x 14.6 A\n"
            "  %(prog)s 5 3 3                  # 360 molecules, ~22.5 x 23.4 x 22.0 A\n"
            "  %(prog)s 7 4 4 -o big_ice.pdb   # 896 molecules, ~31.5 x 31.2 x 29.3 A\n"
        ),
    )
    parser.add_argument("nx", type=int, help="Replications along a (4.497 A)")
    parser.add_argument("ny", type=int, help="Replications along b (7.789 A)")
    parser.add_argument("nz", type=int, help="Replications along c (7.322 A)")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output PDB path (default: ice_Ih_NXxNYxNZ.pdb)",
    )
    args = parser.parse_args()

    if args.nx < 1 or args.ny < 1 or args.nz < 1:
        parser.error("All replication factors must be >= 1")

    out = args.output or f"ice_Ih_{args.nx}x{args.ny}x{args.nz}.pdb"
    n_mol = 8 * args.nx * args.ny * args.nz

    print(f"Building ice Ih crystal: {args.nx}x{args.ny}x{args.nz} supercell")
    print(f"  Molecules: {n_mol}")
    print(f"  Atoms:     {n_mol * 3}")

    unit_cell = build_unit_cell()
    waters, supercell = replicate(unit_cell, args.nx, args.ny, args.nz)
    density = write_pdb(waters, supercell, args.nx, args.ny, args.nz, out)

    print(f"  Box:       {supercell[0]:.3f} x {supercell[1]:.3f} x {supercell[2]:.3f} A")
    print(f"  Density:   {density:.4f} g/cm3")
    print(f"  Output:    {out}")


if __name__ == "__main__":
    main()
