#!/usr/bin/env python3
"""
Prepare ice Ih structure for a cubic periodic box of specified dimension.

Uses GenIce2 (if available) or embedded CIF to generate proton-disordered ice Ih,
then scales the structure to fill exactly an L×L×L nm box at ~0.92 g/cm³ density.

Requirements: ase, (genice2 optional for proton-disordered ice)
  pip install ase
  pip install genice2  # optional, for proton-disordered ice

Usage:
  python prepare_ice_box.py --box 1.0 -o ice.pdb
  python prepare_ice_box.py --box 0.5 --format gro -o ice.gro

Output: PDB or GRO file with ice Ih in cubic box. Use with OpenMM, GROMACS, etc.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

try:
    from ase.io import read, write
except ImportError:
    sys.exit("prepare_ice_box.py requires ASE. Install with: pip install ase")

# Ice Ih density g/cm³ (experimental)
ICE_DENSITY_G_CM3 = 0.92
# Water molar mass g/mol
MOLAR_MASS_WATER = 18.015
# Avogadro
NA = 6.02214076e23

# Embedded ice Ih CIF (12 molecules) - fallback when GenIce2 unavailable
# P 63 c m, a=b=7.82 A, c=7.36 A, gamma=120
_ICE_IH_CIF = """
data_1011023
_cell_length_a 7.82
_cell_length_b 7.82
_cell_length_c 7.36
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
_symmetry_space_group_name_H-M 'P 63 c m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
 O1 O 0.3333 0. 0.0625
 O2 O 0.6667 0. 0.9375
 H1 H 0.3333 0. 0.174
 H2 H 0.438 0. 0.026
 H3 H 0.772 0.105 0.975
"""


def molecules_for_box(box_nm: float, density_g_cm3: float = ICE_DENSITY_G_CM3) -> int:
    """Number of water molecules to fill cubic box L×L×L at given density."""
    vol_cm3 = (box_nm * 1e-7) ** 3  # 1 nm = 1e-7 cm, so 1 nm³ = 1e-21 cm³
    mass_g = density_g_cm3 * vol_cm3
    n_mol = mass_g / MOLAR_MASS_WATER * NA
    return max(1, int(round(n_mol)))


def _get_ice_atoms(num_molecules: int, ice_type: str = "1h") -> "atoms":
    """Get ASE Atoms for ice with at least num_molecules."""
    n_per_cell = 12 if ice_type == "1h" else 8  # ice Ic: 8 molecules per cell (Fd-3m)
    try:
        r = max(1, int(np.ceil((num_molecules / n_per_cell) ** (1 / 3))))
        with tempfile.NamedTemporaryFile(suffix=".y", delete=False) as f:
            outpath = f.name
        result = subprocess.run(
            [
                "genice2",
                ice_type,
                "--rep",
                str(r),
                str(r),
                str(r),
                "-f",
                "y",
                "-o",
                outpath,
                "--seed",
                "42",
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode == 0:
            atoms = read(outpath, format="yaml")
            Path(outpath).unlink(missing_ok=True)
            return atoms
    except Exception:
        pass

    # Fallback: embedded CIF
    from io import StringIO

    atoms = read(StringIO(_ICE_IH_CIF), format="cif")
    n_base = len(atoms) // 3
    if num_molecules <= n_base:
        return atoms
    from ase.build import make_supercell

    r = max(1, int(np.ceil((num_molecules / n_base) ** (1 / 3))))
    P = np.diag([r, r, r])
    return make_supercell(atoms, P)


def _hexagonal_to_orthorhombic(atoms):
    """Transform hexagonal ice Ih cell to orthorhombic. Preserves Cartesian positions."""
    cell = np.array(atoms.get_cell())
    pos = atoms.get_positions()
    cell_ortho = np.array([
        cell[0],
        2.0 * cell[1] - cell[0],
        cell[2],
    ])
    frac_new = pos @ np.linalg.inv(cell_ortho.T)
    frac_wrapped = frac_new % 1.0
    from ase import Atoms

    return Atoms(
        symbols=atoms.get_chemical_symbols(),
        positions=frac_wrapped @ cell_ortho,
        cell=cell_ortho,
        pbc=True,
    )


def prepare_ice_cubic_box(box_nm: float, num_molecules: int | None = None, ice_type: str = "1h") -> "Atoms":
    """
    Prepare ice Ih in a roughly cubic periodic box. Uses orthorhombic transformation
    + uniform scaling to preserve O-O distances (~2.75 Å) and tetrahedral angles.

    Parameters
    ----------
    box_nm : float
        Target box volume = (box_nm)^3 nm³.
    num_molecules : int, optional
        Override molecule count. Default: computed from density.
    ice_type : str
        '1h' (hexagonal) or '1c' (cubic ice, naturally cubic cell).

    Returns
    -------
    atoms : ase.Atoms
        ASE Atoms with ice in orthorhombic/cubic box.
    """
    if num_molecules is None:
        num_molecules = molecules_for_box(box_nm)

    atoms_full = _get_ice_atoms(num_molecules, ice_type=ice_type)
    if ice_type == "1c":
        atoms_ortho = atoms_full.copy()
    else:
        atoms_ortho = _hexagonal_to_orthorhombic(atoms_full)

    from ase.build import make_supercell

    cell = np.array(atoms_ortho.get_cell())
    Lx, Ly, Lz = np.linalg.norm(cell[0]), np.linalg.norm(cell[1]), np.linalg.norm(cell[2])
    n_mol_per_cell = len(atoms_ortho) // 3
    r_min = max(1, int(np.ceil((num_molecules / n_mol_per_cell) ** (1.0 / 3.0))))
    inv_ratios = np.array([1.0 / Lx, 1.0 / Ly, 1.0 / Lz])
    inv_ratios = inv_ratios / inv_ratios.min()
    r_scale = max(r_min, 1)
    na = max(1, int(round(r_scale * inv_ratios[0])))
    nb = max(1, int(round(r_scale * inv_ratios[1])))
    nc = max(1, int(round(r_scale * inv_ratios[2])))
    P = np.diag([na, nb, nc])
    supercell = make_supercell(atoms_ortho, P)

    symbols = supercell.get_chemical_symbols()
    pos_full = supercell.get_positions()
    mol_list = []
    seen = set()
    for i, s in enumerate(symbols):
        if s == "O" and i not in seen:
            oidx = i
            dists = [
                (np.linalg.norm(pos_full[oidx] - pos_full[j]), j)
                for j, sj in enumerate(symbols)
                if sj == "H"
            ]
            dists.sort(key=lambda x: x[0])
            h1, h2 = dists[0][1], dists[1][1]
            mol_list.append([oidx, h1, h2])
            seen.update([oidx, h1, h2])
    mol_list = mol_list[:num_molecules]

    pos_ang = np.array([pos_full[idx] for mol in mol_list for idx in mol])
    cell_super = np.array(supercell.get_cell())
    vol_angstrom3 = np.abs(np.linalg.det(cell_super))

    vol_target = (box_nm * 10.0) ** 3
    scale = (vol_target / vol_angstrom3) ** (1.0 / 3.0)

    pos_scaled = pos_ang * scale
    cell_scaled = cell_super * scale

    from ase import Atoms

    symbols_out = ["O", "H", "H"] * len(mol_list)
    out = Atoms(
        symbols=symbols_out,
        positions=pos_scaled,
        cell=cell_scaled,
        pbc=True,
    )
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Prepare ice Ih for a cubic periodic box"
    )
    ap.add_argument(
        "--box",
        "-b",
        type=float,
        required=True,
        help="Cubic box side length in nm (e.g. 0.5, 1.0, 2.0)",
    )
    ap.add_argument(
        "--output",
        "-o",
        type=str,
        default="ice_box.pdb",
        help="Output file (default: ice_box.pdb)",
    )
    ap.add_argument(
        "--format",
        "-f",
        type=str,
        default="pdb",
        choices=["pdb", "gro", "xyz"],
        help="Output format (default: pdb)",
    )
    ap.add_argument(
        "--molecules",
        "-n",
        type=int,
        default=None,
        help="Override number of molecules (default: from density)",
    )
    ap.add_argument(
        "--ice-type",
        type=str,
        default="1h",
        choices=["1h", "1c"],
        help="Ice type: 1h (hexagonal Ih) or 1c (cubic Ic). Default: 1h",
    )
    args = ap.parse_args()

    atoms = prepare_ice_cubic_box(
        args.box,
        num_molecules=args.molecules,
        ice_type=getattr(args, "ice_type", "1h"),
    )

    ext = {"pdb": ".pdb", "gro": ".gro", "xyz": ".xyz"}[args.format]
    outpath = args.output if args.output.endswith(ext) else args.output + ext
    write(outpath, atoms)  # ASE infers format from extension

    n_mol = len(atoms) // 3
    cell = atoms.get_cell()
    vol_ang3 = np.abs(np.linalg.det(np.array(cell)))
    vol_nm3 = vol_ang3 * 1e-3  # 1 nm³ = 10³ Å³
    density = (n_mol * 18.015 * 1.66054e-24) / (vol_nm3 * 1e-21)
    L_sides = np.linalg.norm(np.array(cell), axis=1) * 0.1  # Å to nm
    print(f"Ice prepared: {L_sides[0]*10:.2f} x {L_sides[1]*10:.2f} x {L_sides[2]*10:.2f} Å box")
    print(f"  Molecules: {n_mol}")
    print(f"  Density: {density:.3f} g/cm³")
    print(f"  Output: {outpath}")


if __name__ == "__main__":
    main()
