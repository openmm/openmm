#!/usr/bin/env python3
"""
Convert HOOMD GSD dimer file (init-0.gsd, Hartree/Bohr/a.u.) to PDB for OpenMM.

Reads a GSD file from cav-hoomd (positions in Bohr, types O/N/L), excludes the
cavity particle L, and writes a PDB with residues OO/NN and atoms A/B so that
OpenMM can load it with diamer_forcefield.xml.

Usage:
  python gsd_to_pdb.py init-0.gsd -o init-0.pdb
  python gsd_to_pdb.py init-0.gsd -o init-0.pdb --frame -1

Conversion:
  - Positions: Bohr to Angstrom (x 0.529177); PDB uses Angstroms.
  - Box: L_bohr to L_Angstrom = L_bohr x 0.529177.
  - Molecular particles only (exclude typeid L); 250 dimers -> 500 atoms.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

try:
    import gsd.hoomd
except ImportError:
    raise SystemExit("gsd_to_pdb.py requires gsd. Install with: pip install gsd")

# 1 Bohr = 0.529177 Angstrom (NIST)
BOHR_TO_ANGSTROM = 0.529177


def _format_atom_line(serial: int, atom_name: str, res_name: str, chain_id: str,
                     res_seq: int, x: float, y: float, z: float, element: str) -> str:
    """Format one ATOM record for PDB (coordinates in Angstroms).
    Element must be in columns 77-78 so OpenMM's reader uses it; otherwise
    atom name 'B' is guessed as Boron and ForceField template matching fails.
    """
    # Build through column 76, then 2-char element at 77-78 (0-indexed: 76:78)
    core = (
        f"ATOM  {serial:5d}  {atom_name:<4s}{res_name:>3s} {chain_id}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
    )
    elem_2 = f"{element:>2s}"  # right-justify in cols 77-78
    line = core.ljust(76) + elem_2 + "\n"
    return line


def gsd_to_pdb(gsd_path: str | Path, pdb_path: str | Path, frame: int = 0) -> None:
    """
    Convert GSD dimer snapshot to PDB for OpenMM.

    Parameters
    ----------
    gsd_path : path
        Input GSD file (e.g. init-0.gsd).
    pdb_path : path
        Output PDB file (e.g. init-0.pdb).
    frame : int
        Frame index (default 0). Use -1 for last frame.
    """
    gsd_path = Path(gsd_path)
    pdb_path = Path(pdb_path)
    if not gsd_path.exists():
        raise FileNotFoundError(f"GSD file not found: {gsd_path}")

    with gsd.hoomd.open(gsd_path, "r") as f:
        snap = f[frame]

    box = snap.configuration.box
    Lx, Ly, Lz = float(box[0]), float(box[1]), float(box[2])
    pos = np.array(snap.particles.position, dtype=np.float64)
    bonds_group = np.array(snap.bonds.group)
    bonds_typeid = np.array(snap.bonds.typeid)

    n_mol = len(bonds_group)
    n_particles_mol = 2 * n_mol
    if snap.particles.N != n_particles_mol and snap.particles.N != n_particles_mol + 1:
        raise ValueError(
            f"GSD has {snap.particles.N} particles, expected {n_particles_mol} or {n_particles_mol}+1 (with cavity)"
        )

    # Box in Angstroms for CRYST1
    Lx_a, Ly_a, Lz_a = Lx * BOHR_TO_ANGSTROM, Ly * BOHR_TO_ANGSTROM, Lz * BOHR_TO_ANGSTROM

    lines = []
    lines.append("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1 \n" % (Lx_a, Ly_a, Lz_a))

    serial = 1
    for b in range(n_mol):
        i, j = int(bonds_group[b][0]), int(bonds_group[b][1])
        tid = int(bonds_typeid[b])
        is_oo = tid == 0
        res_name = "OO" if is_oo else "NN"
        elem = "O" if is_oo else "N"
        res_seq = b + 1
        chain_id = "A"
        # Atom A
        r_i = pos[i] * BOHR_TO_ANGSTROM
        lines.append(_format_atom_line(serial, "A", res_name, chain_id, res_seq, r_i[0], r_i[1], r_i[2], elem))
        serial += 1
        # Atom B
        r_j = pos[j] * BOHR_TO_ANGSTROM
        lines.append(_format_atom_line(serial, "B", res_name, chain_id, res_seq, r_j[0], r_j[1], r_j[2], elem))
        serial += 1

    # CONECT records so OpenMM topology has A-B bonds (required for ForceField template matching)
    for b in range(n_mol):
        a1, a2 = 2 * b + 1, 2 * b + 2
        lines.append(f"CONECT{a1:5d}{a2:5d}\n")
        lines.append(f"CONECT{a2:5d}{a1:5d}\n")

    lines.append("END\n")

    pdb_path.write_text("".join(lines))

    box_nm = (Lx * 0.0529177, Ly * 0.0529177, Lz * 0.0529177)
    print(f"Wrote {pdb_path}")
    print(f"  Molecules: {n_mol} (from bonds)")
    print(f"  Box: {Lx:.2f} x {Ly:.2f} x {Lz:.2f} Bohr -> {box_nm[0]:.4f} x {box_nm[1]:.4f} x {box_nm[2]:.4f} nm")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Convert HOOMD GSD dimer file (Bohr/a.u.) to PDB for OpenMM (diamer_forcefield.xml)."
    )
    ap.add_argument("gsd", type=str, help="Input GSD file (e.g. init-0.gsd)")
    ap.add_argument("-o", "--output", type=str, default=None, help="Output PDB file (default: <gsd_stem>.pdb)")
    ap.add_argument("--frame", type=int, default=0, help="Frame index (default 0; use -1 for last)")
    args = ap.parse_args()

    gsd_path = Path(args.gsd)
    pdb_path = Path(args.output) if args.output else gsd_path.with_suffix(".pdb")
    gsd_to_pdb(gsd_path, pdb_path, frame=args.frame)


if __name__ == "__main__":
    main()
