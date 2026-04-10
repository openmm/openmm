#!/usr/bin/env python3
"""
Wrap PDB trajectory positions into the primary simulation cell (minimum image).

Reads a multi-frame PDB with CRYST1 (orthorhombic box), wraps all atom positions
into [0, L) per axis, and writes a new PDB. Preserves residue/chain info.

Usage:
  python wrap_pdb.py ice_rpmd_uma_T243_b8.pdb -o ice_rpmd_wrapped.pdb
  python wrap_pdb.py input.pdb                  # writes input_wrapped.pdb

Orthorhombic PBC wrapping:
  frac = pos / L
  pos_wrapped = (frac - floor(frac)) * L   =>  pos in [0, L)
"""
import argparse
import sys
from pathlib import Path


def parse_cryst1(line: str) -> tuple[float, float, float]:
    """Parse CRYST1 line; return (a, b, c) in Angstrom."""
    # CRYST1   15.646   14.707   18.073  90.00  90.00  90.00 P 1           1
    parts = line.split()
    if len(parts) < 4:
        raise ValueError(f"Cannot parse CRYST1: {line}")
    return float(parts[1]), float(parts[2]), float(parts[3])


def wrap_pos_ortho(x: float, y: float, z: float, lx: float, ly: float, lz: float) -> tuple[float, float, float]:
    """Wrap Cartesian position into [0, L) for orthorhombic box."""
    import math
    xw = x - lx * math.floor(x / lx)
    yw = y - ly * math.floor(y / ly)
    zw = z - lz * math.floor(z / lz)
    # Avoid 1.0 edge due to float; clamp to just under L
    if xw >= lx:
        xw = 0.0
    if yw >= ly:
        yw = 0.0
    if zw >= lz:
        zw = 0.0
    return xw, yw, zw


# PDB fixed-column format (standard)
# Col 1-6: record, 7-11: serial, 13-16: name, 17-20: resname, 22: chain, 23-26: resseq
# Col 31-38: x, 39-46: y, 47-54: z, 55-60: occ, 61-66: b, 77-78: element


def parse_atom_line(line: str) -> tuple | None:
    """Parse ATOM/HETATM line using PDB column convention. Returns (prefix, serial, name, resname, chain, resseq, x, y, z, occ, b, rest) or None."""
    if len(line) < 54:
        return None
    try:
        prefix = line[:6]
        serial = int(line[6:11])
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21:22] if len(line) > 21 else " "
        resseq = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        occ = float(line[54:60]) if len(line) > 60 else 1.0
        b = float(line[60:66]) if len(line) > 66 else 0.0
        rest = line[66:].rstrip("\n") if len(line) > 66 else ""
        return (prefix, serial, name, resname, chain, resseq, x, y, z, occ, b, rest)
    except (ValueError, IndexError):
        return None


def format_atom_line(prefix: str, serial: int, name: str, resname: str, chain: str, resseq: int,
                    x: float, y: float, z: float, occ: float, b: float, rest: str) -> str:
    """Format ATOM line (PDB standard columns)."""
    return f"{prefix}{serial:5d} {name:<4s}{resname:<3s} {chain}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}          {rest}\n"


def process_pdb(in_path: Path, out_path: Path) -> None:
    """Read PDB, wrap positions, write wrapped PDB."""
    lx, ly, lz = None, None, None
    header_lines = []
    frame_count = 0

    with open(in_path) as f_in, open(out_path, "w") as f_out:
        in_model = False
        for line in f_in:
            if line.startswith("CRYST1"):
                lx, ly, lz = parse_cryst1(line)
                header_lines.append(line)
                continue

            if line.startswith("MODEL"):
                if frame_count == 0 and header_lines:
                    for h in header_lines:
                        f_out.write(h)
                in_model = True
                frame_count += 1
                f_out.write(line)
                continue

            if line.startswith("ENDMDL"):
                in_model = False
                f_out.write(line)
                continue

            if line.startswith("END"):
                f_out.write(line)
                break

            if in_model and (line.startswith("ATOM  ") or line.startswith("HETATM")):
                parsed = parse_atom_line(line)
                if parsed is None:
                    f_out.write(line)
                    continue
                prefix, serial, name, resname, chain, resseq, x, y, z, occ, b, rest = parsed
                xw, yw, zw = wrap_pos_ortho(x, y, z, lx, ly, lz)
                f_out.write(format_atom_line(
                    prefix, serial, name, resname, chain, resseq,
                    xw, yw, zw, occ, b, rest
                ))
                continue

            # REMARK, etc. - keep before first MODEL
            if not in_model and frame_count == 0:
                header_lines.append(line)
            elif in_model:
                f_out.write(line)

    print(f"Wrapped {frame_count} frames from {in_path} -> {out_path}")
    print(f"Box: {lx:.2f} x {ly:.2f} x {lz:.2f} Å")


def main():
    parser = argparse.ArgumentParser(description="Wrap PDB trajectory into primary cell")
    parser.add_argument("input", type=Path, help="Input PDB file")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="Output PDB (default: input_wrapped.pdb)")
    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: {args.input} not found", file=sys.stderr)
        sys.exit(1)

    out = args.output
    if out is None:
        stem = args.input.stem
        out = args.input.parent / f"{stem}_wrapped.pdb"

    process_pdb(args.input, out)


if __name__ == "__main__":
    main()
