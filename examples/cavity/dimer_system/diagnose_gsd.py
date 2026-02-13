#!/usr/bin/env python3
"""
Diagnose GSD content and units for debugging OpenMM simulation instability.

Phase 2: Sanity-check GSD content and units.
Phase 4: Verify bond order and position indexing.

Usage:
  python diagnose_gsd.py [gsd_path] [frame]
  Default: init-0.gsd, frame 0
"""

import argparse
import sys
from pathlib import Path

import numpy as np

BOHR_TO_NM = 0.0529177  # 1 Bohr = 0.0529177 nm
# Expected bond lengths: O-O ~2.28 Bohr (~0.12 nm), N-N ~2.07 Bohr (~0.11 nm)
R0_OO_BOHR = 2.281655158
R0_NN_BOHR = 2.0743522177


def min_image_disp(dr, box):
    """Minimum-image displacement; box is (Lx, Ly, Lz)."""
    L = np.asarray(box[:3], dtype=np.float64)
    return dr - np.round(dr / L) * L


def bond_length(pos, i, j, box):
    """Bond length between particles i and j with minimum-image convention."""
    dr = pos[j] - pos[i]
    d = min_image_disp(dr, box)
    return np.sqrt(np.dot(d, d))


def main():
    parser = argparse.ArgumentParser(description="Diagnose GSD file: box, positions, bonds, units.")
    parser.add_argument("gsd", nargs="?", default="init-0.gsd", help="GSD file path (default: init-0.gsd)")
    parser.add_argument("frame", nargs="?", type=int, default=0, help="Frame index (default: 0)")
    parser.add_argument("--overlap-threshold", type=float, default=0.05,
                        help="Report pairs closer than this (nm) as overlap (default: 0.05)")
    args = parser.parse_args()

    try:
        import gsd.hoomd
    except ImportError:
        print("Error: gsd is required. pip install gsd", file=sys.stderr)
        sys.exit(1)

    path = Path(args.gsd)
    if not path.exists():
        print(f"Error: file not found: {path}", file=sys.stderr)
        sys.exit(1)

    with gsd.hoomd.open(path, "r") as f:
        nframes = len(f)
        if args.frame < 0:
            frame_idx = nframes + args.frame
        else:
            frame_idx = args.frame
        if frame_idx < 0 or frame_idx >= nframes:
            print(f"Error: frame {args.frame} out of range (file has {nframes} frames)", file=sys.stderr)
            sys.exit(1)
        snap = f[frame_idx]

    box = np.array(snap.configuration.box[:3], dtype=np.float64)
    pos = np.array(snap.particles.position, dtype=np.float64)
    N = snap.particles.N
    bonds_group = np.array(snap.bonds.group)
    bonds_typeid = np.array(snap.bonds.typeid) if hasattr(snap.bonds, "typeid") else None
    n_mol = len(bonds_group)

    print("=" * 60)
    print("GSD diagnostic:", path.name, f"frame {frame_idx}")
    print("=" * 60)
    print(f"Particles: {N}")
    print(f"Bonds: {n_mol}")
    print()
    print("--- Box (raw) ---")
    print(f"  configuration.box[:3] = {box}")
    print(f"  L = {box[0]:.6g} (first length)")
    if 30 < box[0] < 50:
        print("  Interpretation: ~40 → consistent with BOHR (40 Bohr ≈ 2.12 nm).")
        print(f"  In nm: L_nm = {box[0] * BOHR_TO_NM:.4f} nm")
    elif 1.5 < box[0] < 5:
        print("  Interpretation: ~2.1 → file may be in NM. If we apply Bohr→nm we would get tiny box!")
        print("  Try --gsd-in-nm when loading if simulation blows up.")
    else:
        print("  Interpretation: unclear; check units (Bohr vs nm).")
    print()
    print("--- Positions (raw) ---")
    for ax, name in enumerate(["x", "y", "z"]):
        lo, hi = pos[:, ax].min(), pos[:, ax].max()
        print(f"  {name}: min = {lo:.6g}, max = {hi:.6g}")
    print()
    print("--- Bond order (Phase 4: first 5 bonds.group) ---")
    for b in range(min(5, n_mol)):
        i, j = int(bonds_group[b][0]), int(bonds_group[b][1])
        typ = int(bonds_typeid[b]) if bonds_typeid is not None else None
        typ_str = "O-O" if typ == 0 else "N-N" if typ == 1 else f"typeid={typ}"
        print(f"  bond {b}: group = [{i}, {j}]  ({typ_str})")
    if n_mol >= 5:
        print("  ...")
    print("  Expected for initlattice: [0,1], [2,3], [4,5], ...")
    print()
    print("--- Bond lengths (minimum-image) ---")
    n_show = min(10, n_mol)
    for b in range(n_show):
        i, j = int(bonds_group[b][0]), int(bonds_group[b][1])
        L_raw = bond_length(pos, i, j, box)
        L_nm = L_raw * BOHR_TO_NM
        typ = int(bonds_typeid[b]) if bonds_typeid is not None else None
        expected_bohr = R0_OO_BOHR if typ == 0 else R0_NN_BOHR if typ == 1 else None
        expected_nm = expected_bohr * BOHR_TO_NM if expected_bohr else None
        print(f"  bond {b}: length = {L_raw:.4f} (raw), {L_nm:.4f} nm  ", end="")
        if expected_nm is not None:
            print(f"(expected ~{expected_nm:.4f} nm for {'O-O' if typ == 0 else 'N-N'})")
        else:
            print()
    if n_mol > n_show:
        print("  ...")
    print()
    print("--- Overlap check ---")
    # Simple check: any pair (not necessarily bonded) closer than threshold in nm
    pos_nm = pos * BOHR_TO_NM
    box_nm = box * BOHR_TO_NM
    threshold = args.overlap_threshold
    n_overlap = 0
    for i in range(N):
        for j in range(i + 1, N):
            dr = pos_nm[j] - pos_nm[i]
            d = np.sqrt(np.dot(min_image_disp(dr, box_nm), min_image_disp(dr, box_nm)))
            if d < threshold:
                n_overlap += 1
                if n_overlap <= 5:
                    print(f"  pair ({i}, {j}): distance = {d:.4f} nm < {threshold} nm")
    if n_overlap == 0:
        print(f"  No pairs closer than {threshold} nm (OK).")
    else:
        print(f"  Total pairs closer than {threshold} nm: {n_overlap}")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
