#!/usr/bin/env python3
"""
Compute Q6 order parameter from i-PI trajectory (centroid or single-bead xyz).

i-PI xyz format per frame:
  Line 1: natoms
  Line 2: # [CELL(abcABC): a b c alpha beta gamma] or other comment
  Lines 3+: symbol x y z (O and H)

When --beads > 1, reads traj_0 .. traj_{N-1}, averages positions for centroid
order (matches OpenMM RPMD centroid convention). Derives bead paths from --traj
(e.g. ipi/ice__i-pi.traj_0.xyz -> ipi/ice__i-pi.traj_{0..N-1}.xyz).

Output CSV matches lammps_order_from_dump schema for plot_rpmd_comparison.py.

Usage:
  cd tests/uma_ice_rpmd
  python ipi_order_from_traj.py --traj ipi/ice_traj.xyz -o pipeline_out/ice_order_ipi_rpmd.csv
  python ipi_order_from_traj.py --traj ipi/ice__i-pi.traj_0.xyz --beads 4 -o pipeline_out/ice_order_ipi_rpmd.csv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from ice_order_parameters import ice_order_metrics


def _abc_to_box_orth(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Convert a,b,c,angles (deg) to (Lx,Ly,Lz) for orthorhombic or triclinic.
    For orthogonal (alpha=beta=gamma=90): Lx=a, Ly=b, Lz=c.
    """
    if abs(alpha - 90) < 0.1 and abs(beta - 90) < 0.1 and abs(gamma - 90) < 0.1:
        return np.array([a, b, c], dtype=np.float64)
    # General triclinic: approximate orthorhombic box lengths
    rad = np.pi / 180.0
    cos_a = np.cos(alpha * rad)
    cos_b = np.cos(beta * rad)
    cos_g = np.cos(gamma * rad)
    vol = a * b * c * np.sqrt(
        1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g
    )
    # Use lengths as rough Lx,Ly,Lz (for orthorhombic projection)
    return np.array([a, b, c], dtype=np.float64)


BOHR_TO_ANG = 0.529177249  # 1 Bohr in Angstrom


def _traj_paths_from_prefix(traj_path: Path, beads: int) -> list[Path]:
    """Infer bead trajectory paths from e.g. ipi/ice__i-pi.traj_0.xyz -> traj_0..traj_{beads-1}."""
    stem = traj_path.stem  # e.g. ice__i-pi.traj_0
    base = re.sub(r"_\d+$", "", stem)  # e.g. ice__i-pi.traj
    parent = traj_path.parent
    suffix = traj_path.suffix
    return [parent / f"{base}_{b}{suffix}" for b in range(beads)]


def _iter_ipi_xyz_frames(traj_path: Path):
    """Yield (step, time_ps, pos_ang, box_orth, z) per frame.
    step/time parsed from comment if present; else frame index.
    Positions/cell in Bohr (atomic_unit) are converted to Angstrom.
    """
    cell_re = re.compile(
        r"#\s*CELL[\(\{\[]?abcABC[\)\}\]]?:\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)"
        r"\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)"
    )
    step_re = re.compile(r"Step:\s*(\d+)")
    time_re = re.compile(r"time\{\w+\}:\s*([\d.eE+-]+)")
    atomic_unit_re = re.compile(r"positions\{atomic_unit\}")

    text = traj_path.read_text()
    lines = text.splitlines()
    i = 0
    frame_idx = 0
    while i < len(lines):
        try:
            natoms = int(lines[i])
        except (ValueError, IndexError):
            i += 1
            continue
        i += 1
        if i >= len(lines):
            break
        comment = lines[i]
        i += 1

        step = frame_idx
        time_ps = float(frame_idx) * 0.5 / 1000.0  # default 0.5 fs timestep
        m = step_re.search(comment)
        if m:
            step = int(m.group(1))
        m = time_re.search(comment)
        if m:
            time_ps = float(m.group(1))
        elif "time" in comment.lower():
            # i-PI may use time{au}: val; convert au to ps if needed
            pass  # keep default

        box_orth = np.array([15.6, 14.7, 18.1], dtype=np.float64)  # fallback from init.xyz
        m = cell_re.search(comment)
        if m:
            a, b, c = float(m.group(1)), float(m.group(2)), float(m.group(3))
            alpha = float(m.group(4)) if len(m.groups()) >= 6 else 90.0
            beta = float(m.group(5)) if len(m.groups()) >= 6 else 90.0
            gamma = float(m.group(6)) if len(m.groups()) >= 6 else 90.0
            box_orth = _abc_to_box_orth(a, b, c, alpha, beta, gamma)

        in_bohr = bool(atomic_unit_re.search(comment)) or "atomic_unit" in comment.lower()
        # i-PI often outputs cell in Å, positions in Bohr; convert positions only

        pos_list = []
        z_list = []
        for _ in range(natoms):
            if i >= len(lines):
                break
            parts = lines[i].split()
            i += 1
            if len(parts) >= 4:
                sym = parts[0]
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                pos_list.append([x, y, z])
                z_list.append(8 if sym.upper() == "O" else 1)

        if len(pos_list) == natoms:
            pos_ang = np.array(pos_list, dtype=np.float64)
            if in_bohr:
                pos_ang = pos_ang * BOHR_TO_ANG
            z = np.array(z_list, dtype=np.int64)
            yield step, time_ps, pos_ang, box_orth, z
        frame_idx += 1


def _mic_centroid_ang(positions_per_bead: list[np.ndarray], box_orth: np.ndarray) -> np.ndarray:
    """MIC-aware centroid averaging over RPMD beads (positions in Angstrom).

    Uses bead 0 as reference; folds displacements via minimum-image convention
    before averaging, preventing PBC-wrap artifacts at cell boundaries.
    """
    ref = positions_per_bead[0]
    accum = np.zeros_like(ref)
    for b in range(len(positions_per_bead)):
        delta = positions_per_bead[b] - ref
        for k in range(3):
            delta[:, k] -= box_orth[k] * np.round(delta[:, k] / box_orth[k])
        accum += delta
    return ref + accum / len(positions_per_bead)


def _iter_ipi_xyz_frames_centroid(traj_path: Path, beads: int):
    """Yield (step, time_ps, centroid_pos_ang, box_orth, z) per frame.
    Reads traj_0..traj_{beads-1}, averages positions for centroid order
    using MIC-aware averaging to handle PBC correctly.
    """
    paths = _traj_paths_from_prefix(traj_path, beads)
    for p in paths:
        if not p.is_file():
            raise FileNotFoundError(f"Bead trajectory not found: {p}")
    bead_frames: list[list[tuple[int, float, np.ndarray, np.ndarray, np.ndarray]]] = []
    for p in paths:
        bead_frames.append(list(_iter_ipi_xyz_frames(p)))
    n_frames = min(len(f) for f in bead_frames)
    if n_frames == 0:
        print(f"Warning: no frames parsed from bead trajectories", file=sys.stderr)
    for i in range(n_frames):
        step, time_ps, _foo, box_orth, z = bead_frames[0][i]
        poses = [bead_frames[b][i][2] for b in range(beads)]
        centroid = _mic_centroid_ang(poses, box_orth)
        yield step, time_ps, centroid, box_orth, z


def main() -> None:
    ap = argparse.ArgumentParser(description="i-PI trajectory → ice order CSV")
    ap.add_argument("--traj", type=Path, default=_SCRIPT_DIR / "ipi" / "ice_traj.xyz")
    ap.add_argument("-o", "--output", type=Path, default=None)
    ap.add_argument("--every", type=int, default=1, help="Write every N frames")
    ap.add_argument("--dt-fs", type=float, default=0.5, help="Timestep in fs for time_ps")
    ap.add_argument(
        "--beads",
        type=int,
        default=1,
        help="RPMD beads; when >1, read traj_0..traj_{N-1} and use centroid positions",
    )
    args = ap.parse_args()

    if not args.traj.is_file():
        print(f"Missing trajectory: {args.traj}", file=sys.stderr)
        sys.exit(1)

    out_path = args.output or args.traj.with_name(args.traj.stem + "_order.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.beads <= 1:
        frame_iter = _iter_ipi_xyz_frames(args.traj)
    else:
        frame_iter = _iter_ipi_xyz_frames_centroid(args.traj, args.beads)

    with out_path.open("w") as f:
        f.write("step,time_ps,T_K,PE_kj_mol,q6_mean,q6_std,q_tet_mean,q_tet_std\n")
        for frame_idx, (step, time_ps, pos_ang, box_orth, z) in enumerate(frame_iter):
            if frame_idx % args.every != 0:
                continue
            res = ice_order_metrics(pos_ang, box_orth, z=z)
            f.write(
                f"{step},{time_ps:.6f},,,{res.q6_mean:.6f},{res.q6_std:.6f},"
                f"{res.q_tet_mean:.6f},{res.q_tet_std:.6f}\n"
            )
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
