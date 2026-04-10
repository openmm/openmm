#!/usr/bin/env python3
"""
Post-process LAMMPS custom dump → ice order CSV (same schema as OpenMM ice_order.csv).

Reads orthorhombic box from data.ice_uma, parses dump (id type mass x y z), builds
positions in id order with type 1=O/8, 2=H/1, and calls ice_order_metrics per frame.

Use ``--thermo-log`` with the LAMMPS ``-log`` file from ``run_lammps_uma_ice.py`` (default
``lammps/lammps_uma_run.log``) to fill **T_K** and **PE_kj_mol** by matching **Step** to each
dump frame (same schema as OpenMM order CSVs).

Usage:
  cd tests/uma_ice_rpmd
  python lammps_order_from_dump.py --data lammps/data.ice_uma --dump lammps/dump.openmm_match_long.lammpstrj -o ice_order_lammps.csv
  python lammps_order_from_dump.py --data lammps/data.ice_uma --dump lammps/dump.foo.lammpstrj -o ice_order_lammps.csv --thermo-log lammps/lammps_uma_run.log --dt-fs 0.1
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data
from ice_order_parameters import ice_order_metrics
from lammps_thermo_log import parse_lammps_thermo_log


def _read_box_orth_ang(data_path: Path) -> np.ndarray:
    """(Lx, Ly, Lz) in Å from data file (orthorhombic diagonal)."""
    _, cell, _ = parse_lammps_data(data_path)
    return np.array([cell[0, 0], cell[1, 1], cell[2, 2]], dtype=np.float64)


def _iter_frames(dump_path: Path, dt_fs: float = 1.0):
    """Yield (step, time_ps, pos_ang (n,3), z (n,)) per frame. Dump must have id type [mass] x y z."""
    lines = dump_path.read_text().splitlines()
    n_lines = len(lines)
    i = 0
    while i < n_lines:
        if not lines[i].startswith("ITEM: TIMESTEP"):
            i += 1
            continue
        step = int(lines[i + 1])
        i += 2
        if i >= n_lines or not lines[i].startswith("ITEM: NUMBER OF ATOMS"):
            continue
        nat = int(lines[i + 1])
        i += 2
        while i < n_lines and not lines[i].startswith("ITEM: ATOMS"):
            i += 1
        if i >= n_lines:
            break
        header = lines[i].split()[2:]
        idx_id = header.index("id")
        idx_type = header.index("type")
        idx_x = header.index("x")
        i += 1
        rows = []
        for _ in range(nat):
            p = lines[i].split()
            aid = int(p[idx_id])
            typ = int(p[idx_type])
            x, y, z = float(p[idx_x]), float(p[idx_x + 1]), float(p[idx_x + 2])
            rows.append((aid, typ, x, y, z))
            i += 1
        rows.sort(key=lambda r: r[0])
        pos_ang = np.array([[r[2], r[3], r[4]] for r in rows], dtype=np.float64)
        z = np.array([8 if r[1] == 1 else 1 for r in rows], dtype=np.int64)
        time_ps = step * dt_fs / 1000.0
        yield step, time_ps, pos_ang, z


def main() -> None:
    ap = argparse.ArgumentParser(description="LAMMPS dump → ice order CSV (OpenMM schema)")
    ap.add_argument("--data", type=Path, default=_SCRIPT_DIR / "lammps" / "data.ice_uma", help="LAMMPS data file for box")
    ap.add_argument("--dump", type=Path, required=True, help="LAMMPS custom lammpstrj (id type mass x y z)")
    ap.add_argument("-o", "--output", type=Path, default=None, help="Output CSV (default: dump stem + _order.csv)")
    ap.add_argument("--every", type=int, default=1, help="Write CSV row every N frames (default 1)")
    ap.add_argument("--dt-fs", type=float, default=1.0, help="Timestep in fs (metal 0.001 = 1 fs)")
    ap.add_argument(
        "--thermo-log",
        type=Path,
        default=None,
        help=(
            "LAMMPS -log file with thermo (Step, Temp, PotEng). "
            "Fills T_K and PE_kj_mol columns by matching timestep to dump frames."
        ),
    )
    args = ap.parse_args()

    if not args.dump.is_file():
        print(f"Missing dump file: {args.dump}", file=sys.stderr)
        sys.exit(1)
    if not args.data.is_file():
        print(f"Missing data file: {args.data}", file=sys.stderr)
        sys.exit(1)

    box_orth = _read_box_orth_ang(args.data)
    out_path = args.output or args.dump.with_name(args.dump.stem + "_order.csv")

    thermo_by_step: dict[int, tuple[float, float]] = {}
    if args.thermo_log is not None:
        thermo_by_step = parse_lammps_thermo_log(args.thermo_log)
        if not thermo_by_step and args.thermo_log.is_file():
            print(
                f"Warning: no thermo rows parsed from {args.thermo_log} (T_K/PE_kj_mol empty).",
                file=sys.stderr,
            )

    with out_path.open("w") as f:
        f.write("step,time_ps,T_K,PE_kj_mol,q6_mean,q6_std,q_tet_mean,q_tet_std\n")
        for frame_idx, (step, time_ps, pos_ang, z) in enumerate(_iter_frames(args.dump, args.dt_fs)):
            if frame_idx % args.every != 0:
                continue
            res = ice_order_metrics(pos_ang, box_orth, z=z)
            t_k = pe_kj = ""
            if step in thermo_by_step:
                T_k, pe = thermo_by_step[step]
                t_k = f"{T_k:.4f}"
                pe_kj = f"{pe:.6f}"
            f.write(
                f"{step},{time_ps:.6f},{t_k},{pe_kj},{res.q6_mean:.6f},{res.q6_std:.6f},"
                f"{res.q_tet_mean:.6f},{res.q_tet_std:.6f}\n"
            )
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
