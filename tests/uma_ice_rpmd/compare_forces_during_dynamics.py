#!/usr/bin/env python3
"""Compare UMA forces at positions extracted from LAMMPS trajectory.

Tests whether the OpenMM UMA force computation (direct AtomicData constructor)
produces the same forces as the LAMMPS-style computation for configurations
encountered during dynamics, not just the initial minimized structure.

Usage:
    python compare_forces_during_dynamics.py \
        --dump lammps/dump.pipeline_1ps_02fs.lammpstrj \
        --data pipeline_out/shared_data.ice_uma \
        --frames 0 20 50 100
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import torch
from ase.geometry import wrap_positions

_SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPT_DIR))
from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data


def read_dump_frame(dump_path: str, frame_idx: int):
    """Read a specific frame from a LAMMPS dump file. Returns (step, pos, box)."""
    with open(dump_path) as f:
        current_frame = -1
        step = None
        n_atoms = 0
        box_bounds = []
        atoms = []
        reading_atoms = False
        expect_step = False
        expect_natoms = False

        for line in f:
            line = line.strip()
            if line == "ITEM: TIMESTEP":
                if current_frame == frame_idx and atoms:
                    break
                current_frame += 1
                reading_atoms = False
                expect_step = True
                atoms = []
                box_bounds = []
                step = None
            elif expect_step:
                expect_step = False
                step = int(line)
                expect_natoms = False
            elif current_frame == frame_idx:
                if "ITEM: NUMBER" in line:
                    expect_natoms = True
                elif expect_natoms:
                    n_atoms = int(line)
                    expect_natoms = False
                elif "ITEM: BOX" in line:
                    box_bounds = []
                elif len(box_bounds) < 3 and not line.startswith("ITEM:"):
                    parts = line.split()
                    if len(parts) >= 2:
                        box_bounds.append((float(parts[0]), float(parts[1])))
                elif "ITEM: ATOMS" in line:
                    reading_atoms = True
                elif reading_atoms:
                    parts = line.split()
                    if len(parts) >= 6:
                        atoms.append(
                            (int(parts[0]), int(parts[1]), float(parts[2]),
                             float(parts[3]), float(parts[4]), float(parts[5]))
                        )

    if not atoms:
        raise ValueError(f"Frame {frame_idx} not found in {dump_path}")

    atoms.sort(key=lambda a: a[0])
    pos = np.array([[a[3], a[4], a[5]] for a in atoms], dtype=np.float64)
    lx = box_bounds[0][1] - box_bounds[0][0]
    ly = box_bounds[1][1] - box_bounds[1][0]
    lz = box_bounds[2][1] - box_bounds[2][0]
    cell = np.array([[lx, 0, 0], [0, ly, 0], [0, 0, lz]], dtype=np.float64)
    return step, pos, cell


def compute_forces_direct(pos_ang, cell, atomic_numbers, predict_unit, task_name,
                          charge=0, spin=1, device="cuda"):
    """Compute forces using direct AtomicData construction (OpenMM-style)."""
    from fairchem.core.datasets.atomic_data import AtomicData

    pos_wrapped = wrap_positions(pos_ang, cell=cell, pbc=True, eps=0)
    n = len(pos_wrapped)

    data = AtomicData(
        pos=torch.tensor(pos_wrapped, dtype=torch.float32),
        atomic_numbers=torch.tensor(atomic_numbers),
        cell=torch.tensor(cell, dtype=torch.float32).unsqueeze(0),
        pbc=torch.tensor([True, True, True], dtype=torch.bool).unsqueeze(0),
        natoms=torch.tensor([n], dtype=torch.long),
        edge_index=torch.empty((2, 0), dtype=torch.long),
        cell_offsets=torch.empty((0, 3), dtype=torch.float32),
        nedges=torch.tensor([0], dtype=torch.long),
        charge=torch.LongTensor([charge]),
        spin=torch.LongTensor([spin]),
        fixed=torch.zeros(n, dtype=torch.long),
        tags=torch.zeros(n, dtype=torch.long),
        batch=torch.zeros(n, dtype=torch.long),
        dataset=[task_name],
    )

    with torch.no_grad():
        data_dev = data.to(device)
        pred = predict_unit.predict(data_dev)
        energy = float(pred["energy"].cpu().item())
        forces = pred["forces"].cpu().numpy().copy()

    return energy, forces


def compute_forces_lammps_style(pos_ang, cell, atomic_numbers, predict_unit, task_name,
                                charge=0, spin=1, device="cuda"):
    """Compute forces using atomic_data_from_lammps_data (LAMMPS callback style)."""
    from fairchem.lammps.lammps_fc import atomic_data_from_lammps_data

    pos_wrapped = wrap_positions(pos_ang, cell=cell, pbc=True, eps=0)
    n = len(pos_wrapped)

    cell_t = torch.tensor(cell, dtype=torch.float32).unsqueeze(0)
    data = atomic_data_from_lammps_data(
        pos_wrapped, atomic_numbers, n, cell_t,
        [True, True, True], task_name,
        charge=charge, spin=spin,
    )

    with torch.no_grad():
        data_dev = data.to(device)
        pred = predict_unit.predict(data_dev)
        energy = float(pred["energy"].cpu().item())
        forces = pred["forces"].cpu().numpy().copy()

    return energy, forces


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", type=Path, required=True)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--frames", type=int, nargs="+", default=[0, 10, 20, 50, 100])
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--model", default="uma-s-1p1")
    args = ap.parse_args()

    _, _, Z = parse_lammps_data(args.data)
    atomic_numbers = Z.astype(np.int64)

    from fairchem.core.calculate import pretrained_mlip
    predict_unit = pretrained_mlip.get_predict_unit(args.model, device=args.device)

    print(f"{'Frame':>6s} {'Step':>8s} {'E_direct':>14s} {'E_lammps':>14s} {'dE(eV)':>12s} "
          f"{'max_dF':>10s} {'mean_dF':>10s} {'max_F':>10s}")

    for fi in args.frames:
        try:
            step, pos, cell = read_dump_frame(str(args.dump), fi)
        except ValueError as e:
            print(f"{fi:6d} {str(e)}")
            continue

        e_d, f_d = compute_forces_direct(pos, cell, atomic_numbers, predict_unit, "omol",
                                         device=args.device)
        e_l, f_l = compute_forces_lammps_style(pos, cell, atomic_numbers, predict_unit, "omol",
                                               device=args.device)

        df = np.abs(f_d - f_l)
        fmag = np.linalg.norm(f_d, axis=1)
        print(f"{fi:6d} {step:8d} {e_d:14.6f} {e_l:14.6f} {e_d - e_l:12.8f} "
              f"{df.max():10.6f} {df.mean():10.6f} {fmag.max():10.4f}")


if __name__ == "__main__":
    main()
