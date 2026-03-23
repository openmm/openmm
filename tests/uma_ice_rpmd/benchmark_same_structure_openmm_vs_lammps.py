#!/usr/bin/env python3
"""
Apples-to-apples UMA benchmark: same atoms, same box, same positions as LAMMPS data file.

Reads tests/uma_ice_rpmd/lammps/data.ice_uma (or path you pass) — the same file
run_lammps_uma_ice.py uses — so the structure matches LAMMPS exactly.

Times three paths (all call the same Fairchem uma-s-1p1 predict on CUDA):

1. **LAMMPS_pipeline** — wrap_positions + atomic_data_from_lammps_data + predict
   (mirrors fairchem.lammps.FixExternalCallback).

2. **OpenMM_pipeline** — molecular water wrap + AtomicData.from_ase once (warmup),
   then per-iter: update pos/cell + atomicdata_list_to_batch + predict
   (mirrors openmmml batch path, 1 bead).

3. **predict_only** — single frozen AtomicData on GPU, predict only (ML lower bound).

Usage:
  python benchmark_same_structure_openmm_vs_lammps.py
  python benchmark_same_structure_openmm_vs_lammps.py --data lammps/data.ice_uma --reps 80
"""
from __future__ import annotations

import argparse
import statistics
import sys
import time
from pathlib import Path

import numpy as np
import torch

_SCRIPT_DIR = Path(__file__).resolve().parent
_LAMMPS_DIR = _SCRIPT_DIR / "lammps"


def parse_lammps_data(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse LAMMPS data: orthorhombic box, Atoms section.
    Returns pos_angstrom (n,3), cell (3,3), atomic_numbers (n,) int64.
    """
    text = path.read_text()
    lines = text.splitlines()
    lx = ly = lz = 0.0
    xy = xz = yz = 0.0
    n_atoms = 0
    in_atoms = False
    coords = []
    types = []
    for line in lines:
        s = line.strip()
        if s.startswith("Atoms"):
            in_atoms = True
            continue
        if in_atoms:
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 5:
                continue
            try:
                mol_id = int(parts[0])
                typ = int(parts[1])
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            except ValueError:
                continue
            coords.append([x, y, z])
            types.append(typ)
    # parse box from header
    for line in lines:
        if "xlo xhi" in line and "ylo" not in line:
            a = line.split()
            lx = float(a[1]) - float(a[0])
        elif "ylo yhi" in line:
            a = line.split()
            ly = float(a[1]) - float(a[0])
        elif "zlo zhi" in line:
            a = line.split()
            lz = float(a[1]) - float(a[0])
        elif "xy xz yz" in line:
            a = line.split()
            xy, xz, yz = float(a[0]), float(a[1]), float(a[2])
        elif "atoms" in line and line.strip()[0].isdigit():
            n_atoms = int(line.split()[0])
    pos = np.array(coords, dtype=np.float64)
    assert pos.shape[0] == n_atoms, f"atom count {pos.shape[0]} vs header {n_atoms}"
    cell = np.array([[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]], dtype=np.float64)
    # type 1 = O (8), 2 = H (1)
    Z = np.array([8 if t == 1 else 1 for t in types], dtype=np.int64)
    return pos, cell, Z


def main() -> None:
    ap = argparse.ArgumentParser(description="Same-structure UMA: LAMMPS vs OpenMM pipeline timing")
    ap.add_argument(
        "--data",
        type=Path,
        default=_LAMMPS_DIR / "data.ice_uma",
        help="LAMMPS data file (same as LAMMPS run)",
    )
    ap.add_argument("--model", default="uma-s-1p1", help="Fairchem model name")
    ap.add_argument("--device", default="cuda", choices=["cuda", "cpu"])
    ap.add_argument("--warmup", type=int, default=15)
    ap.add_argument("--reps", type=int, default=50)
    args = ap.parse_args()

    if not args.data.is_file():
        print(
            f"Missing {args.data}. Build with: python {_LAMMPS_DIR / 'build_lammps_ice_data.py'} "
            f"--nx 2 --ny 2 --nz 4 -o {args.data}"
        )
        sys.exit(1)

    pos, cell, atomic_numbers = parse_lammps_data(args.data)
    n = pos.shape[0]
    periodicity = np.array([True, True, True])
    task_name = "omol"

    from fairchem.core.calculate import pretrained_mlip
    from fairchem.core.datasets.atomic_data import AtomicData, atomicdata_list_to_batch
    from fairchem.lammps.lammps_fc import atomic_data_from_lammps_data, restricted_cell_from_lammps_box
    from ase.geometry import wrap_positions
    from ase import Atoms

    # Import openmmml molecular wrap (same as MD)
    from openmmml.models.umapotential_pythonforce_batch import _wrap_water_molecules_per_bead

    device = args.device
    if device == "cuda" and not torch.cuda.is_available():
        print("CUDA unavailable; use --device cpu")
        sys.exit(1)

    predictor = pretrained_mlip.get_predict_unit(args.model, device=device)
    cell_torch = torch.from_numpy(np.asarray(cell, dtype=np.float32)).float().unsqueeze(0)

    # --- LAMMPS box tuple for restricted_cell_from_lammps_box
    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    boxlo = np.array([0.0, 0.0, 0.0])
    boxhi = np.array([lx, ly, lz])
    cell_lmp = restricted_cell_from_lammps_box(boxlo, boxhi, xy, yz, xz)

    def run_lammps_pipeline() -> None:
        x_wrapped = wrap_positions(pos, cell=cell, pbc=periodicity, eps=0)
        ad = atomic_data_from_lammps_data(
            x_wrapped,
            atomic_numbers,
            n,
            cell_lmp,
            periodicity,
            task_name,
            charge=0,
            spin=1,
        )
        ad = ad.to(device)
        predictor.predict(ad)

    def build_openmm_template():
        pos_mol = _wrap_water_molecules_per_bead(pos.copy(), cell)
        atoms = Atoms(numbers=atomic_numbers, positions=pos_mol, cell=cell, pbc=True)
        atoms.info["charge"] = 0
        atoms.info["spin"] = 1
        data = AtomicData.from_ase(
            atoms, task_name=task_name, r_edges=False, r_data_keys=["spin", "charge"]
        )
        data.sid = ["bead-0"]
        return data

    template = build_openmm_template()

    def run_openmm_pipeline() -> None:
        pos_mol = _wrap_water_molecules_per_bead(pos.copy(), cell)
        template.pos = torch.from_numpy(np.ascontiguousarray(pos_mol, dtype=np.float32)).float()
        template.cell = torch.from_numpy(np.ascontiguousarray(cell, dtype=np.float32)).float().unsqueeze(0)
        co = getattr(template, "cell_offsets", None)
        if co is not None and co.numel() > 0:
            template.cell_offsets = co.float()
        batch = atomicdata_list_to_batch([template])
        batch = batch.to(device, non_blocking=True)
        predictor.predict(batch)

    # predict_only: clone template to GPU once, only swap pos if needed (already set)
    batch_frozen = atomicdata_list_to_batch([build_openmm_template()])
    batch_frozen = batch_frozen.to(device)

    def run_predict_only() -> None:
        with torch.no_grad():
            predictor.predict(batch_frozen)

    def bench(name: str, fn, warmup: int, reps: int) -> list[float]:
        for _ in range(warmup):
            fn()
        if device == "cuda":
            torch.cuda.synchronize()
        times = []
        for _ in range(reps):
            t0 = time.perf_counter()
            fn()
            if device == "cuda":
                torch.cuda.synchronize()
            times.append((time.perf_counter() - t0) * 1000.0)
        return times

    print("=" * 72)
    print("Same structure UMA benchmark")
    print("=" * 72)
    print(f"  Data file:    {args.data}")
    print(f"  Atoms:        {n}")
    print(f"  Box (Å):      lx={lx:.4f} ly={ly:.4f} lz={lz:.4f}  (same as LAMMPS)")
    print(f"  Model:        {args.model}  device={device}")
    print(f"  Warmup / rep: {args.warmup} / {args.reps}")
    print()

    t_lmp = bench("LAMMPS_pipeline", run_lammps_pipeline, args.warmup, args.reps)
    t_omm = bench("OpenMM_pipeline", run_openmm_pipeline, args.warmup, args.reps)
    t_pred = bench("predict_only", run_predict_only, args.warmup, args.reps)

    def stat(xs: list[float]) -> str:
        return f"mean={statistics.mean(xs):.2f} ms  std={statistics.stdev(xs):.2f} ms  min={min(xs):.2f}  max={max(xs):.2f}"

    print(f"  LAMMPS_pipeline (wrap + AtomicData_lmp + predict):  {stat(t_lmp)}")
    print(f"  OpenMM_pipeline (mol_wrap + batch + predict):       {stat(t_omm)}")
    print(f"  predict_only (same graph, no Python prep):          {stat(t_pred)}")
    print()
    m_lmp, m_omm, m_pred = statistics.mean(t_lmp), statistics.mean(t_omm), statistics.mean(t_pred)
    print(f"  Interpretation: ~{m_pred:.0f} ms is almost entirely one UMA forward ({n} atoms).")
    print(f"  LAMMPS vs OpenMM pipeline prep differs by ~{abs(m_omm - m_lmp):.2f} ms (noise at this scale).")
    print(f"  OpenMM − LAMMPS pipeline: {m_omm - m_lmp:+.2f} ms per eval")
    print("=" * 72)


if __name__ == "__main__":
    main()
