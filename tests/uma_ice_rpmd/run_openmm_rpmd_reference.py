#!/usr/bin/env python3
"""
Run OpenMM RPMD reference for i-PI+LAMMPS benchmark comparison.

Parameters match i-PI input: 243 K, PILE thermostat, 0.5 fs timestep.
Default supercell 2×2×2 (64 molecules), 4 beads. Structure from ``generate_ice_ih`` via ``--nx/--ny/--nz``.

Uses the same LAMMPS data file as the i-PI+LAMMPS run (via conversion to
XYZ with convert_lammps_to_ipi_xyz.py) so both engines start from an
identical initial structure. If the data file is missing, it is built automatically with
``--nx/--ny/--nz`` (default 2×2×2 → 64 molecules).

Produces ice_order_openmm_rpmd.csv for plot_rpmd_comparison.py.

Usage:
  cd tests/uma_ice_rpmd
  python run_openmm_rpmd_reference.py
  python run_openmm_rpmd_reference.py --nx 2 --ny 2 --nz 2 --beads 4 --steps 100   # short test
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_LAMMPS_DIR = _SCRIPT_DIR / "lammps"


def _ensure_lammps_data(nx: int, ny: int, nz: int) -> Path:
    """Build LAMMPS data file if it doesn't exist. Return its path."""
    molecules = 8 * nx * ny * nz
    data = _LAMMPS_DIR / f"data.ice_uma_{molecules}"
    if data.is_file():
        return data
    builder = _LAMMPS_DIR / "build_lammps_ice_data.py"
    if not builder.is_file():
        raise FileNotFoundError(f"Missing builder: {builder}")
    subprocess.check_call(
        [
            sys.executable,
            str(builder),
            "--nx",
            str(nx),
            "--ny",
            str(ny),
            "--nz",
            str(nz),
            "-o",
            str(data),
        ],
        cwd=str(_SCRIPT_DIR),
    )
    return data


def _lammps_data_to_xyz(data_path: Path, molecules: int) -> Path:
    """Convert LAMMPS data to ASE-readable extended XYZ with cell info."""
    import sys as _sys
    _sys.path.insert(0, str(_SCRIPT_DIR))
    from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data
    from ase import Atoms
    from ase.io import write as ase_write

    pos_ang, cell, Z = parse_lammps_data(data_path)
    symbols = ["O" if z == 8 else "H" for z in Z]
    atoms = Atoms(symbols=symbols, positions=pos_ang, cell=cell, pbc=True)

    xyz_path = _SCRIPT_DIR / "pipeline_out" / f"init_openmm_rpmd_{molecules}.xyz"
    xyz_path.parent.mkdir(parents=True, exist_ok=True)
    ase_write(str(xyz_path), atoms, format="extxyz")
    return xyz_path


def main() -> None:
    ap = argparse.ArgumentParser(description="OpenMM RPMD reference (match i-PI params)")
    ap.add_argument("--steps", type=int, default=None, help="Override: MD steps (default: from --prod and --dt)")
    ap.add_argument("--prod", type=float, default=1.0, help="Production time in ps (default: 1.0)")
    ap.add_argument("--dt", type=float, default=0.1, help="Timestep in fs (default: 0.1)")
    ap.add_argument(
        "--nx",
        type=int,
        default=2,
        help="Ice Ih supercell replication along a (default 2; 2×2×2 → 64 molecules)",
    )
    ap.add_argument("--ny", type=int, default=2, help="Supercell along b (default 2)")
    ap.add_argument("--nz", type=int, default=2, help="Supercell along c (default 2)")
    ap.add_argument("--beads", type=int, default=4, help="RPMD beads (default 4, scaled from 8)")
    ap.add_argument("--seed", type=int, default=284759, help="Random seed for thermostat/velocities")
    ap.add_argument("--order-csv", type=Path, default=_SCRIPT_DIR / "pipeline_out" / "ice_order_openmm_rpmd.csv")
    ap.add_argument("--platform", default="cuda")
    ap.add_argument(
        "--ml-device",
        default=None,
        help="UMA compute device: cuda (default with CUDA platform), cpu (slower; use if 32 beads OOM on GPU).",
    )
    ap.add_argument(
        "--inference-turbo-rpmd",
        action="store_true",
        help="Forward to test_uma_ice_rpmd: Fairchem 'turbo' preset. Incompatible with batched RPMD on some stacks (merge_mole); prefer --optimize-inference-tf32-only for 32+ beads.",
    )
    ap.add_argument(
        "--optimize-inference-tf32-only",
        action="store_true",
        help="Forward to test_uma_ice_rpmd: TF32 without merge_mole (fits batched RPMD; use for high bead count / VRAM limits).",
    )
    ap.add_argument(
        "--rpmd-thermostat",
        type=str,
        default="pile-g",
        choices=["pile-g", "pile", "none"],
        help="OpenMM RPMD bath: pile-g = PILE_G Bussi centroid + PILE internal (default); "
        "pile = PILE on all modes; none = no thermostat",
    )
    ap.add_argument(
        "--rpmd-friction",
        type=float,
        default=1.0,
        help="PILE internal friction 1/ps (default: 1.0)",
    )
    ap.add_argument(
        "--rpmd-centroid-friction",
        type=float,
        default=0.5,
        help="PILE_G centroid (Bussi) coupling 1/ps (default: 0.5)",
    )
    args = ap.parse_args()

    molecules = 8 * args.nx * args.ny * args.nz
    # Ensure identical initial structure: LAMMPS data → XYZ → --input
    data_path = _ensure_lammps_data(args.nx, args.ny, args.nz)
    xyz_path = _lammps_data_to_xyz(data_path, molecules)
    print(
        f"Using shared structure (ice Ih {args.nx}x{args.ny}x{args.nz} = {molecules} mol): "
        f"{data_path} → {xyz_path}"
    )

    dt_fs = args.dt
    # test_uma_ice_rpmd uses --prod (ps), not step count; derive prod from --steps when given
    if args.steps is not None:
        prod_ps = args.steps * dt_fs / 1000.0
    else:
        prod_ps = args.prod
    report_ps = max(0.01, min(0.1, max(prod_ps / 5.0, 0.02)))
    cmd = [
        sys.executable,
        str(_SCRIPT_DIR / "test_uma_ice_rpmd.py"),
        "--input", str(xyz_path),
        "--beads", str(args.beads),
        "--temperature", "243",
        "--dt", str(dt_fs),
        "--equil", "0",
        "--prod", str(prod_ps),
        "--order-csv", str(args.order_csv),
        "--report-interval", str(report_ps),
        "--order-every", "10",
        "--seed", str(args.seed),
        "--platform", args.platform,
        "--model", "uma-s-1p1-pythonforce-batch",
        "--rpmd-thermostat",
        args.rpmd_thermostat,
        "--rpmd-friction",
        str(args.rpmd_friction),
        "--rpmd-centroid-friction",
        str(args.rpmd_centroid_friction),
    ]
    if args.ml_device:
        cmd.extend(["--ml-device", args.ml_device])
    if args.inference_turbo_rpmd:
        cmd.append("--inference-turbo-rpmd")
    if args.optimize_inference_tf32_only:
        cmd.append("--optimize-inference-tf32-only")
    args.order_csv.parent.mkdir(parents=True, exist_ok=True)
    subprocess.check_call(cmd, cwd=str(_SCRIPT_DIR))
    print(f"Order CSV: {args.order_csv}")


if __name__ == "__main__":
    main()
