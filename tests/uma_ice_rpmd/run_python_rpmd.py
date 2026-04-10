#!/usr/bin/env python3
"""
run_python_rpmd.py – Run RPMD with UMA, no MPI, no LAMMPS, no i-PI.

A single Python process manages all ring-polymer beads.  Forces are
evaluated once per timestep via a batched UMA forward pass (all beads in one
GPU call), so GPU memory scales with model size × 1, not model size × N_beads.

Memory comparison on an RTX 4070 (11.7 GB):
  MPI / LAMMPS approach : N_beads × ~1.1 GB model  → 32 beads needs ~35 GB
  This script           : 1 × ~1.1 GB model        → 32 beads fit in ~3–4 GB

Usage
-----
# Smoke test: 8 beads, 8-molecule ice, 20 steps
python run_python_rpmd.py \\
    --data  lammps/data.ice_uma_8 \\
    --beads 8 --temp 243 --dt 0.1 \\
    --equil 0.0 --prod 0.02

# Match OpenMM reference: 32 beads, 64-molecule ice, 10 ps production
python run_python_rpmd.py \\
    --data  lammps/data.ice_uma_64 \\
    --beads 32 --temp 243 --dt 0.1 \\
    --equil 2.0 --prod 10.0 \\
    --gamma-centroid 5e-4 \\
    --report 0.1 \\
    --out pipeline_out/python_rpmd_32b
"""
from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent


# ── Structure loading ──────────────────────────────────────────────────────────

def _parse_lammps_data(path: Path):
    """
    Minimal LAMMPS data file parser (atom_style atomic, orthogonal box).

    Returns
    -------
    positions : ndarray (n_atoms, 3) Å
    cell      : ndarray (3, 3) Å
    Z_list    : list[int]  atomic numbers (8=O, 1=H, 0=unknown)
    """
    with open(path) as fh:
        lines = fh.readlines()

    n_atoms   = 0
    xlo, xhi  = 0.0, 0.0
    ylo, yhi  = 0.0, 0.0
    zlo, zhi  = 0.0, 0.0
    mass_map: dict[int, float] = {}
    positions: list = []
    z_list:    list[int] = []
    section   = ""

    for raw in lines:
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        if "atoms" in s and "atom types" not in s:
            n_atoms = int(s.split()[0])
        elif s.endswith("xlo xhi"):
            xlo, xhi = map(float, s.split()[:2])
        elif s.endswith("ylo yhi"):
            ylo, yhi = map(float, s.split()[:2])
        elif s.endswith("zlo zhi"):
            zlo, zhi = map(float, s.split()[:2])
        elif s == "Masses":
            section = "masses"
        elif s == "Atoms":
            section = "atoms"
        elif section == "masses" and len(s.split()) == 2:
            t, m = int(s.split()[0]), float(s.split()[1])
            mass_map[t] = m
        elif section == "atoms" and len(s.split()) >= 5:
            parts  = s.split()
            atype  = int(parts[1])
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            positions.append([x, y, z])
            m = mass_map.get(atype, 0.0)
            if abs(m - 15.999) < 1.0:
                z_list.append(8)
            elif abs(m - 1.008) < 0.5:
                z_list.append(1)
            else:
                z_list.append(0)

    cell = np.array([
        [xhi - xlo, 0.0,       0.0],
        [0.0,       yhi - ylo, 0.0],
        [0.0,       0.0,       zhi - zlo],
    ])
    return np.array(positions), cell, z_list


def _load_structure(path: Path):
    """
    Load initial structure as an ase.Atoms object.
    Supports extended XYZ and LAMMPS data files.
    """
    from ase import Atoms
    from ase.io import read as ase_read

    p = Path(path)
    suffix = p.suffix.lower()

    # ASE can read many formats directly
    if suffix in (".xyz", ".extxyz"):
        return ase_read(str(p), format="extxyz" if "extxyz" in p.name else "xyz")

    try:
        at = ase_read(str(p))
        if at is not None:
            return at
    except Exception:
        pass

    # Custom parser for LAMMPS atomic data files
    pos, cell, Z = _parse_lammps_data(p)
    symbols = ["O" if z == 8 else ("H" if z == 1 else "X") for z in Z]
    return Atoms(symbols=symbols, positions=pos, cell=cell, pbc=True)


# ── Entry point ────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Pure-Python RPMD with UMA (no MPI, no LAMMPS, no i-PI).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Structure
    ap.add_argument(
        "--data", type=Path,
        default=_SCRIPT_DIR / "lammps" / "data.ice_uma_8",
        metavar="FILE",
        help="Initial structure: LAMMPS data file or extended XYZ",
    )

    # RPMD parameters
    ap.add_argument("--beads", type=int, default=8,    help="Number of ring-polymer beads")
    ap.add_argument("--temp",  type=float, default=243.0, help="Temperature [K]")
    ap.add_argument("--dt",    type=float, default=0.1,   help="Timestep [fs]")

    # Thermostat
    ap.add_argument(
        "--gamma-centroid", type=float, default=5.0e-4,
        metavar="GAMMA",
        help=(
            "Centroid Langevin friction [fs⁻¹]  "
            "(default: 5e-4 ≈ 0.5 ps⁻¹, matching OpenMM pile-g centroid). "
            "Internal modes use PILE-L critical damping automatically."
        ),
    )

    # Run length
    ap.add_argument(
        "--equil", type=float, default=1.0,
        help="Equilibration time [ps]. Use 0 for quick smoke tests.",
    )
    ap.add_argument("--prod",  type=float, default=2.0,  help="Production time [ps]")

    # ML model
    ap.add_argument("--model",  default="uma-s-1p1", help="FairChem model name")
    ap.add_argument("--device", default="cuda",       help="Compute device (cuda / cpu)")
    ap.add_argument("--task",   default="omol",       help="FairChem task name")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin",   type=int, default=1)

    # Performance
    ap.add_argument(
        "--chunk", type=int, default=None,
        metavar="N",
        help=(
            "Split bead batch into chunks of N for GPU memory control "
            "(default: all beads in one call). Useful for N_beads > 16 on <12 GB GPUs."
        ),
    )

    # Output
    ap.add_argument(
        "--report", type=float, default=0.1,
        help="Console/CSV reporting interval [ps]",
    )
    ap.add_argument(
        "--traj-every", type=int, default=10,
        help="Write centroid trajectory frame every N production steps",
    )
    ap.add_argument(
        "--out", type=Path,
        default=_SCRIPT_DIR / "pipeline_out" / "python_rpmd",
        help="Output prefix (directory/basename). Suffixes _centroid.xyz and _thermo.csv are appended.",
    )
    ap.add_argument("--seed", type=int, default=42, help="Random seed")
    ap.add_argument("-v", "--verbose", action="store_true")

    args = ap.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    log = logging.getLogger(__name__)

    # ── Import engine ──────────────────────────────────────────────────────────
    sys.path.insert(0, str(_SCRIPT_DIR))
    from rpmd_engine import RPMDConfig, RPMDEngine

    # ── Output files ───────────────────────────────────────────────────────────
    out = args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    traj_path   = Path(str(out) + "_centroid.xyz")
    thermo_path = Path(str(out) + "_thermo.csv")

    # ── Structure ──────────────────────────────────────────────────────────────
    log.info("Loading structure: %s", args.data)
    atoms = _load_structure(args.data)
    cell_diag = np.array(atoms.get_cell()).diagonal()
    log.info("  %d atoms, cell [%.3f, %.3f, %.3f] Å", len(atoms), *cell_diag)

    # ── Engine config ──────────────────────────────────────────────────────────
    cfg = RPMDConfig(
        n_beads        = args.beads,
        temperature    = args.temp,
        dt             = args.dt,
        model_name     = args.model,
        device         = args.device,
        task_name      = args.task,
        charge         = args.charge,
        spin           = args.spin,
        seed           = args.seed,
        gamma_centroid = args.gamma_centroid,
        chunk_size     = args.chunk,
    )

    log.info(
        "RPMD config: %d beads, T=%.1f K, dt=%.3f fs, "
        "γ_centroid=%.2e fs⁻¹, device=%s",
        cfg.n_beads, cfg.temperature, cfg.dt,
        cfg.gamma_centroid, cfg.device,
    )

    engine = RPMDEngine(atoms, cfg)

    # ── Initial forces ─────────────────────────────────────────────────────────
    engine.initialise_forces()

    # ── Helpers ────────────────────────────────────────────────────────────────
    dt_fs        = args.dt
    report_steps = max(1, round(args.report * 1000.0 / dt_fs))

    def _run_phase(n_steps: int, phase: str, write_output: bool = False) -> None:
        from ase import Atoms as _Atoms
        from ase.io import write as ase_write

        for local_step in range(1, n_steps + 1):
            engine.step()

            if local_step % report_steps == 0 or local_step == n_steps:
                t_ps  = engine.step_count * dt_fs * 1e-3
                E_pe  = engine.potential_energy_ev
                E_ke  = engine.kinetic_energy_ev
                E_sp  = engine.ring_spring_energy_ev
                T_k   = engine.kinetic_temperature_k
                log.info(
                    "[%s] step %6d  t=%.3f ps  "
                    "E_pe=%+.4f  E_ke=%.4f  E_spring=%.4f  T=%.1f K",
                    phase, engine.step_count, t_ps,
                    E_pe, E_ke, E_sp, T_k,
                )
                if write_output:
                    with open(thermo_path, "a", newline="") as fh:
                        csv.writer(fh).writerow([
                            engine.step_count,
                            f"{t_ps:.5f}",
                            f"{E_pe:.6f}",
                            f"{E_ke:.6f}",
                            f"{E_sp:.6f}",
                            f"{T_k:.3f}",
                        ])

            if write_output and local_step % args.traj_every == 0:
                cent = _Atoms(
                    symbols   = engine.symbols,
                    positions = engine.centroid_positions,
                    cell      = engine.cell,
                    pbc       = engine.pbc,
                )
                cent.info["step"] = engine.step_count
                cent.info["t_ps"] = f"{engine.step_count * dt_fs * 1e-3:.5f}"
                ase_write(str(traj_path), cent, format="extxyz", append=True)

    # ── CSV header ─────────────────────────────────────────────────────────────
    with open(thermo_path, "w", newline="") as fh:
        csv.writer(fh).writerow(
            ["step", "time_ps", "E_phys_eV", "E_kin_eV", "E_spring_eV", "T_inst_K"]
        )

    # ── Equilibration ──────────────────────────────────────────────────────────
    equil_steps = round(args.equil * 1000.0 / dt_fs)
    if equil_steps > 0:
        log.info(
            "=== Equilibration: %d steps (%.3f ps) ===", equil_steps, args.equil
        )
        _run_phase(equil_steps, phase="EQUIL", write_output=False)

    # ── Production ─────────────────────────────────────────────────────────────
    prod_steps = round(args.prod * 1000.0 / dt_fs)
    log.info("=== Production: %d steps (%.3f ps) ===", prod_steps, args.prod)
    _run_phase(prod_steps, phase="PROD", write_output=True)

    log.info("Done.")
    log.info("  Centroid trajectory : %s", traj_path)
    log.info("  Thermo data         : %s", thermo_path)


if __name__ == "__main__":
    main()
