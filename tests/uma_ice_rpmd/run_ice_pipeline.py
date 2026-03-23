#!/usr/bin/env python3
"""
End-to-end ice pipeline: build Ih supercell (``--nx/--ny/--nz``) → 1 ps TIP4P/2005f (OpenMM) → 1 ps UMA (OpenMM)
→ 1 ps UMA (LAMMPS) → plot order parameters. Runs strictly sequentially (no parallel sims)
to avoid GPU OOM. Same thermostat (Langevin), same timestep (--dt-fs), same T and friction for both.

Units (parity):
  Temperature: K (OpenMM unit.kelvin; LAMMPS metal = K).
  Langevin: OpenMM friction in 1/ps; LAMMPS fix langevin damp in time units (metal → ps), damp = 1/friction (e.g. friction 10 → damp 0.1 ps).
  Timestep: fs; LAMMPS metal timestep in ps = dt_fs * 0.001.
Order/trajectory every 10 fs.

Usage:
  cd tests/uma_ice_rpmd
  # Default: proton-ordered ice Ih 2×2×2 (same as RPMD benchmark) → 64 molecules
  python run_ice_pipeline.py --out-dir pipeline_out
  # Other supercells:
  python run_ice_pipeline.py --out-dir pipeline_out --nx 1 --ny 2 --nz 2
  # Background:
  nohup python run_ice_pipeline.py --out-dir pipeline_out >> pipeline.log 2>&1 &
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_LAMMPS_DIR = _SCRIPT_DIR / "lammps"


def _build_data_supercell(work_dir: Path, nx: int, ny: int, nz: int) -> Path:
    """Proton-ordered ice Ih from ``generate_ice_ih`` replication (same as RPMD ``--nx/--ny/--nz``)."""
    molecules = 8 * nx * ny * nz
    data_path = work_dir / "lammps" / f"data.ice_uma_{molecules}"
    data_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(_LAMMPS_DIR / "build_lammps_ice_data.py"),
        "--nx",
        str(nx),
        "--ny",
        str(ny),
        "--nz",
        str(nz),
        "-o",
        str(data_path),
    ]
    subprocess.check_call(cmd, cwd=str(work_dir))
    if not data_path.is_file():
        print(f"Build did not create {data_path}", file=sys.stderr)
        sys.exit(1)
    n_atoms = 0
    for line in data_path.read_text().splitlines():
        parts = line.strip().split()
        if len(parts) == 2 and parts[1] == "atoms":
            n_atoms = int(parts[0])
            break
    if n_atoms != molecules * 3:
        print(f"Expected {molecules * 3} atoms, got {n_atoms}", file=sys.stderr)
        sys.exit(1)
    print(f"Built {data_path} ({n_atoms} atoms, Ih supercell {nx}×{ny}×{nz})")
    return data_path


def _run_openmm_classical(
    work_dir: Path,
    out_dir: Path,
    steps: int,
    dt_fs: float,
    output_every: int,
    data_path: Path,
    platform: str = "cpu",
    variable_step: bool = False,
    cutoff_nm: float | None = None,
) -> None:
    """Run TIP4P/2005f OpenMM; wait for completion."""
    cmd = [
        sys.executable,
        str(_SCRIPT_DIR / "run_openmm_ice_classical_flex.py"),
        "--data", str(data_path),
        "--steps", str(steps),
        "--dt-fs", str(dt_fs),
        "--friction", "10",
        "--report-every", str(output_every),
        "--traj-every", str(output_every),
        "--order-every", str(output_every),
        "--order-csv", str(out_dir / "ice_order_classical.csv"),
        "--trajectory", str(out_dir / "traj_classical.dcd"),
        "--platform", platform,
    ]
    if variable_step:
        cmd.append("--variable-step")
    if cutoff_nm is not None:
        cmd.extend(["--cutoff", str(cutoff_nm)])
    subprocess.check_call(cmd, cwd=str(work_dir))
    print("OpenMM TIP4P/2005f done.")


def _run_openmm_uma(
    work_dir: Path,
    out_dir: Path,
    steps: int,
    dt_fs: float,
    output_every: int,
    data_path: Path,
    skip_minimize: bool = False,
    temperature: float = 243.0,
    friction: float = 10.0,
    integrator: str = "langevin-middle",
    no_cm_removal: bool = False,
) -> None:
    """Run UMA OpenMM; wait for completion. Temperature in K, friction in 1/ps (Langevin)."""
    data_path = data_path or (work_dir / "lammps" / "data.ice_uma")
    csv_suffix = f"_openmm_{integrator}" if integrator != "langevin-middle" else "_openmm"
    csv_name = f"ice_order_uma{csv_suffix.replace('-', '_')}.csv"
    traj_name = f"traj_uma{csv_suffix.replace('-', '_')}.dcd"
    cmd = [
        sys.executable,
        str(_SCRIPT_DIR / "run_openmm_ice_lammps_match.py"),
        "--data", str(data_path),
        "--steps", str(steps),
        "--dt-fs", str(dt_fs),
        "--temperature", str(temperature),
        "--friction", str(friction),
        "--report-every", str(output_every),
        "--traj-every", str(output_every),
        "--order-every", str(output_every),
        "--order-csv", str(out_dir / csv_name),
        "--trajectory", str(out_dir / traj_name),
        "--integrator", integrator,
    ]
    if skip_minimize:
        cmd.append("--skip-minimize")
    if no_cm_removal:
        cmd.append("--no-cm-removal")
    subprocess.check_call(cmd, cwd=str(work_dir))
    print(f"OpenMM UMA done (integrator={integrator}).")


def _run_lammps_uma(
    work_dir: Path,
    out_dir: Path,
    steps: int,
    output_every: int,
    dt_fs: float,
    data_path: Path | None = None,
    skip_minimize: bool = False,
    temperature: float = 243.0,
    friction: float = 10.0,
) -> None:
    """Run LAMMPS UMA then post-process dump to order CSV; wait for completion.
    LAMMPS timestep is set from dt_fs so total simulated time matches OpenMM (same steps * dt_fs).
    Temperature in K, friction in 1/ps; LAMMPS uses damp (ps) = 1/friction (units metal).
    """
    infile = _LAMMPS_DIR / "in.ice_uma_pipeline_1ps_02fs.lmp"
    data_path = data_path or (work_dir / "lammps" / "data.ice_uma")
    cmd = [
        sys.executable,
        str(_LAMMPS_DIR / "run_lammps_uma_ice.py"),
        "--infile", str(infile),
        "--data", str(data_path),
        "--dt-fs", str(dt_fs),
        "--temperature", str(temperature),
        "--friction", str(friction),
        "--steps", str(steps),
        "--dump-every", str(output_every),
        "--thermo-every", str(output_every),
    ]
    if skip_minimize:
        cmd.append("--skip-minimize")
    subprocess.check_call(cmd, cwd=str(work_dir))
    dump_path = work_dir / "lammps" / "dump.pipeline_1ps_02fs.lammpstrj"
    if not dump_path.is_file():
        print(f"LAMMPS did not produce {dump_path}", file=sys.stderr)
        sys.exit(1)
    # Use same dt_fs as LAMMPS run (set via --dt-fs above) so time_ps in CSV matches actual simulated time.
    thermo_log = work_dir / "lammps" / "lammps_uma_run.log"
    cmd2 = [
        sys.executable,
        str(_SCRIPT_DIR / "lammps_order_from_dump.py"),
        "--data", str(data_path),
        "--dump", str(dump_path),
        "-o", str(out_dir / "ice_order_uma_lammps.csv"),
        "--dt-fs", str(dt_fs),
        "--thermo-log", str(thermo_log),
    ]
    subprocess.check_call(cmd2, cwd=str(work_dir))
    print("LAMMPS UMA + order post-process done.")


def plot_order_vs_time(
    csv_classical: Path,
    csv_uma_openmm: Path,
    csv_uma_lammps: Path,
    out_png: Path,
) -> None:
    """Plot Q6 and q_tet vs time_ps for all three runs. Uses matplotlib Agg (headless)."""
    import csv as csv_module
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    def _read_time_q6_qtet(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        t, q6, qt = [], [], []
        with path.open() as f:
            r = csv_module.DictReader(f)
            for row in r:
                t.append(float(row["time_ps"]))
                q6.append(float(row["q6_mean"]))
                qt.append(float(row["q_tet_mean"]))
        return np.array(t), np.array(q6), np.array(qt)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
    for label, path, color in [
        ("TIP4P/2005f (OpenMM)", csv_classical, "C0"),
        ("UMA (OpenMM)", csv_uma_openmm, "C1"),
        ("UMA (LAMMPS)", csv_uma_lammps, "C2"),
    ]:
        if not path.is_file():
            continue
        time_ps, q6, q_tet = _read_time_q6_qtet(path)
        ax1.plot(time_ps, q6, label=label, color=color)
        ax2.plot(time_ps, q_tet, label=label, color=color)
    ax1.set_ylabel("Q6 (Steinhardt)")
    ax2.set_ylabel("q_tet (tetrahedral)")
    ax2.set_xlabel("Time (ps)")
    ax1.legend()
    ax2.legend()
    ax1.set_title("Ice order parameters")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close()
    print(f"Saved {out_png}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="End-to-end ice pipeline: structure → 1 ps × 3 (TIP4P/2005f, UMA OpenMM, UMA LAMMPS) → plot"
    )
    ap.add_argument("--out-dir", type=Path, default=_SCRIPT_DIR / "pipeline_out", help="Output directory for CSVs, trajectories, plot")
    ap.add_argument("--work-dir", type=Path, default=_SCRIPT_DIR, help="Work directory (lammps/ data and outputs)")
    ap.add_argument("--dt-fs", type=float, default=0.1, help="Timestep in fs for all runs (default 0.1 for TIP4P/2005f stability; 1 ps = 10000 steps)")
    ap.add_argument("--steps", type=int, default=None, help="MD steps per run (default: 1000/dt_fs for 1 ps, e.g. 5000 @ 0.2 fs)")
    ap.add_argument("--temperature", type=float, default=243.0, help="Temperature in K for NVT (OpenMM and LAMMPS, default 243)")
    ap.add_argument("--friction", type=float, default=10.0, help="Langevin friction in 1/ps (OpenMM and LAMMPS). LAMMPS damp (ps) = 1/friction (default 10 -> damp 0.1 ps)")
    ap.add_argument(
        "--nx",
        type=int,
        default=2,
        metavar="NX",
        help="Ice Ih supercell replication along a (default 2; with ny=nz=2 → 64 H2O). Uses generate_ice_ih replication.",
    )
    ap.add_argument(
        "--ny",
        type=int,
        default=2,
        metavar="NY",
        help="Supercell along b (default 2).",
    )
    ap.add_argument(
        "--nz",
        type=int,
        default=2,
        metavar="NZ",
        help="Supercell along c (default 2).",
    )
    ap.add_argument(
        "--skip-classical",
        action="store_true",
        help="Skip TIP4P/2005f (optional; use only if debugging).",
    )
    ap.add_argument("--skip-lammps", action="store_true", help="Skip LAMMPS run (still plot if OpenMM CSVs exist)")
    ap.add_argument(
        "--only-plot",
        action="store_true",
        help="Only run plotting from existing CSVs in out-dir (or from --order-csv-* paths)",
    )
    ap.add_argument("--order-csv-classical", type=Path, default=None)
    ap.add_argument("--order-csv-uma-openmm", type=Path, default=None)
    ap.add_argument("--order-csv-uma-lammps", type=Path, default=None)
    ap.add_argument(
        "--platform-classical",
        default="cuda",
        choices=["cpu", "cuda"],
        help="Platform for TIP4P/2005f run (default: cuda)",
    )
    ap.add_argument(
        "--variable-step-classical",
        action="store_true",
        help="Use VariableLangevinIntegrator for classical run (adaptive dt for stability)",
    )
    ap.add_argument(
        "--integrator",
        default="langevin-middle",
        choices=["langevin-middle", "bbk"],
        help="OpenMM UMA integration scheme: langevin-middle = BAOAB (default), "
        "bbk = LAMMPS-equivalent VV + Langevin (fix nve + fix langevin)",
    )
    ap.add_argument(
        "--no-cm-removal",
        action="store_true",
        help="Disable CMMotionRemover in OpenMM UMA run (LAMMPS uses 'fix langevin zero yes')",
    )
    _SHARED_DEFAULT = object()  # sentinel for "use default path"
    ap.add_argument(
        "--shared-initial-structure",
        type=Path,
        nargs="?",
        default=None,
        const=_SHARED_DEFAULT,
        metavar="PATH",
        help="Use same initial structure for both UMA runs: run OpenMM UMA minimize once, write LAMMPS data to PATH, then both UMA MD from that config (minimization skipped). If PATH is omitted, uses <out-dir>/shared_data.ice_uma.",
    )
    args = ap.parse_args()

    if args.nx < 1 or args.ny < 1 or args.nz < 1:
        ap.error("--nx, --ny, and --nz must be >= 1")

    work_dir = args.work_dir.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.shared_initial_structure is _SHARED_DEFAULT:
        args.shared_initial_structure = out_dir / "shared_data.ice_uma"

    dt_fs = args.dt_fs
    steps = args.steps if args.steps is not None else int(1000.0 / dt_fs)  # 1 ps
    output_every = max(1, int(10.0 / dt_fs))  # every 10 fs

    if args.only_plot:
        csv_c = args.order_csv_classical or (out_dir / "ice_order_classical.csv")
        csv_o = args.order_csv_uma_openmm or (out_dir / "ice_order_uma_openmm.csv")
        csv_l = args.order_csv_uma_lammps or (out_dir / "ice_order_uma_lammps.csv")
        plot_order_vs_time(csv_c, csv_o, csv_l, out_dir / "order_vs_time.png")
        return

    # Step 0–1: proton-ordered ice Ih supercell (same builder as RPMD benchmark)
    data_path = _build_data_supercell(work_dir, args.nx, args.ny, args.nz)
    n_molecules_effective = 8 * args.nx * args.ny * args.nz

    shared_data_path: Path | None = None
    skip_minimize = False
    if args.shared_initial_structure is not None:
        shared_path = args.shared_initial_structure if args.shared_initial_structure.is_absolute() else work_dir / args.shared_initial_structure
        print(f"Creating shared initial structure: OpenMM UMA minimize -> {shared_path}")
        subprocess.check_call(
            [
                sys.executable,
                str(_SCRIPT_DIR / "export_minimized_ice.py"),
                "--data", str(data_path),
                "-o", str(shared_path),
            ],
            cwd=str(work_dir),
        )
        shared_data_path = shared_path
        skip_minimize = True

    # Step 2: OpenMM classical (use smaller cutoff for small boxes)
    if not args.skip_classical:
        # Orthorhombic ice supercells can have a side < ~1 nm; default 0.7 nm cutoff then violates PBC.
        cutoff_classical = 0.4 if n_molecules_effective <= 64 else None
        _run_openmm_classical(
            work_dir,
            out_dir,
            steps,
            dt_fs,
            output_every,
            data_path=data_path,
            platform=args.platform_classical,
            variable_step=args.variable_step_classical,
            cutoff_nm=cutoff_classical,
        )
    else:
        print("Skipping TIP4P/2005f (--skip-classical).")
    # Step 3: OpenMM UMA
    _run_openmm_uma(
        work_dir, out_dir, steps, dt_fs, output_every,
        data_path=shared_data_path or data_path, skip_minimize=skip_minimize,
        temperature=args.temperature, friction=args.friction,
        integrator=args.integrator, no_cm_removal=args.no_cm_removal,
    )
    # Step 4: LAMMPS UMA (optional)
    if not args.skip_lammps:
        _run_lammps_uma(
            work_dir, out_dir, steps, output_every, dt_fs,
            data_path=shared_data_path or data_path, skip_minimize=skip_minimize,
            temperature=args.temperature, friction=args.friction,
        )
    else:
        print("Skipping LAMMPS (--skip-lammps).")
    # Step 5: plot (use whatever CSVs exist)
    csv_c = out_dir / "ice_order_classical.csv"
    csv_suffix = f"_openmm_{args.integrator}" if args.integrator != "langevin-middle" else "_openmm"
    csv_o = out_dir / f"ice_order_uma{csv_suffix.replace('-', '_')}.csv"
    csv_l = out_dir / "ice_order_uma_lammps.csv"
    if csv_l.is_file() or csv_c.is_file() or csv_o.is_file():
        plot_order_vs_time(csv_c, csv_o, csv_l, out_dir / "order_vs_time.png")
    else:
        print("No order CSVs found; skipping plot.", file=sys.stderr)
    print("Pipeline done.")


if __name__ == "__main__":
    main()
