#!/usr/bin/env python3
"""
LAMMPS RPMD with UMA potential — no i-PI server required.

Uses LAMMPS's built-in fix pimd/langevin (REPLICA package) for ring-polymer
molecular dynamics, with forces from FairChem UMA via fix external pf/callback.

Architecture
------------
    N MPI ranks  →  N LAMMPS partitions  →  N ring-polymer beads
    Each partition evaluates UMA forces independently via fix external.
    fix pimd/langevin handles inter-bead springs, PILE_L thermostat, and
    velocity-Verlet integration with OBABO Trotter splitting.

    No i-PI server, no socket communication — everything runs in-process.

Requirements
------------
    conda activate lammps-pimd   # Python 3.12, mpich 5, mpi4py, pip lammps Jul-2025
    # or: conda activate lammps-env  (conda-forge lammps Aug-2024, openmpi)

Usage
-----
    mpirun -np <nbeads> python run_lammps_pimd_uma.py [options]

    The number of MPI ranks MUST equal --beads (one partition per bead).

Equivalent OpenMM commands
--------------------------
    OpenMM reference (run_openmm_rpmd_reference.py default):
        python run_openmm_rpmd_reference.py --nx 2 --ny 2 --nz 2 --beads 32 \\
            --dt 0.5 --prod 1.0 --equil 2.0 --rpmd-thermostat pile-g \\
            --rpmd-friction 1.0 --rpmd-centroid-friction 0.5

    Equivalent LAMMPS PIMD (lammps-pimd env):
        mpirun -np 32 python run_lammps_pimd_uma.py --beads 32 \\
            --dt-fs 0.5 --prod 1.0 --equil 2.0 \\
            --tau 2.0 --method nmpimd --data ../lammps/data.ice_uma_64

    Parameter mapping:
        OpenMM  --dt 0.5 fs              →  LAMMPS  --dt-fs 0.5
        OpenMM  --prod 1.0 ps            →  LAMMPS  --prod 1.0
        OpenMM  --equil 2.0 ps           →  LAMMPS  --equil 2.0
        OpenMM  --beads 32               →  LAMMPS  --beads 32  (+  mpirun -np 32)
        OpenMM  --temp 243               →  LAMMPS  --temp 243
        OpenMM  --rpmd-friction 1.0 /ps  →  LAMMPS  --tau 1.0 ps  (tau = 1/friction)
        OpenMM  pile-g centroid friction →  LAMMPS  --tau <1/centroid_friction> ps
        OpenMM  --seed 284759            →  LAMMPS  --seed 284759

    Thermostat note:
        OpenMM pile-g = Bussi/SVR stochastic velocity rescaling on centroid +
            Langevin (PILE) on internal modes.
        LAMMPS PILE_L = Langevin on ALL modes including centroid.
        These produce the same canonical ensemble but different centroid dynamics.
        For thermodynamic properties (kinetic energy estimators, density) results
        should agree; for centroid diffusion / time-correlation functions they differ.
        PILE_L is the standard choice for structural/thermodynamic PIMD sampling.

GPU Memory
----------
    Each MPI rank loads its own UMA predictor (~2 GB on GPU for uma-s + activations).
    With --beads 4:  ~8 GB.  With --beads 32: need --multi-gpu across 4+ GPUs.
    Single 12 GB GPU: up to ~4 beads (with no other GPU processes running).

Thermo Output (fix pimd/langevin global vector, NVT)
-----------------------------------------------------
    f_rpmd[1]  = kinetic energy of this bead/mode
    f_rpmd[2]  = spring elastic energy of this bead/mode
    f_rpmd[3]  = potential energy of this bead
    f_rpmd[4]  = total energy of all beads (conserved in NVE)
    f_rpmd[5]  = primitive kinetic energy estimator  (quantum KE)
    f_rpmd[6]  = virial energy estimator             (quantum KE, better convergence)
    f_rpmd[7]  = centroid-virial energy estimator     (quantum KE)
    f_rpmd[8]  = primitive pressure estimator
    f_rpmd[9]  = thermodynamic pressure estimator
    f_rpmd[10] = centroid-virial pressure estimator

    Items [1]-[3] differ per partition log; [4]-[10] are the same across partitions.
    Best quantum KE estimator: f_rpmd[7] (centroid-virial, lowest variance).

References
----------
    LAMMPS fix pimd/langevin:
        https://docs.lammps.org/fix_pimd.html
    PILE_L thermostat:
        Ceriotti, Parrinello, Markland, Manolopoulos, JCP 133, 124104 (2010)
    OpenMM PILE-G vs PILE_L comparison:
        Ceriotti, More, Manolopoulos, CPC 185, 1019 (2014)
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent

# Minimal stub: LAMMPS requires -in with -partition, but we issue all
# commands from Python after the constructor.  A stub file satisfies the
# LAMMPS requirement without defining any simulation state.
_STUB_CONTENTS = "# stub for -partition (actual setup via Python API)\n"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="LAMMPS PIMD/RPMD + UMA (no i-PI server)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--beads", type=int, required=True,
        help="Number of ring-polymer beads (must equal MPI ranks)",
    )
    ap.add_argument("--temp", type=float, default=243.0, help="Temperature in K (default: 243)")
    ap.add_argument("--steps", type=int, default=None, help="MD steps (overrides --prod)")
    ap.add_argument("--prod", type=float, default=None, help="Production time in ps")
    ap.add_argument(
        "--equil", type=float, default=0.0,
        help="Equilibration time in ps before production, no output (default: 0). "
             "Match OpenMM --equil 2.0 for comparable equilibrated runs.",
    )
    ap.add_argument("--dt-fs", type=float, default=0.5, help="Timestep in fs (default: 0.5)")
    ap.add_argument(
        "--data", type=Path, default=None,
        help="LAMMPS data file (default: data.ice_uma_<molecules>)",
    )
    ap.add_argument("--molecules", type=int, default=64, help="Water molecules for data file naming (default: 64)")
    ap.add_argument("--model", default="uma-s-1p1", help="FairChem model name")
    ap.add_argument("--device", default="cuda", help="PyTorch device (default: cuda)")
    ap.add_argument("--multi-gpu", action="store_true", help="Distribute beads across multiple GPUs")
    ap.add_argument("--ngpus", type=int, default=None, help="Number of GPUs (with --multi-gpu)")
    ap.add_argument("--task-name", default="omol", help="UMA task name")
    ap.add_argument("--seed", type=int, default=284759, help="Base random seed")
    ap.add_argument(
        "--tau", type=float, default=1.0,
        help="PILE_L thermostat damping time for centroid mode in ps (default: 1.0)",
    )
    ap.add_argument(
        "--scale", type=float, default=1.0,
        help="Scaling factor for non-centroid mode damping (default: 1.0)",
    )
    ap.add_argument(
        "--method", default="pimd", choices=["pimd", "nmpimd"],
        help="PIMD integration: pimd (Cartesian, default) or nmpimd (normal-mode)",
    )
    ap.add_argument(
        "--integrator", default="obabo", choices=["obabo", "baoab"],
        help="Trotter splitting scheme (default: obabo)",
    )
    ap.add_argument("--thermo-every", type=int, default=10, help="Thermo output interval (default: 10)")
    ap.add_argument("--dump-every", type=int, default=100, help="Trajectory dump interval (default: 100)")
    ap.add_argument("--outdir", type=Path, default=None, help="Output directory (default: pimd_output/)")
    ap.add_argument("--charge", type=int, default=0, help="System charge")
    ap.add_argument("--spin", type=int, default=1, help="System spin (UMA omol default: 1)")
    args = ap.parse_args()

    # ── MPI setup ──────────────────────────────────────────────────────────
    try:
        from mpi4py import MPI
    except ImportError:
        print(
            "ERROR: mpi4py is required.\n"
            "  pip install mpi4py\n"
            "  Then: mpirun -np <nbeads> python run_lammps_pimd_uma.py ...",
            file=sys.stderr,
        )
        sys.exit(1)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()

    if nranks != args.beads:
        if rank == 0:
            print(
                f"ERROR: MPI ranks ({nranks}) != --beads ({args.beads}).\n"
                f"  Run: mpirun -np {args.beads} python {Path(sys.argv[0]).name} --beads {args.beads} ...",
                file=sys.stderr,
            )
        comm.Abort(1)

    # ── GPU assignment ─────────────────────────────────────────────────────
    if args.multi_gpu:
        import torch
        ngpus = args.ngpus or torch.cuda.device_count()
        device = f"cuda:{rank % ngpus}"
        if rank == 0:
            print(f"Multi-GPU: {ngpus} GPUs, assigning beads round-robin")
    else:
        device = args.device

    os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

    # ── Import LAMMPS and FairChem ─────────────────────────────────────────
    try:
        from lammps import lammps
    except ImportError:
        if rank == 0:
            print("ERROR: LAMMPS Python interface required. pip install lammps", file=sys.stderr)
        comm.Abort(1)

    from fairchem.lammps.lammps_fc import (
        FIX_EXT_ID,
        FIX_EXTERNAL_CMD,
        FixExternalCallback,
    )
    try:
        from fairchem.core import pretrained_mlip
    except ImportError:
        from fairchem.core.calculate import pretrained_mlip

    # ── Resolve paths ──────────────────────────────────────────────────────
    data_file = args.data
    if data_file is None:
        candidates = [
            _SCRIPT_DIR / f"data.ice_uma_{args.molecules}",
            _SCRIPT_DIR / "data.ice_uma",
        ]
        for c in candidates:
            if c.is_file():
                data_file = c
                break
    if data_file is None or not data_file.is_file():
        if rank == 0:
            searched = [str(c) for c in (candidates if args.data is None else [args.data])]
            print(
                f"ERROR: Data file not found. Searched: {searched}\n"
                f"  Build one: python build_lammps_ice_data.py --nx 2 --ny 2 --nz 2 -o data.ice_uma_64",
                file=sys.stderr,
            )
        comm.Abort(1)
    data_file = data_file.resolve()

    outdir = (args.outdir or (_SCRIPT_DIR / "pimd_output")).resolve()
    if rank == 0:
        outdir.mkdir(parents=True, exist_ok=True)
    comm.Barrier()

    # ── Resolve steps ──────────────────────────────────────────────────────
    if args.steps is not None:
        steps = args.steps
    elif args.prod is not None:
        steps = int(args.prod * 1000.0 / args.dt_fs)
    else:
        steps = 1000

    equil_steps = int(args.equil * 1000.0 / args.dt_fs) if args.equil > 0 else 0
    dt_ps = args.dt_fs * 0.001  # metal units: ps

    if rank == 0:
        print("=" * 60)
        print("LAMMPS PIMD + UMA  (no i-PI)")
        print("=" * 60)
        print(f"  Beads:        {args.beads}")
        print(f"  Temperature:  {args.temp} K")
        print(f"  Timestep:     {args.dt_fs} fs  ({dt_ps:.6f} ps)")
        print(f"  Steps:        {steps}  ({steps * args.dt_fs / 1000:.3f} ps)")
        if equil_steps:
            print(f"  Equil:        {equil_steps} steps ({args.equil:.3f} ps, no output)")
        print(f"  Method:       {args.method}")
        print(f"  Integrator:   {args.integrator}")
        print(f"  Thermostat:   PILE_L  tau={args.tau} ps  scale={args.scale}")
        print(f"  Model:        {args.model}")
        print(f"  Device:       {device}")
        print(f"  Data:         {data_file}")
        print(f"  Output:       {outdir}")
        print("=" * 60)
        sys.stdout.flush()

    # ── GPU memory check ──────────────────────────────────────────────────
    if "cuda" in device:
        import torch
        if torch.cuda.is_available():
            gpu_idx = torch.cuda.current_device() if ":" not in device else int(device.split(":")[1])
            total_mem = torch.cuda.get_device_properties(gpu_idx).total_memory / 1e9
            free_mem_approx = total_mem * 0.9  # rough estimate after driver overhead
            per_bead_est = 2.0  # ~2 GB per UMA-s instance (model + activations)
            needed = per_bead_est * (args.beads if not args.multi_gpu else 1)
            if rank == 0:
                print(f"  GPU {gpu_idx}: {total_mem:.1f} GB total, ~{per_bead_est:.0f} GB per bead")
                if needed > free_mem_approx:
                    print(
                        f"  WARNING: {args.beads} beads need ~{needed:.0f} GB but GPU has ~{free_mem_approx:.0f} GB.\n"
                        f"    Options: fewer beads, --multi-gpu, --device cpu, or free GPU memory.\n"
                        f"    Check: nvidia-smi --query-compute-apps=pid,used_memory --format=csv"
                    )

    # ── Load UMA predictor (each rank loads its own copy) ──────────────────
    if rank == 0:
        print(f"Loading UMA predictor ({args.model}) on {nranks} ranks...")
        sys.stdout.flush()

    predictor = pretrained_mlip.get_predict_unit(args.model, device=device)
    comm.Barrier()

    if rank == 0:
        print("All ranks loaded UMA predictor.")
        sys.stdout.flush()

    # ── Create LAMMPS with partitions ──────────────────────────────────────
    # LAMMPS requires -in with -partition.  We use a stub file (comment only)
    # so that the constructor completes with a clean state.  All simulation
    # setup is then issued from Python, which avoids the state-clearing that
    # -in causes when it processes a real input file.
    stub_path = outdir / "_pimd_stub.lmp"
    if rank == 0:
        stub_path.write_text(_STUB_CONTENTS)
    comm.Barrier()

    log_file = str(outdir / f"log.pimd.{rank:03d}")
    machine = os.environ.get("LAMMPS_MACHINE_NAME")
    lmp = lammps(
        name=machine,
        cmdargs=[
            "-partition", f"{args.beads}x1",
            "-in", str(stub_path),
            "-nocite",
            "-log", log_file,
            "-echo", "screen",
        ],
    )
    lmp._predictor = predictor
    lmp._task_name = args.task_name

    # ── System setup (executed on each partition independently) ─────────────
    lmp.commands_list([
        "units metal",
        "atom_style atomic",
        "atom_modify map yes",
        "boundary p p p",
        f"read_data {data_file}",
    ])

    # UMA forces via fix external (no pair_style needed)
    lmp.command(FIX_EXTERNAL_CMD)
    callback = FixExternalCallback(charge=args.charge, spin=args.spin)
    lmp.set_fix_external_callback(FIX_EXT_ID, callback, lmp)

    # Each bead gets unique velocity seed (LAMMPS PIMD docs recommendation)
    vseed = args.seed + (rank + 1) * 13579
    lmp.command(f"velocity all create {args.temp} {vseed} rot yes dist gaussian")

    # ── RPMD integration: fix pimd/langevin ────────────────────────────────
    pimd_cmd = (
        f"fix rpmd all pimd/langevin "
        f"method {args.method} "
        f"integrator {args.integrator} "
        f"ensemble nvt "
        f"temp {args.temp} "
        f"thermostat PILE_L {args.seed} "
        f"tau {args.tau}"
    )
    if args.scale != 1.0:
        pimd_cmd += f" scale {args.scale}"
    lmp.command(pimd_cmd)

    lmp.command(f"timestep {dt_ps}")

    # ── Output ─────────────────────────────────────────────────────────────
    lmp.command(
        "thermo_style custom step time temp pe ke "
        "f_rpmd[4] f_rpmd[5] f_rpmd[6] f_rpmd[7]"
    )
    lmp.command(f"thermo {args.thermo_every}")

    # Per-bead trajectory dump.  Positions are Cartesian for all methods;
    # velocities are in normal-mode space for method=nmpimd.
    bead_tag = f"{rank:03d}"
    dump_file = str(outdir / f"dump.bead_{bead_tag}.lammpstrj")
    lmp.command(
        f"dump traj all custom {args.dump_every} {dump_file} "
        f"id type x y z vx vy vz ix iy iz"
    )
    lmp.command("dump_modify traj sort id flush yes")

    # ── Equilibration (no dump/thermo output) ─────────────────────────────
    if equil_steps > 0:
        if rank == 0:
            print(f"\nEquilibrating {equil_steps} steps ({args.equil:.3f} ps), no output...")
            sys.stdout.flush()
        lmp.command("thermo 0")
        lmp.command(f"run {equil_steps}")
        lmp.command(f"thermo {args.thermo_every}")

    # ── Production run ─────────────────────────────────────────────────────
    if rank == 0:
        print(f"\nRunning {steps} steps ({steps * args.dt_fs / 1000:.3f} ps)...")
        sys.stdout.flush()

    lmp.command(f"run {steps}")

    # ── Cleanup ────────────────────────────────────────────────────────────
    if rank == 0:
        print("\nDone.")
        print(f"  Logs:    {outdir}/log.pimd.XXX  (one per bead)")
        print(f"  Dumps:   {outdir}/dump.bead_XXX.lammpstrj")
        print(f"  Quantum estimators in thermo: f_rpmd[5] (prim KE), f_rpmd[6] (virial KE)")
        sys.stdout.flush()

    del lmp._predictor
    lmp.close()


if __name__ == "__main__":
    main()
