#!/usr/bin/env python3
"""
Orchestrate i-PI + LAMMPS UMA for RPMD benchmark vs OpenMM.

Starts i-PI (PIMD/RPMD driver), waits for socket, then runs LAMMPS as force
client with fix ipi + fix external (UMA). Parameters match OpenMM RPMD:
243 K, PILE-L, 0.5 fs timestep. Default 32 molecules, 4 beads (scaled from 128×8).

Prerequisites:
  pip install ipi fairchem-lammps fairchem-core lammps
  LAMMPS built with MISC package (fix ipi)

Usage:
  cd tests/uma_ice_rpmd
  python run_ipi_lammps_uma_rpmd.py [--steps 1000] [--device cuda]
  python run_ipi_lammps_uma_rpmd.py --molecules 32 --beads 4 --steps 100   # short test
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_ROOT = _SCRIPT_DIR
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

IPI_DIR = _SCRIPT_DIR / "ipi"
LAMMPS_DIR = _SCRIPT_DIR / "lammps"
IPI_PORT_DEFAULT = 65535


def _wait_for_port(host: str, port: int, timeout_s: float = 60.0, interval_s: float = 0.2) -> bool:
    """Return True when port is accepting connections."""
    import socket
    t0 = time.time()
    while time.time() - t0 < timeout_s:
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.settimeout(1.0)
                if s.connect_ex((host, port)) == 0:
                    return True
        except OSError:
            pass
        time.sleep(interval_s)
    return False


def _write_ipi_input(
    output_path: Path,
    molecules: int,
    beads: int,
    total_steps: int,
    *,
    dt_fs: float = 0.5,
    seed: int = 284759,
    port: int = IPI_PORT_DEFAULT,
) -> None:
    """Write i-PI input.xml with natoms=3*molecules, nbeads=beads."""
    natoms = 3 * molecules
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<simulation verbosity="medium">
  <!-- RPMD benchmark: {molecules} mol, {beads} beads, 243 K, PILE-L, {dt_fs} fs -->
  <total_steps>{total_steps}</total_steps>
  <prng>
    <seed>{seed}</seed>
  </prng>
  <output>
    <trajectory filename="traj" stride="10" format="xyz">positions</trajectory>
    <checkpoint filename="RESTART" stride="{total_steps}"/>
  </output>

  <ffsocket mode="inet" name="lammps">
    <address>127.0.0.1</address>
    <port>{port}</port>
    <timeout>3600</timeout>
  </ffsocket>

  <system prefix="ice_">
    <forces>
      <force forcefield="lammps"/>
    </forces>
    <ensemble>
      <temperature units="kelvin">243</temperature>
    </ensemble>
    <motion mode="dynamics">
      <dynamics mode="nvt">
        <timestep units="femtosecond">{dt_fs}</timestep>
        <thermostat mode="pile_l">
          <tau units="femtosecond">1000</tau>
        </thermostat>
      </dynamics>
    </motion>
    <beads nbeads="{beads}" natoms="{natoms}"/>
    <initialize nbeads="{beads}">
      <file mode="xyz" units="angstrom">init.xyz</file>
    </initialize>
  </system>
</simulation>
"""
    output_path.write_text(xml)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Run i-PI + LAMMPS UMA RPMD benchmark"
    )
    ap.add_argument("--steps", type=int, default=None, help="Override: MD steps (default: from --prod and --dt-fs)")
    ap.add_argument("--prod", type=float, default=None, help="Production time in ps (e.g. 100); used with --dt-fs to set steps")
    ap.add_argument("--dt-fs", type=float, default=0.5, help="Timestep in fs (default 0.5)")
    ap.add_argument("--seed", type=int, default=284759, help="i-PI PRNG seed")
    ap.add_argument(
        "--port",
        type=int,
        default=None,
        help=f"i-PI socket port (default {IPI_PORT_DEFAULT}, use for SLURM array isolation)",
    )
    ap.add_argument("--molecules", type=int, default=32, help="Water molecules (default 32, scaled from 128)")
    ap.add_argument("--beads", type=int, default=4, help="RPMD beads (default 4, scaled from 8)")
    ap.add_argument("--device", default="cuda", help="Device for UMA (cuda/cpu)")
    ap.add_argument("--model", default="uma-s-1p1")
    ap.add_argument(
        "--data",
        type=Path,
        default=None,
        help="LAMMPS data file (default: lammps/data.ice_uma_{molecules})",
    )
    ap.add_argument("--ipi-input", type=Path, default=IPI_DIR / "input.xml")
    ap.add_argument("--lammps-in", type=Path, default=LAMMPS_DIR / "in.ice_uma_ipi_client.lmp")
    args = ap.parse_args()

    # Resolve data path: use scaled default when --data not given
    data_path = args.data or (LAMMPS_DIR / f"data.ice_uma_{args.molecules}")
    if not data_path.is_file():
        print(
            f"Missing data file: {data_path}. Run: python lammps/build_lammps_ice_data.py -n {args.molecules} -o {data_path}",
            file=sys.stderr,
        )
        sys.exit(1)
    if not args.lammps_in.is_file():
        print(f"Missing LAMMPS input: {args.lammps_in}", file=sys.stderr)
        sys.exit(1)

    # Convert LAMMPS data → i-PI init.xyz
    init_xyz = IPI_DIR / "init.xyz"
    converter = _SCRIPT_DIR / "ipi" / "convert_lammps_to_ipi_xyz.py"
    if not converter.is_file():
        print(f"Missing converter: {converter}", file=sys.stderr)
        sys.exit(1)
    subprocess.check_call(
        [sys.executable, str(converter), "--data", str(data_path), "-o", str(init_xyz)],
        cwd=str(_SCRIPT_DIR),
    )
    print(f"Wrote {init_xyz} ({3 * args.molecules} atoms)")

    # Clean stale i-PI output files from previous runs (different molecule count etc.)
    for stale in IPI_DIR.glob("ice_*i-pi*"):
        stale.unlink()
        print(f"  Removed stale: {stale.name}")
    for stale in IPI_DIR.glob("RESTART*"):
        stale.unlink()
        print(f"  Removed stale: {stale.name}")

    port = args.port if args.port is not None else IPI_PORT_DEFAULT

    # Resolve steps: from --prod/--dt-fs or explicit --steps
    if args.steps is not None:
        steps = args.steps
    elif args.prod is not None:
        steps = int(args.prod * 1000.0 / args.dt_fs)
    else:
        steps = 1000  # default short test

    # Write i-PI input.xml with natoms/nbeads
    _write_ipi_input(
        args.ipi_input,
        molecules=args.molecules,
        beads=args.beads,
        total_steps=steps,
        dt_fs=args.dt_fs,
        seed=args.seed,
        port=port,
    )
    print(f"Wrote {args.ipi_input} ({args.molecules} mol, {args.beads} beads)")

    # Start i-PI (use i-pi executable; python -m ipi may not work)
    ipi_proc = subprocess.Popen(
        ["i-pi", str(args.ipi_input)],
        cwd=str(IPI_DIR),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    print(f"Started i-PI (pid {ipi_proc.pid}) from {IPI_DIR}")

    try:
        if not _wait_for_port("127.0.0.1", port, timeout_s=60):
            print("i-PI port did not become ready in time", file=sys.stderr)
            ipi_proc.terminate()
            sys.exit(1)
        print(f"i-PI ready: 127.0.0.1:{port}")
        time.sleep(5.0)  # Give i-PI time to finish init and enter accept loop

        # Run LAMMPS client with fix external + fix ipi
        try:
            from fairchem.lammps.lammps_fc import (
                FIX_EXT_ID,
                FIX_EXTERNAL_CMD,
                FixExternalCallback,
            )
            from fairchem.core import pretrained_mlip
        except ImportError:
            from fairchem.core.calculate import pretrained_mlip

        from lammps import lammps

        predictor = pretrained_mlip.get_predict_unit(args.model, device=args.device)
        lmp = lammps(cmdargs=["-nocite", "-log", "none", "-echo", "screen"])
        lmp._predictor = predictor
        lmp._task_name = "omol"

        in_text = args.lammps_in.read_text()
        # Resolve read_data path
        data_resolved = data_path.resolve()
        if "read_data       data.ice_uma" in in_text and data_resolved.name != "data.ice_uma":
            in_text = in_text.replace(
                "read_data       data.ice_uma",
                f"read_data       {data_resolved}",
            )

        # Parse commands: skip comments, inject fix external before fix ipi
        commands_pre: list[str] = []
        commands_post: list[str] = []
        in_post = False
        for line in in_text.splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if "fix" in s and "ipi" in s:
                in_post = True
                # Replace port (last token) with configured port for SLURM isolation
                parts = s.split()
                if len(parts) >= 5:
                    parts[-1] = str(port)
                    s = " ".join(parts)
                commands_post.append(s)
                continue
            if s.startswith("run "):
                commands_post.append(f"run             {steps}")
                continue
            if in_post:
                commands_post.append(s)
            else:
                commands_pre.append(s)

        os.chdir(LAMMPS_DIR)
        for cmd in commands_pre:
            lmp.command(cmd)
        lmp.command(FIX_EXTERNAL_CMD)
        lmp.set_fix_external_callback(
            FIX_EXT_ID,
            FixExternalCallback(charge=0, spin=1),
            lmp,
        )
        for cmd in commands_post:
            lmp.command(cmd)

        del lmp._predictor
        del lmp
        print("LAMMPS client finished")

    finally:
        ipi_proc.terminate()
        ipi_proc.wait(timeout=10)
        print("i-PI stopped")

    # Report outputs
    traj = IPI_DIR / "ice_traj.xyz"
    props = IPI_DIR / "ice_out"
    if traj.exists():
        print(f"Trajectory: {traj}")
    if props.exists():
        print(f"Properties: {props}")


if __name__ == "__main__":
    main()
