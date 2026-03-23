#!/usr/bin/env python3
"""
Run LAMMPS + UMA on ice using FAIRChem's fix-external driver (fairchem-lammps).

Energy minimization runs *after* fix external is registered (UMA forces required).

Prerequisites:
  pip install fairchem-lammps fairchem-core lammps

Usage:
  cd tests/uma_ice_rpmd/lammps
  python build_lammps_ice_data.py -i ../ice.cif -o data.ice_uma
  python run_lammps_uma_ice.py
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_ROOT = _SCRIPT_DIR.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def _split_minimize_then_md(script: str) -> tuple[list[str], list[str]]:
    """Lines before first 'minimize' → pre (no ML forces). Rest → post (ML + min + MD)."""
    lines = script.splitlines()
    pre: list[str] = []
    post: list[str] = []
    seen_min = False
    for line in lines:
        if line.strip().startswith("#"):
            (post if seen_min else pre).append(line)
            continue
        if not seen_min and line.strip().lower().startswith("minimize"):
            seen_min = True
        if seen_min:
            post.append(line)
        else:
            pre.append(line)
    if not seen_min:
        # No minimize: behave like fairchem (only "run" after external — we still need external first)
        pre = []
        post = []
        for line in lines:
            s = line.strip()
            if s.startswith("run ") or s == "run":
                post.append(line)
            else:
                pre.append(line)
    return pre, post


def _rewrite_post_for_parity(
    post: list[str],
    steps: int | None = None,
    dump_every: int | None = None,
    thermo_every: int | None = None,
    min_eval: int | None = None,
    skip_minimize: bool = False,
) -> list[str]:
    """Replace run, thermo, dump interval, and minimize max eval in post block for OpenMM parity."""
    out: list[str] = []
    for line in post:
        s = line.strip()
        if steps is not None and (s.startswith("run ") or s == "run"):
            out.append(f"run             {steps}\n")
            continue
        if thermo_every is not None and s.startswith("thermo "):
            out.append(f"thermo          {thermo_every}\n")
            continue
        if dump_every is not None and s.startswith("dump "):
            parts = line.split()
            if len(parts) >= 4:
                # dump ID group style N file args...
                parts[4] = str(dump_every)
                out.append(" ".join(parts) + "\n")
            else:
                out.append(line)
            continue
        if skip_minimize and s.lower().startswith("minimize "):
            out.append("minimize        0.0 1.0e-6 0 0\n")
            continue
        if min_eval is not None and s.lower().startswith("minimize "):
            parts = s.split()
            if len(parts) >= 5:
                # minimize etol ftol maxiter maxeval
                parts[4] = str(min_eval)
                out.append(" ".join(parts) + "\n")
            else:
                out.append(line)
            continue
        out.append(line)
    return out


def main():
    ap = argparse.ArgumentParser(description="LAMMPS + UMA ice (fairchem-lammps)")
    ap.add_argument("--input", "-i", type=Path, default=_ROOT / "ice.cif")
    ap.add_argument("--molecules", "-n", type=int, default=None)
    ap.add_argument("--model", default="uma-s-1p1")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--data", type=Path, default=_SCRIPT_DIR / "data.ice_uma")
    ap.add_argument("--infile", type=Path, default=_SCRIPT_DIR / "in.ice_uma.lmp")
    ap.add_argument("--spin", type=int, default=1, help="omol spin (match OpenMM)")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--steps", type=int, default=None, help="Override run length (match OpenMM --steps)")
    ap.add_argument("--dump-every", type=int, default=None, help="Dump every N steps (match OpenMM --traj-every)")
    ap.add_argument("--thermo-every", type=int, default=None, help="Thermo every N steps (match OpenMM --report-every)")
    ap.add_argument("--min-eval", type=int, default=None, help="Minimize max force evaluations (match OpenMM --minimize-iters)")
    ap.add_argument(
        "--dt-fs",
        type=float,
        default=None,
        help="Timestep in fs (metal: timestep in ps = dt_fs * 0.001). If set, overrides input script timestep for pipeline parity.",
    )
    ap.add_argument(
        "--skip-minimize",
        action="store_true",
        help="Skip minimization (0 iterations). Use when starting from an already-minimized structure (e.g. --shared-initial-structure in pipeline).",
    )
    ap.add_argument(
        "--temperature",
        type=float,
        default=None,
        help="Temperature in K (overrides input script; must match OpenMM for parity). units metal: K.",
    )
    ap.add_argument(
        "--friction",
        type=float,
        default=None,
        help="Langevin friction in 1/ps (overrides input script). LAMMPS damp (ps) = 1/friction, e.g. 10 -> damp 0.1.",
    )
    ap.add_argument("--seed", type=int, default=284759, help="Random seed for velocity and Langevin (match OpenMM)")
    ap.add_argument(
        "--log-file",
        type=Path,
        default=_SCRIPT_DIR / "lammps_uma_run.log",
        help=(
            "LAMMPS -log path (thermo lines for post-processing). "
            "Default: lammps/lammps_uma_run.log. Pass with lammps_order_from_dump.py --thermo-log."
        ),
    )
    ap.add_argument(
        "--no-log-file",
        action="store_true",
        help="Use LAMMPS -log none (no thermo file; order CSV T_K/PE stay empty).",
    )
    args = ap.parse_args()

    try:
        from fairchem.lammps.lammps_fc import (
            FIX_EXT_ID,
            FIX_EXTERNAL_CMD,
            FixExternalCallback,
            check_input_script,
        )
        from lammps import lammps
    except ImportError as e:
        print("Install: pip install fairchem-lammps fairchem-core lammps", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)

    try:
        from fairchem.core import pretrained_mlip
    except ImportError:
        from fairchem.core.calculate import pretrained_mlip

    if not args.data.exists() or args.molecules is not None:
        build = _SCRIPT_DIR / "build_lammps_ice_data.py"
        cmd = [sys.executable, str(build), "-i", str(args.input), "-o", str(args.data)]
        if args.molecules:
            cmd += ["-n", str(args.molecules)]
        print(" ", " ".join(cmd))
        subprocess.check_call(cmd)

    in_text = args.infile.read_text()
    # Resolve before os.chdir(_SCRIPT_DIR): relative paths are relative to invocation cwd.
    data_path_resolved = Path(args.data).expanduser().resolve()
    if "read_data       data.ice_uma" in in_text and data_path_resolved.name != "data.ice_uma":
        # Use absolute path so LAMMPS finds the file when cwd is _SCRIPT_DIR (e.g. shared_data.ice_uma in pipeline_out).
        in_text = in_text.replace("read_data       data.ice_uma", f"read_data       {data_path_resolved}")
    if args.dt_fs is not None:
        # Override timestep for pipeline parity (same dt as OpenMM). Metal units: 1 fs = 0.001 ps.
        # Replace only the actual LAMMPS command line (not comment lines), since count=1 would
        # otherwise hit the first "timestep" in a comment and leave the real command unchanged.
        dt_ps = args.dt_fs * 0.001
        timestep_re = re.compile(r"^(\s*timestep\s+)[\d.e+-]+")
        new_lines = []
        for line in in_text.splitlines():
            stripped = line.strip()
            if stripped and not stripped.startswith("#") and timestep_re.match(line):
                line = timestep_re.sub(rf"\g<1>{dt_ps}", line)
            new_lines.append(line)
        in_text = "\n".join(new_lines) + ("\n" if in_text.endswith("\n") else "")
    # Temperature/friction/seed parity with OpenMM. Units: metal → T in K, time in ps; LAMMPS damp (ps) = 1/friction (1/ps).
    if args.temperature is not None or args.friction is not None:
        T = args.temperature if args.temperature is not None else 243.0
        damp_ps = (1.0 / args.friction) if args.friction is not None else 0.1
        seed = args.seed
        if args.temperature is not None:
            in_text = re.sub(
                r"velocity\s+all\s+create\s+[\d.]+\s+\d+(\s+dist\s+gaussian)?",
                f"velocity        all create {T} {seed}\\1",
                in_text,
                count=1,
            )
            in_text = re.sub(r"velocity\s+all\s+scale\s+[\d.]+", f"velocity        all scale {T}", in_text, count=1)
        in_text = re.sub(
            r"(fix\s+\S+\s+all\s+langevin)\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+\d+",
            f"\\1 {T} {T} {damp_ps} {seed}",
            in_text,
            count=1,
        )
    if "read_data       data.ice_uma" in in_text and data_path_resolved.name != "data.ice_uma":
        tmp_in = _SCRIPT_DIR / "_in_active.lmp"
        tmp_in.write_text(in_text)
        in_path = str(tmp_in)
    else:
        in_path = str(args.infile)

    os.chdir(_SCRIPT_DIR)
    args.data = data_path_resolved
    check_input_script(in_text)

    dump_custom = None
    dump_xyz = None
    for line in in_text.splitlines():
        s = line.strip()
        if s.startswith("dump ") and "custom" in s:
            parts = s.split()
            if len(parts) >= 6:
                dump_custom = parts[5]
        if s.startswith("dump ") and " xyz " in s:
            parts = s.split()
            if len(parts) >= 6:
                dump_xyz = parts[5]

    print(f"Model {args.model} device={args.device} task=omol spin={args.spin}")
    if dump_custom:
        print(f"  lammpstrj (id type mass x y z): {dump_custom}")
    if dump_xyz:
        print(f"  xyz (O/H for Ovito):            {dump_xyz}")

    pre, post = _split_minimize_then_md(in_text)
    if any(x is not None for x in (args.steps, args.dump_every, args.thermo_every, args.min_eval)) or args.skip_minimize:
        post = _rewrite_post_for_parity(
            post,
            steps=args.steps,
            dump_every=args.dump_every,
            thermo_every=args.thermo_every,
            min_eval=args.min_eval,
            skip_minimize=args.skip_minimize,
        )
    log_arg = "none"
    if not args.no_log_file and args.log_file is not None:
        log_arg = str(args.log_file.resolve())
        args.log_file.parent.mkdir(parents=True, exist_ok=True)

    machine = os.environ.get("LAMMPS_MACHINE_NAME")
    lmp = lammps(name=machine, cmdargs=["-nocite", "-log", log_arg, "-echo", "screen"])
    lmp._predictor = pretrained_mlip.get_predict_unit(args.model, device=args.device)
    lmp._task_name = "omol"

    lmp.commands_list([x for x in pre if x.strip()])
    lmp.command(FIX_EXTERNAL_CMD)
    lmp.set_fix_external_callback(FIX_EXT_ID, FixExternalCallback(charge=args.charge, spin=args.spin), lmp)
    print("--- UMA fix external on → minimize + MD ---")
    lmp.commands_list([x for x in post if x.strip() or x == ""])

    del lmp._predictor
    print("Done.")
    if log_arg != "none":
        print(f"  Thermo log (for lammps_order_from_dump --thermo-log): {log_arg}")
    if dump_custom and (_SCRIPT_DIR / dump_custom).is_file():
        ext = _SCRIPT_DIR / (Path(dump_custom).stem + ".extxyz")
        try:
            subprocess.run(
                [
                    sys.executable,
                    str(_SCRIPT_DIR / "lammpstrj_to_extxyz.py"),
                    dump_custom,
                    str(data_path_resolved),
                    "-o",
                    str(ext),
                ],
                cwd=str(_SCRIPT_DIR),
                check=True,
            )
            print(f"  Ovito (best): {ext}  — Lattice + O/H species")
        except Exception as e:
            print(f"  extxyz conversion skipped: {e}")
    if dump_xyz and (_SCRIPT_DIR / dump_xyz).is_file():
        print(f"  Ovito (xyz):  {_SCRIPT_DIR / dump_xyz}  — element x y z (set cell in Ovito if needed)")


if __name__ == "__main__":
    main()
