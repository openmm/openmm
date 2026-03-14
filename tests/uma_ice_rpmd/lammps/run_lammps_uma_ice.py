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
    if "read_data       data.ice_uma" in in_text and args.data.name != "data.ice_uma":
        in_text = in_text.replace("read_data       data.ice_uma", f"read_data       {args.data.name}")
        tmp_in = _SCRIPT_DIR / "_in_active.lmp"
        tmp_in.write_text(in_text)
        in_path = str(tmp_in)
    else:
        in_path = str(args.infile)

    os.chdir(_SCRIPT_DIR)
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
    machine = os.environ.get("LAMMPS_MACHINE_NAME")
    lmp = lammps(name=machine, cmdargs=["-nocite", "-log", "none", "-echo", "screen"])
    lmp._predictor = pretrained_mlip.get_predict_unit(args.model, device=args.device)
    lmp._task_name = "omol"

    lmp.commands_list([x for x in pre if x.strip()])
    lmp.command(FIX_EXTERNAL_CMD)
    lmp.set_fix_external_callback(FIX_EXT_ID, FixExternalCallback(charge=args.charge, spin=args.spin), lmp)
    print("--- UMA fix external on → minimize + MD ---")
    lmp.commands_list([x for x in post if x.strip() or x == ""])

    del lmp._predictor
    print("Done.")
    data_for_box = args.data.name if args.data.name == "data.ice_uma" else str(args.data)
    if dump_custom and (_SCRIPT_DIR / dump_custom).is_file():
        ext = _SCRIPT_DIR / (Path(dump_custom).stem + ".extxyz")
        try:
            subprocess.run(
                [sys.executable, str(_SCRIPT_DIR / "lammpstrj_to_extxyz.py"), dump_custom, data_for_box, "-o", str(ext)],
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
