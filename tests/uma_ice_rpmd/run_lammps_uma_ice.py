#!/usr/bin/env python3
"""Launcher: run from tests/uma_ice_rpmd/ — real script is lammps/run_lammps_uma_ice.py"""
import subprocess
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_script = _here / "lammps" / "run_lammps_uma_ice.py"
if not _script.is_file():
    sys.exit(f"Missing {_script} — use: cd lammps && python run_lammps_uma_ice.py")
subprocess.check_call([sys.executable, str(_script)] + sys.argv[1:], cwd=str(_here / "lammps"))
