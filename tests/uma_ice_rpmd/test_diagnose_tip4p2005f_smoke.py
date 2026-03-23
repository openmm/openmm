"""Smoke test: diagnostic script runs without error on repo LAMMPS data."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

_SCRIPT_DIR = Path(__file__).resolve().parent


@pytest.mark.timeout(120)
def test_diagnose_tip4p2005f_smoke() -> None:
    data = _SCRIPT_DIR / "lammps" / "data.ice_uma_32"
    if not data.is_file():
        pytest.skip(f"Missing {data}")
    r = subprocess.run(
        [
            sys.executable,
            str(_SCRIPT_DIR / "diagnose_tip4p2005f.py"),
            "--skip-scaling",
            "--skip-tip4pew",
            "--data",
            str(data),
        ],
        cwd=str(_SCRIPT_DIR),
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr + r.stdout
    assert "MISSING intramolecular exclusions" not in r.stdout
