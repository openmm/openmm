"""Smoke test: diagnostic script runs without error on repo LAMMPS data."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

_SCRIPT_DIR = Path(__file__).resolve().parent


@pytest.mark.timeout(120)
def test_diagnose_tip4p2005f_smoke() -> None:
    data: Path | None = None
    for name in ("data.ice_uma_64", "data.ice_uma_32"):
        candidate = _SCRIPT_DIR / "lammps" / name
        if candidate.is_file():
            data = candidate
            break
    if data is None:
        pytest.skip("Missing lammps/data.ice_uma_64 or data.ice_uma_32")
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
    out = (r.stderr or "") + (r.stdout or "")
    if r.returncode != 0 and (
        "symbol lookup error" in out or "undefined symbol" in out
    ):
        pytest.skip(
            "Mixed conda Python with in-tree OpenMM CUDA plugin (use one consistent OpenMM build)."
        )
    assert r.returncode == 0, out
    assert "MISSING intramolecular exclusions" not in r.stdout
