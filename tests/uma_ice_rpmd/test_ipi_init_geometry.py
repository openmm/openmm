"""Tests for i-PI init geometry: symmetric cell (fix ipi parity) and diagnostics."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_SCRIPT_DIR = Path(__file__).resolve().parent
_CONVERTER = _SCRIPT_DIR / "ipi" / "convert_lammps_to_ipi_xyz.py"
_DIAGNOSE = _SCRIPT_DIR / "diagnose_fix_ipi_box.py"
_IPI_DIR = _SCRIPT_DIR / "ipi"

if str(_IPI_DIR) not in sys.path:
    sys.path.insert(0, str(_IPI_DIR))


def test_shift_orthorhombic_centers_to_origin() -> None:
    import convert_lammps_to_ipi_xyz as conv

    shift_orthorhombic_positions_to_symmetric_cell = conv.shift_orthorhombic_positions_to_symmetric_cell

    pos = np.array([[0.0, 0.0, 0.0], [8.994, 15.578, 14.644]], dtype=np.float64)
    out = shift_orthorhombic_positions_to_symmetric_cell(
        pos,
        xlo=0.0,
        xhi=8.994,
        ylo=0.0,
        yhi=15.578,
        zlo=0.0,
        zhi=14.644,
        xy=0.0,
        xz=0.0,
        yz=0.0,
    )
    assert out[0, 0] == pytest.approx(-4.497)
    assert out[1, 0] == pytest.approx(4.497)
    assert out[0, 1] == pytest.approx(-7.789)
    assert out[1, 2] == pytest.approx(7.322)


def test_shift_skips_triclinic_tilt() -> None:
    import convert_lammps_to_ipi_xyz as conv

    shift_orthorhombic_positions_to_symmetric_cell = conv.shift_orthorhombic_positions_to_symmetric_cell

    pos = np.array([[1.0, 2.0, 3.0]], dtype=np.float64)
    out = shift_orthorhombic_positions_to_symmetric_cell(
        pos,
        xlo=0.0,
        xhi=10.0,
        ylo=0.0,
        yhi=10.0,
        zlo=0.0,
        zhi=10.0,
        xy=0.5,
        xz=0.0,
        yz=0.0,
    )
    assert np.allclose(out, pos)


def test_converter_output_bounds_for_ice_data(tmp_path: Path) -> None:
    """Orthorhombic UMA ice data -> init coords inside [-L/2, L/2] per axis."""
    data = _SCRIPT_DIR / "lammps" / "data.ice_uma_64"
    if not data.is_file():
        pytest.skip("lammps/data.ice_uma_64 not in tree")
    out_xyz = tmp_path / "init.xyz"
    subprocess.run(
        [sys.executable, str(_CONVERTER), "--data", str(data), "-o", str(out_xyz)],
        check=True,
        cwd=str(_SCRIPT_DIR),
    )
    lines = out_xyz.read_text().splitlines()
    nat = int(lines[0])
    lx, ly, lz = 8.994, 15.578, 14.644  # from data file (Å)
    half = np.array([lx / 2, ly / 2, lz / 2])
    for ln in lines[2 : 2 + nat]:
        parts = ln.split()
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        assert np.all(np.abs(np.array([x, y, z])) <= half + 1e-6)


@pytest.mark.skipif(
    not _DIAGNOSE.is_file(),
    reason="diagnose_fix_ipi_box.py missing",
)
def test_diagnose_fix_ipi_box_reports_invariant_energy() -> None:
    """Integration: UMA energy unchanged by fix_ipi-style box move (needs fairchem env)."""
    try:
        proc = subprocess.run(
            [
                sys.executable,
                str(_DIAGNOSE),
                "--quick",
                "--task-name",
                "omol",
            ],
            cwd=str(_SCRIPT_DIR),
            capture_output=True,
            text=True,
            timeout=180,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pytest.skip("diagnostic subprocess failed to start or timed out")
    if proc.returncode != 0:
        # Missing CUDA / broken torch / no lammps — skip rather than fail CI
        combined = proc.stdout + proc.stderr
        if any(
            x in combined
            for x in (
                "Requires fairchem",
                "ImportError",
                "undefined symbol",
                "No module named",
            )
        ):
            pytest.skip("UMA/LAMMPS stack not available: " + combined[:200])
        pytest.fail(combined)

    text = proc.stdout + proc.stderr
    m = re.search(r"ASCII_DELTA_etotal_eV = ([\d.eE+-]+)", text)
    assert m is not None, text[-2000:]
    delta = float(m.group(1))
    assert abs(delta) < 1e-2, f"Box move should not change UMA energy materially; got Δ={delta} eV"
