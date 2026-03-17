"""Tests for lammps_order_from_dump (LAMMPS dump → ice order CSV)."""
from __future__ import annotations

import csv
import tempfile
from pathlib import Path

import pytest

# Import after path setup
import sys
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))


def _minimal_data_file(path: Path, lx: float = 12.0, ly: float = 12.0, lz: float = 12.0) -> None:
    """Write minimal LAMMPS data with 24 atoms (8 water), orthorhombic box so O have 4+ neighbors."""
    atoms = []
    aid = 1
    for m in range(8):
        ox = 1.5 + (m % 2) * 3.0
        oy = 1.5 + (m // 2) * 3.0
        oz = 1.5 + (m // 4) * 3.0
        atoms.append(f"{aid} 1 {ox:.2f} {oy:.2f} {oz:.2f}")
        aid += 1
        atoms.append(f"{aid} 2 {ox+0.76:.2f} {oy:.2f} {oz:.2f}")
        aid += 1
        atoms.append(f"{aid} 2 {ox:.2f} {oy+0.76:.2f} {oz:.2f}")
        aid += 1
    path.write_text(f"""LAMMPS data
24 atoms
2 atom types
0.0 {lx} xlo xhi
0.0 {ly} ylo yhi
0.0 {lz} zlo zhi

Atoms

""" + "\n".join(atoms) + "\n")


def _minimal_lammpstrj(path: Path, n_frames: int = 2) -> None:
    """Write minimal dump with id type mass x y z, 24 atoms (8 water), ortho box."""
    lines = []
    for step in [0, 10][:n_frames]:
        lines.append("ITEM: TIMESTEP")
        lines.append(str(step))
        lines.append("ITEM: NUMBER OF ATOMS")
        lines.append("24")
        lines.append("ITEM: BOX BOUNDS pp pp pp")
        lines.append("0.0 12.0")
        lines.append("0.0 12.0")
        lines.append("0.0 12.0")
        lines.append("ITEM: ATOMS id type mass x y z")
        aid = 1
        for m in range(8):
            ox = 1.5 + (m % 2) * 3.0
            oy = 1.5 + (m // 2) * 3.0
            oz = 1.5 + (m // 4) * 3.0
            lines.append(f"{aid} 1 15.999 {ox:.2f} {oy:.2f} {oz:.2f}")
            aid += 1
            lines.append(f"{aid} 2 1.008 {ox+0.76:.2f} {oy:.2f} {oz:.2f}")
            aid += 1
            lines.append(f"{aid} 2 1.008 {ox:.2f} {oy+0.76:.2f} {oz:.2f}")
            aid += 1
    path.write_text("\n".join(lines) + "\n")


def test_lammps_order_from_dump_smoke():
    """Run lammps_order_from_dump on minimal dump + data; check CSV schema and finite order params."""
    pytest.importorskip("scipy")
    from lammps_order_from_dump import main as run_script

    with tempfile.TemporaryDirectory() as tmp:
        data_path = Path(tmp) / "data.ice_uma"
        dump_path = Path(tmp) / "dump.lammpstrj"
        out_path = Path(tmp) / "order.csv"
        _minimal_data_file(data_path)
        _minimal_lammpstrj(dump_path, n_frames=2)

        argv_orig = sys.argv
        try:
            sys.argv = [
                "lammps_order_from_dump.py",
                "--data", str(data_path),
                "--dump", str(dump_path),
                "-o", str(out_path),
                "--every", "1",
            ]
            run_script()
        finally:
            sys.argv = argv_orig

        assert out_path.is_file()
        with out_path.open() as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 2
        assert "step" in rows[0] and "q6_mean" in rows[0] and "q_tet_mean" in rows[0]
        for r in rows:
            q6 = float(r["q6_mean"]) if r["q6_mean"] else float("nan")
            qt = float(r["q_tet_mean"]) if r["q_tet_mean"] else float("nan")
            assert (q6 >= 0 and q6 <= 2.0) or (qt >= 0 and qt <= 1.0), "at least one order param finite"
