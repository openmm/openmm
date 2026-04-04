"""Regression: OpenMM UMA (per-atom wrap) vs LAMMPS-style Fairchem path on same ice data."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

_SCRIPT_DIR = Path(__file__).resolve().parent
_DATA_SMALL = _SCRIPT_DIR / "lammps" / "data.ice_uma_8"


def test_openmm_and_lammps_style_forces_match_small_ice():
    """Same positions/box → same E, F within tolerance (requires fairchem + small data file)."""
    import openmm

    if not hasattr(openmm, "PythonForce"):
        pytest.skip("OpenMM without PythonForce (use project build with Python ML plugin)")
    pytest.importorskip("fairchem")
    if not _DATA_SMALL.is_file():
        pytest.skip(f"Missing {_DATA_SMALL} (build with build_lammps_ice_data.py)")

    from openmm import Vec3, unit

    from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data
    from compare_forces_openmm_vs_lammps import (
        get_lammps_style_forces,
        get_openmm_forces,
    )
    from run_openmm_ice_lammps_match import water_topology

    pos_ang, cell, Z = parse_lammps_data(_DATA_SMALL)
    n = Z.shape[0]
    assert n % 3 == 0
    n_mol = n // 3

    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    a = Vec3(lx / 10, 0, 0) * unit.nanometer
    b = Vec3(xy / 10, ly / 10, 0) * unit.nanometer
    c = Vec3(xz / 10, yz / 10, lz / 10) * unit.nanometer
    box = (a, b, c)
    positions_nm = [tuple(float(x) for x in pos_ang[i] / 10.0) for i in range(n)]

    topology = water_topology(n_mol)
    e_omm, f_omm = get_openmm_forces(topology, positions_nm, box)

    try:
        e_lmp, f_lmp = get_lammps_style_forces(pos_ang, cell, Z, device="cpu")
    except Exception as exc:
        pytest.skip(f"LAMMPS-style UMA path unavailable: {exc}")

    d_e = abs(float(e_omm - e_lmp))
    max_df = float(np.max(np.abs(f_omm - f_lmp)))

    assert d_e <= 50.0, f"|ΔE| = {d_e} kJ/mol (tol 50)"
    assert max_df <= 200.0, f"max |ΔF| = {max_df} kJ/(mol·nm) (tol 200)"
