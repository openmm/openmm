"""Tests for i-PI xyz parsing (cell/position units) in ipi_order_from_traj."""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(_SCRIPT_DIR))

from benchmark_same_structure_openmm_vs_lammps import parse_lammps_data
from ice_order_parameters import wrap_cartesian_orthorhombic
from ipi_order_from_traj import (
    BOHR_TO_ANG,
    _cell_abc_to_angstrom,
    _iter_ipi_xyz_frames,
    _positions_are_bohr,
)


def test_cell_abc_atomic_unit_converts_bohr_to_angstrom():
    comment = "# CELL(abcABC): 16.99620 29.43815 27.67315 90 90 90  cell{atomic_unit}"
    a, b, c = _cell_abc_to_angstrom(comment, 16.99620, 29.43815, 27.67315)
    np.testing.assert_allclose([a, b, c], [8.994004, 15.578004, 14.644004], rtol=1e-5)


def test_cell_abc_angstrom_tag_unchanged():
    comment = "# CELL(abcABC): 8.994 15.578 14.644 90 90 90  cell{angstrom}"
    a, b, c = _cell_abc_to_angstrom(comment, 8.994, 15.578, 14.644)
    assert (a, b, c) == (8.994, 15.578, 14.644)


def test_time_ps_from_step_and_dt_fs_matches_ipi():
    """time_ps uses Step: × dt_fs / 1000 (not frame index × wrong default)."""
    lines = [
        "3",
        "# CELL(abcABC): 10 10 10 90 90 90 Step: 10000 Bead: 0 positions{atomic_unit} cell{atomic_unit}",
        "O 0 0 0",
        "H 0 0 0",
        "H 0 0 0",
    ]
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
        f.write("\n".join(lines) + "\n")
        p = Path(f.name)
    try:
        step, t_ps, *_ = next(iter(_iter_ipi_xyz_frames(p, dt_fs=0.1)))
        assert step == 10000
        assert abs(t_ps - 1.0) < 1e-9
        step2, t2, *_ = next(iter(_iter_ipi_xyz_frames(p, dt_fs=0.5)))
        assert abs(t2 - 5.0) < 1e-9
    finally:
        p.unlink(missing_ok=True)


def test_iter_ipi_xyz_first_frame_box_matches_lammps_cell():
    """Synthetic frame: Bohr positions + cell{atomic_unit} → Å box like ice supercell."""
    # One H2O: O at origin, H dummy (ignored for box)
    o_bohr = np.array([0.0, 4.90635886, 0.86063493])  # ~ first O from small ice in Bohr
    lines = [
        "3",
        f"# CELL(abcABC): 16.99620 29.43815 27.67315 90.00000 90.00000 90.00000  Step: 0  Bead: 0 positions{{atomic_unit}}  cell{{atomic_unit}}",
        f"O {o_bohr[0]:.8e} {o_bohr[1]:.8e} {o_bohr[2]:.8e}",
        "H 0.0 0.0 0.0",
        "H 0.0 0.0 0.0",
    ]
    with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
        f.write("\n".join(lines) + "\n")
        p = Path(f.name)
    try:
        step, _t, pos, box, z = next(iter(_iter_ipi_xyz_frames(p)))
        assert step == 0
        np.testing.assert_allclose(box, [8.994004, 15.578004, 14.644004], rtol=1e-4)
        np.testing.assert_allclose(pos[0], o_bohr * BOHR_TO_ANG, rtol=1e-6)
        assert z[0] == 8
    finally:
        p.unlink(missing_ok=True)


@pytest.mark.skipif(
    not (_SCRIPT_DIR / "ipi" / "ice__i-pi.traj_00.xyz").is_file(),
    reason="trajectory fixture not present",
)
def test_real_traj_header_box_orth_in_angstrom():
    traj = _SCRIPT_DIR / "ipi" / "ice__i-pi.traj_00.xyz"
    step, _t, _pos, box, _z = next(iter(_iter_ipi_xyz_frames(traj)))
    assert step == 0
    np.testing.assert_allclose(box, [8.994, 15.578, 14.644], rtol=1e-3)


@pytest.mark.skipif(
    not (_SCRIPT_DIR / "lammps" / "data.ice_uma_64").is_file(),
    reason="data.ice_uma_64 not present",
)
def test_traj_frame0_oxygens_match_lammps_data_after_wrap():
    """LAMMPS data vs i-PI bead-0 frame 0 after per-atom wrap: O sites match (MIC, label perm)."""
    pytest.importorskip("scipy")
    data = _SCRIPT_DIR / "lammps" / "data.ice_uma_64"
    traj = _SCRIPT_DIR / "ipi" / "ice__i-pi.traj_00.xyz"
    if not traj.is_file():
        pytest.skip("trajectory not present")

    pos_data, cell_mat, z_data = parse_lammps_data(data)
    box_data = np.array(
        [cell_mat[0, 0], cell_mat[1, 1], cell_mat[2, 2]], dtype=np.float64
    )
    o_mask = z_data == 8
    o_ref = pos_data[o_mask].copy()

    step, _t, pos_traj, box_traj, z_traj = next(iter(_iter_ipi_xyz_frames(traj)))
    assert step == 0
    np.testing.assert_allclose(box_traj, box_data, rtol=1e-4)

    pos_wrapped = wrap_cartesian_orthorhombic(pos_traj, box_traj)
    o_traj = pos_wrapped[z_traj == 8]

    assert o_ref.shape == o_traj.shape
    # Match each O to nearest periodic image of some ref O (label permutation)
    L = box_traj
    best = np.full(o_traj.shape[0], np.inf)
    for i in range(o_traj.shape[0]):
        for j in range(o_ref.shape[0]):
            d = o_traj[i] - o_ref[j]
            for k in range(3):
                d[k] -= L[k] * np.round(d[k] / L[k])
            best[i] = min(best[i], np.linalg.norm(d))
    # i-PI vs LAMMPS reference can differ by ~0.2 Å in label/MIC pairing for some checkpoints.
    assert np.max(best) < 0.25, f"max MIC distance O mismatch {np.max(best)} Å"


def test_positions_bohr_explicit_tag():
    assert _positions_are_bohr("positions{atomic_unit} cell{atomic_unit}")
    assert not _positions_are_bohr("# CELL cell{angstrom} no positions tag")


def test_positions_legacy_atomic_without_positions_brace():
    assert _positions_are_bohr("foo atomic_unit bar") is True
