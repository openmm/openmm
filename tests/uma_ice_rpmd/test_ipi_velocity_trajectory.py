"""i-PI velocity trajectory parsing and centroid T_K parity helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ipi_order_from_traj import (
    AU_VELOCITY_TO_NM_PS,
    _iter_ipi_xyz_velocity_frames,
    _pos_traj_to_vel_traj_path,
)
from rpmd_thermo_utils import R_KJ_PER_MOL_K, centroid_kinetic_energy_and_temperature


def test_pos_traj_to_vel_traj_path() -> None:
    p = Path("ipi/ice__i-pi.traj_00.xyz")
    assert _pos_traj_to_vel_traj_path(p).name == "ice__i-pi.vtraj_00.xyz"


def test_velocity_frame_atomic_unit_to_nm_ps(tmp_path: Path) -> None:
    v_au = 1.0e-6
    xyz = f"""1
# Step: 0  Bead: 0 velocities{{atomic_unit}}  time{{picosecond}}: 0.0
O {v_au} 0.0 0.0
"""
    p = tmp_path / "b0.xyz"
    p.write_text(xyz, encoding="utf-8")
    rows = list(_iter_ipi_xyz_velocity_frames(p, dt_fs=0.5))
    assert len(rows) == 1
    assert rows[0][2][0, 0] == pytest.approx(v_au * AU_VELOCITY_TO_NM_PS)


def test_centroid_temperature_two_beads_one_oxygen_matches_formula() -> None:
    """Same estimator as OpenMM: mean bead velocity, then ½m|v|² and equipartition."""
    masses = np.array([15.99943], dtype=np.float64)
    ndof = 3
    # v_centroid = (1,0,0) nm/ps from v0=(2,0,0), v1=(0,0,0)
    vel = np.array([[[2.0, 0.0, 0.0]], [[0.0, 0.0, 0.0]]], dtype=np.float64)
    _ke, t_k = centroid_kinetic_energy_and_temperature(
        vel, masses, ndof_centroid_ke=ndof
    )
    ke_exp = 0.5 * masses[0] * 1.0**2
    t_exp = (2.0 * ke_exp) / (ndof * R_KJ_PER_MOL_K)
    assert _ke == pytest.approx(ke_exp)
    assert t_k == pytest.approx(t_exp)
