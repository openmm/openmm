"""Unit tests for i-PI multi-bead trajectory path discovery."""
from __future__ import annotations

from pathlib import Path

import pytest

from ipi_order_from_traj import _traj_paths_from_prefix


def test_traj_paths_zero_padded() -> None:
    base = Path("ipi/ice__i-pi.traj_00.xyz")
    paths = _traj_paths_from_prefix(base, 4)
    names = [p.name for p in paths]
    assert names == [
        "ice__i-pi.traj_00.xyz",
        "ice__i-pi.traj_01.xyz",
        "ice__i-pi.traj_02.xyz",
        "ice__i-pi.traj_03.xyz",
    ]


def test_traj_paths_unpadded() -> None:
    base = Path("ipi/ice__i-pi.traj_0.xyz")
    paths = _traj_paths_from_prefix(base, 3)
    names = [p.name for p in paths]
    assert names == [
        "ice__i-pi.traj_0.xyz",
        "ice__i-pi.traj_1.xyz",
        "ice__i-pi.traj_2.xyz",
    ]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
