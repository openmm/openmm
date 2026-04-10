"""Molecular PBC wrap for UMA: per-atom wrap breaks O–H bonds in unwrapped OpenMM coords."""
from __future__ import annotations

import numpy as np
import pytest


def wrap_water_molecules_per_bead(pos_angstrom: np.ndarray, cell_3x3: np.ndarray) -> np.ndarray:
    """Mirror of openmmml _wrap_water_molecules_per_bead (one bead’s config; O in cell, H min-image)."""
    n = pos_angstrom.shape[0]
    assert n % 3 == 0
    cell = np.asarray(cell_3x3, dtype=np.float64)
    inv = np.linalg.inv(cell)
    out = np.array(pos_angstrom, dtype=np.float64, copy=True)
    for m in range(n // 3):
        o = out[3 * m]
        fo = o @ inv.T
        fo = fo - np.floor(fo)
        o_n = fo[0] * cell[0] + fo[1] * cell[1] + fo[2] * cell[2]
        out[3 * m] = o_n
        for k in (1, 2):
            d = out[3 * m + k] - o
            df = d @ inv.T
            df = df - np.round(df)
            out[3 * m + k] = o_n + df[0] * cell[0] + df[1] * cell[1] + df[2] * cell[2]
    return out


def test_intra_molecular_oh_distance_unchanged():
    cell = np.diag([10.0, 10.0, 10.0])
    o = np.array([14.5, 1.0, 1.0])
    h1 = np.array([15.0, 1.0, 1.0])
    h2 = np.array([14.5, 1.5, 1.0])
    pos = np.vstack([o, h1, h2])
    d_before = np.linalg.norm(pos[1] - pos[0])
    w = wrap_water_molecules_per_bead(pos, cell)
    d_after = np.linalg.norm(w[1] - w[0])
    assert abs(d_before - d_after) < 1e-10


def test_rpmd_two_beads_independent():
    """Each bead wrapped with its own coords only (no cross-bead mixing)."""
    cell = np.diag([10.0, 10.0, 10.0])
    # Bead 0: one water unwrapped
    b0 = np.vstack(
        [
            np.array([14.5, 1.0, 1.0]),
            np.array([15.0, 1.0, 1.0]),
            np.array([14.5, 1.5, 1.0]),
        ]
    )
    # Bead 1: different unwrap (same molecule topology)
    b1 = b0 + np.array([50.0, 0.0, 0.0])
    w0 = wrap_water_molecules_per_bead(b0, cell)
    w1 = wrap_water_molecules_per_bead(b1, cell)
    assert np.allclose(w0, w1, atol=1e-9), "same relative geometry → same wrapped cluster"


def test_per_atom_wrap_would_break_bond():
    cell = np.diag([10.0, 10.0, 10.0])
    from ase.geometry import wrap_positions
    o = np.array([0.5, 5.0, 5.0])
    h1 = np.array([9.8, 5.0, 5.0])  # min-image from O is 0.3 Å in x
    h2 = np.array([0.5, 5.5, 5.0])
    pos = np.vstack([o, h1, h2])
    bad = wrap_positions(pos, cell, pbc=True, eps=0)
    assert np.linalg.norm(bad[1] - bad[0]) > 5.0
    good = wrap_water_molecules_per_bead(pos, cell)
    assert np.linalg.norm(good[1] - good[0]) < 1.0
    assert np.linalg.norm(good[2] - good[0]) < 1.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
