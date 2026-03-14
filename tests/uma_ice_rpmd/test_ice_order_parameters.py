"""Unit tests for ice_order_parameters (Q6, q_tet)."""
from __future__ import annotations

import numpy as np
import pytest

from ice_order_parameters import (
    _mic,
    ice_order_metrics,
    steinhardt_q6_per_particle,
    tetrahedral_order_per_particle,
)


def test_mic_wrap():
    dr = np.array([[6.0, 0.0, 0.0]])
    box = np.array([10.0, 10.0, 10.0])
    assert np.allclose(_mic(dr, box), [[-4.0, 0.0, 0.0]])


def test_tetrahedral_order_perfect():
    pytest.importorskip("scipy")
    # Perfect tetrahedron around origin; fourth vertex from 3 known
    a = np.array(
        [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ],
        dtype=np.float64,
    )
    a /= np.linalg.norm(a[0])
    center = np.array([[0.0, 0.0, 0.0]])
    pos = np.vstack([center, a + center])
    box = np.array([20.0, 20.0, 20.0])
    q = tetrahedral_order_per_particle(pos, box, n_neigh=4)
    assert np.isfinite(q[0])
    assert q[0] > 0.99


def test_q6_finite_small_crystal():
    pytest.importorskip("scipy")
    # 2x2x2 simple cubic O sublattice spacing 3 Å (not ice, but Q6 should be finite)
    pts = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                pts.append([i * 3.0, j * 3.0, k * 3.0])
    pos = np.array(pts, dtype=np.float64)
    box = np.array([20.0, 20.0, 20.0])
    # 3 Å grid: each corner has only 3 face neighbors within 3.5 Å; need r_cut ≥ body diag
    q6 = steinhardt_q6_per_particle(pos, box, r_cut=6.0, r_min=0.5)
    assert np.all(np.isfinite(q6))
    assert np.mean(q6) > 0.1


def test_ice_order_metrics_ohh():
    pytest.importorskip("scipy")
    n_mol = 8
    pos = np.zeros((n_mol * 3, 3))
    for m in range(n_mol):
        pos[3 * m] = [m * 2.8, 0.0, 0.0]
        pos[3 * m + 1] = pos[3 * m] + [0.76, 0.0, 0.0]
        pos[3 * m + 2] = pos[3 * m] + [0.0, 0.76, 0.0]
    box = np.array([50.0, 50.0, 50.0])
    r = ice_order_metrics(pos, box)
    assert r.n_oxygen == n_mol
    assert r.n_q_tet_valid >= 0
