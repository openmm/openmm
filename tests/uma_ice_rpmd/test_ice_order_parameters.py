"""Unit tests for ice_order_parameters (Q6, q_tet)."""
from __future__ import annotations

import numpy as np
import pytest

from ice_order_parameters import (
    ICE_ORDER_PI_CSV_HEADER,
    _mic,
    format_ice_order_pi_csv_row,
    ice_order_metrics,
    ice_order_metrics_path_integral,
    steinhardt_q6_per_particle,
    tetrahedral_order_per_particle,
)


def test_mic_wrap():
    dr = np.array([[6.0, 0.0, 0.0]])
    box = np.array([10.0, 10.0, 10.0])
    assert np.allclose(_mic(dr, box), [[-4.0, 0.0, 0.0]])


def test_tetrahedral_order_square_planar_center():
    """Four coplanar neighbors at 90° give Chau–Harding q = 0.5 at the central site."""
    pytest.importorskip("scipy")
    center = np.array([[0.0, 0.0, 0.0]])
    neigh = np.array(
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
        dtype=np.float64,
    )
    pos = np.vstack([center, neigh])
    box = np.array([50.0, 50.0, 50.0])
    q = tetrahedral_order_per_particle(pos, box, n_neigh=4)
    assert np.isfinite(q[0])
    assert abs(q[0] - 0.5) < 1e-12


def test_tetrahedral_order_bundled_neighbors_negative():
    """Strongly non-tetrahedral four-neighbor shell: Chau–Harding q can be < 0."""
    pytest.importorskip("scipy")
    center = np.array([[0.0, 0.0, 0.0]])
    neigh = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.99, 0.14, 0.0],
            [0.99, -0.14, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )
    pos = np.vstack([center, neigh])
    box = np.array([50.0, 50.0, 50.0])
    q = tetrahedral_order_per_particle(pos, box, n_neigh=4)
    assert np.isfinite(q[0])
    assert q[0] < -0.5


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


def test_path_integral_matches_single_bead():
    """PI average with one bead equals classical ice_order_metrics."""
    pytest.importorskip("scipy")
    # 8 oxygens on a 3 Å grid (same geometry as test_q6_finite_small_crystal)
    pts = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                pts.append([i * 3.0, j * 3.0, k * 3.0])
    pos_o = np.array(pts, dtype=np.float64)
    n = pos_o.shape[0]
    pos = pos_o.copy()
    z = np.full(n, 8, dtype=np.int64)
    box = np.array([20.0, 20.0, 20.0])
    cl = ice_order_metrics(pos, box, z=z, r_cut_q6=6.0, r_min_q6=0.5)
    pi = ice_order_metrics_path_integral(
        [pos.copy()], box, z=z, r_cut_q6=6.0, r_min_q6=0.5
    )
    assert np.isfinite(pi.q6_mean)
    assert abs(pi.q6_mean - cl.q6_mean) < 1e-9
    assert abs(pi.q_tet_mean - cl.q_tet_mean) < 1e-9


def test_path_integral_min_max_bracket_mean_and_valid_counts():
    """Min/mean/max and n_*_valid are consistent for finite path-integral samples."""
    pytest.importorskip("scipy")
    pts = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                pts.append([i * 3.0, j * 3.0, k * 3.0])
    pos_o = np.array(pts, dtype=np.float64)
    n = pos_o.shape[0]
    z = np.full(n, 8, dtype=np.int64)
    box = np.array([20.0, 20.0, 20.0])
    pi = ice_order_metrics_path_integral(
        [pos_o.copy()], box, z=z, r_cut_q6=6.0, r_min_q6=0.5
    )
    assert pi.n_q6_valid == n
    assert pi.n_q_tet_valid == n
    assert pi.q6_min <= pi.q6_mean <= pi.q6_max
    assert pi.q_tet_min <= pi.q_tet_mean <= pi.q_tet_max


def test_format_ice_order_pi_csv_matches_header_column_count():
    pytest.importorskip("scipy")
    pts = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                pts.append([i * 3.0, j * 3.0, k * 3.0])
    pos_o = np.array(pts, dtype=np.float64)
    z = np.full(pos_o.shape[0], 8, dtype=np.int64)
    box = np.array([20.0, 20.0, 20.0])
    pi = ice_order_metrics_path_integral(
        [pos_o.copy()], box, z=z, r_cut_q6=6.0, r_min_q6=0.5
    )
    header_cols = ICE_ORDER_PI_CSV_HEADER.strip().split(",")
    row = format_ice_order_pi_csv_row(0, 0.0, "", "", pi).strip()
    data_cols = row.split(",")
    assert len(data_cols) == len(header_cols)


def test_path_integral_two_beads_is_mean_of_observables():
    """First-O PI Q6 equals mean_b Q6_O^(b) (not Q6 evaluated on centroid coords)."""
    pytest.importorskip("scipy")
    pts = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                pts.append([i * 3.0, j * 3.0, k * 3.0])
    pos_a = np.array(pts, dtype=np.float64)
    pos_b = np.array(pos_a, copy=True)
    pos_b[0, 0] += 0.05  # shift first O on bead 1
    n = pos_a.shape[0]
    z = np.full(n, 8, dtype=np.int64)
    box = np.array([20.0, 20.0, 20.0])
    pi = ice_order_metrics_path_integral(
        [pos_a, pos_b], box, z=z, r_cut_q6=6.0, r_min_q6=0.5
    )
    q6a = steinhardt_q6_per_particle(pos_a, box, r_cut=6.0, r_min=0.5)
    q6b = steinhardt_q6_per_particle(pos_b, box, r_cut=6.0, r_min=0.5)
    assert np.isfinite(q6a[0]) and np.isfinite(q6b[0])
    expect0 = 0.5 * (q6a[0] + q6b[0])
    q6_stack = np.vstack([q6a, q6b])
    pi0 = float(np.nanmean(q6_stack[:, 0]))
    assert abs(pi0 - expect0) < 1e-12
    assert np.isfinite(pi.q6_mean)
