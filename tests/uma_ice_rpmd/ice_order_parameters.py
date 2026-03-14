#!/usr/bin/env python3
"""
Ice / water structural order parameters (oxygen-only, PBC).

Useful to track **melting or amorphization** during MD: ordered ice tends to have
higher **mean Q6** (Steinhardt bond-orientational order) and higher **mean q_tet**
(tetrahedral order) than disordered liquid-like regions.

References
----------
- Steinhardt et al., PRB 28, 784 (1983) — Q_l bond-orientational order.
- Chau & Harding, J. Phys. Chem. B 110, 16731 (2006) — tetrahedral order q.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

def _sph_harm_batch(l: int, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Y_l^m for m=-l..l; shape (2l+1, n). SciPy >=1.14: sph_harm_y(n,m,theta,phi)."""
    try:
        from scipy.special import sph_harm_y
    except ImportError as e:
        raise ImportError("ice_order_parameters Q6 needs scipy.special.sph_harm_y") from e
    out = np.zeros((2 * l + 1, len(theta)), dtype=np.complex128)
    for k, m in enumerate(range(-l, l + 1)):
        out[k] = sph_harm_y(l, m, theta, phi)
    return out


def _mic(dr: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Minimum-image convention; box = (Lx, Ly, Lz) orthorhombic."""
    out = np.asarray(dr, dtype=np.float64).copy()
    for k in range(3):
        L = box[k]
        out[..., k] -= L * np.round(out[..., k] / L)
    return out


def oxygen_positions(pos_ang: np.ndarray, z: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Oxygen Cartesian coordinates (Å). Assumes O,H,H per molecule if z not given.
    """
    pos_ang = np.asarray(pos_ang, dtype=np.float64)
    n = pos_ang.shape[0]
    if z is not None:
        z = np.asarray(z)
        o_mask = z == 8
        return pos_ang[o_mask]
    if n % 3 != 0:
        raise ValueError("Expected N divisible by 3 (O,H,H) or pass z for O only")
    return pos_ang[::3].copy()


@dataclass
class IceOrderResult:
    """Mean/std over oxygen atoms (nan if no valid sites)."""

    q6_mean: float
    q6_std: float
    q_tet_mean: float
    q_tet_std: float
    n_oxygen: int
    n_q6_valid: int
    n_q_tet_valid: int


def steinhardt_q6_per_particle(
    pos_o_ang: np.ndarray,
    box_ang: np.ndarray,
    r_cut: float = 3.7,
    r_min: float = 2.0,
    l: int = 6,
) -> np.ndarray:
    """
    Steinhardt Q_l for each oxygen (l=6 sensitive to crystalline order).

    Parameters
    ----------
    pos_o_ang
        (N, 3) O positions in Å.
    box_ang
        (3,) orthorhombic lengths Lx, Ly, Lz in Å.
    r_cut, r_min
        First-shell window (Å). Ice Ih first O–O peak ~2.76 Å; exclude r < r_min
        to avoid double-counting bonded H neighbors if ever mixed.
    """
    pos = np.asarray(pos_o_ang, dtype=np.float64)
    box = np.asarray(box_ang, dtype=np.float64).ravel()[:3]
    n = pos.shape[0]
    q_out = np.full(n, np.nan, dtype=np.float64)
    if n < 2:
        return q_out
    # Precompute all pairwise MIC distances
    for i in range(n):
        dr = _mic(pos - pos[i], box)
        d = np.linalg.norm(dr, axis=1)
        mask = (d > r_min) & (d < r_cut)
        mask[i] = False
        idx = np.where(mask)[0]
        if idx.size < 4:
            continue
        qh = dr[idx]
        dist = np.linalg.norm(qh, axis=1, keepdims=True)
        dist[dist < 1e-12] = 1.0
        rhat = qh / dist
        # Colatitude theta from +z, phi from +x
        x, y, zc = rhat[:, 0], rhat[:, 1], rhat[:, 2]
        theta = np.arccos(np.clip(zc, -1.0, 1.0))
        phi = np.arctan2(y, x)
        ylm_block = _sph_harm_batch(l, theta, phi)
        acc = np.mean(ylm_block, axis=1)
        pref = 4.0 * np.pi / (2 * l + 1)
        q_out[i] = float(np.sqrt(pref * np.sum(np.abs(acc) ** 2)))
    return q_out


def tetrahedral_order_per_particle(
    pos_o_ang: np.ndarray,
    box_ang: np.ndarray,
    n_neigh: int = 4,
) -> np.ndarray:
    """
    Chau–Harding tetrahedral order q in [0, 1]; 1 = perfect tetrahedron of 4 neighbors.

    Uses the **n_neigh** nearest oxygen neighbors (minimum image).
    """
    pos = np.asarray(pos_o_ang, dtype=np.float64)
    box = np.asarray(box_ang, dtype=np.float64).ravel()[:3]
    n = pos.shape[0]
    q_out = np.full(n, np.nan, dtype=np.float64)
    if n <= n_neigh:
        return q_out
    for i in range(n):
        dr = _mic(pos - pos[i], box)
        d = np.linalg.norm(dr, axis=1)
        d[i] = np.inf
        idx = np.argpartition(d, n_neigh)[:n_neigh]
        if np.any(d[idx] >= 1e6):
            continue
        vecs = dr[idx]
        nv = np.linalg.norm(vecs, axis=1, keepdims=True)
        nv[nv < 1e-12] = 1.0
        u = vecs / nv
        s = 0.0
        for a in range(n_neigh):
            for b in range(a + 1, n_neigh):
                c = float(np.clip(np.dot(u[a], u[b]), -1.0, 1.0))
                s += (c + 1.0 / 3.0) ** 2
        q_out[i] = 1.0 - 3.0 / 8.0 * s
    return q_out


def ice_order_metrics(
    pos_all_ang: np.ndarray,
    box_orth_ang: np.ndarray,
    z: Optional[np.ndarray] = None,
    r_cut_q6: float = 3.7,
    r_min_q6: float = 2.0,
) -> IceOrderResult:
    """
    Aggregate Q6 and q_tet over all oxygens with enough neighbors.

    box_orth_ang : (Lx, Ly, Lz) in Å.
    """
    pos_o = oxygen_positions(pos_all_ang, z)
    box = np.asarray(box_orth_ang, dtype=np.float64).ravel()[:3]
    q6 = steinhardt_q6_per_particle(pos_o, box, r_cut=r_cut_q6, r_min=r_min_q6)
    qt = tetrahedral_order_per_particle(pos_o, box, n_neigh=4)
    m6 = np.isfinite(q6)
    mt = np.isfinite(qt)
    return IceOrderResult(
        q6_mean=float(np.nanmean(q6)) if np.any(m6) else float("nan"),
        q6_std=float(np.nanstd(q6)) if np.any(m6) else float("nan"),
        q_tet_mean=float(np.nanmean(qt)) if np.any(mt) else float("nan"),
        q_tet_std=float(np.nanstd(qt)) if np.any(mt) else float("nan"),
        n_oxygen=int(pos_o.shape[0]),
        n_q6_valid=int(np.sum(m6)),
        n_q_tet_valid=int(np.sum(mt)),
    )


def positions_nm_to_oxygen_angstrom(pos_nm) -> Tuple[np.ndarray, np.ndarray]:
    """OpenMM positions (nm, all atoms O,H,H) → full (n,3) Å + orthorhombic box from state."""
    from openmm import unit

    pos_ang = np.array(
        [[p[i].value_in_unit(unit.nanometer) * 10.0 for i in range(3)] for p in pos_nm],
        dtype=np.float64,
    )
    return pos_ang
