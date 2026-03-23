#!/usr/bin/env python3
"""
Ice / water structural order parameters (oxygen-only, PBC).

Useful to track **melting or amorphization** during MD: ordered ice tends to have
higher **mean Q6** (Steinhardt bond-orientational order) and higher **mean q_tet**
(Chau–Harding tetrahedral order) than disordered liquid-like regions.

**q_tet range:** The implemented Chau–Harding formula is ``q = 1 - (3/8) Σ (cos θ + 1/3)²``
over the six angles between the four O–O neighbor unit vectors. It equals **1** for a
perfect tetrahedron (cos θ = -1/3) but is **not** bounded below by 0: strongly
non-tetrahedral four-neighbor geometries can yield **negative** *q*. Path-integral
CSVs therefore may show negative ``q_tet_p10`` (or min) when some oxygens have
bead-averaged ⟨q⟩ well below zero.

**RPMD / path integrals:** use :func:`ice_order_metrics_path_integral` — for each
oxygen *i*, compute the observable on **each bead**, then average over beads
⟨O⟩_i = (1/P) Σ_b O_i^(b), and report statistics over *i*. Do **not** evaluate
Q6 or q_tet on centroid coordinates (that biases nonlinear observables).

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


@dataclass
class IceOrderPiResult:
    """
    Path-integral / RPMD estimator for local order parameters.

    For each oxygen *i*, first compute Q6 and q_tet on **each bead** *b*, then
    form the bead average ⟨O_i⟩ = (1/P) Σ_b O_i^(b). **Do not** evaluate the
    observable on centroid coordinates (that would bias observables nonlinear
    in positions).

    Reported ``*_mean`` / ``*_std`` are over oxygens of those ⟨·⟩_bead values.
    Percentiles (p10/p50/p90) and ``*_min`` / ``*_max`` describe the **distribution
    across oxygens** of ⟨Q6⟩_bead,i and ⟨q_tet⟩_bead,i at this time slice. Because
    Chau–Harding *q* can be negative, ``q_tet_p10`` or ``q_tet_min`` may be below
    zero without indicating a code error.

    ``n_q6_valid`` / ``n_q_tet_valid`` count oxygens with finite bead-averaged Q6 /
    q_tet (same as the sample size for those statistics).
    """

    q6_mean: float
    q6_std: float
    q6_p10: float
    q6_p50: float
    q6_p90: float
    q6_min: float
    q6_max: float
    q_tet_mean: float
    q_tet_std: float
    q_tet_p10: float
    q_tet_p50: float
    q_tet_p90: float
    q_tet_min: float
    q_tet_max: float
    n_oxygen: int
    n_q6_valid: int
    n_q_tet_valid: int


# CSV header for RPMD path-integral order output (see ice_order_metrics_path_integral)
ICE_ORDER_PI_CSV_HEADER = (
    "step,time_ps,T_K,PE_kj_mol,"
    "q6_mean,q6_std,q6_p10,q6_p50,q6_p90,q6_min,q6_max,"
    "q_tet_mean,q_tet_std,q_tet_p10,q_tet_p50,q_tet_p90,q_tet_min,q_tet_max,"
    "n_q6_valid,n_q_tet_valid,n_oxygen\n"
)


def format_ice_order_pi_csv_row(
    step: int,
    time_ps: float,
    t_k: str,
    pe_kj_mol: str,
    res: IceOrderPiResult,
) -> str:
    """One CSV data line matching :data:`ICE_ORDER_PI_CSV_HEADER`."""
    return (
        f"{step},{time_ps:.6f},{t_k},{pe_kj_mol},"
        f"{res.q6_mean:.6f},{res.q6_std:.6f},{res.q6_p10:.6f},{res.q6_p50:.6f},{res.q6_p90:.6f},"
        f"{res.q6_min:.6f},{res.q6_max:.6f},"
        f"{res.q_tet_mean:.6f},{res.q_tet_std:.6f},{res.q_tet_p10:.6f},{res.q_tet_p50:.6f},{res.q_tet_p90:.6f},"
        f"{res.q_tet_min:.6f},{res.q_tet_max:.6f},"
        f"{res.n_q6_valid},{res.n_q_tet_valid},{res.n_oxygen}\n"
    )


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
    Chau–Harding tetrahedral order *q* from four nearest O–O neighbor directions.

    Uses ``q = 1 - (3/8) Σ_{a<b} (cos θ_ab + 1/3)²`` (minimum-image **n_neigh**
    nearest neighbors, default 4). For a **perfect tetrahedron**, cos θ = -1/3 for
    all pairs and *q* → **1**. For strongly **non-tetrahedral** arrangements (e.g.
    nearly coplanar or compressed neighbor shells), the sum can exceed 8/3 and
    **q can be negative**; there is no lower bound in this definition.

    Returns NaN for sites with fewer than ``n_neigh`` finite neighbors.
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


def _mean_over_beads_per_site(q_bead: np.ndarray) -> np.ndarray:
    """
    Mean along axis 0 (beads) using only finite values per oxygen.

    Columns that are all-NaN stay NaN **without** NumPy's "Mean of empty slice"
    warning (that happens when ``np.nanmean`` sees no finite entries in a slice).
    """
    q_bead = np.asarray(q_bead, dtype=np.float64)
    _, n_o = q_bead.shape
    out = np.full(n_o, np.nan, dtype=np.float64)
    for i in range(n_o):
        col = q_bead[:, i]
        ok = np.isfinite(col)
        if np.any(ok):
            out[i] = float(np.mean(col[ok]))
    return out


def ice_order_metrics_path_integral(
    bead_full_positions_ang: list[np.ndarray],
    box_orth_ang: np.ndarray,
    z: Optional[np.ndarray] = None,
    r_cut_q6: float = 3.7,
    r_min_q6: float = 2.0,
) -> IceOrderPiResult:
    """
    RPMD / path-integral estimator: **average observable over beads, then** summarize.

    For each bead *b*, compute per-oxygen Q6 and q_tet from that bead's positions.
    For each oxygen *i*, ⟨O⟩_i = mean_b O_i^(b) (ignoring NaNs per bead). Then
    ``q6_mean`` = mean_i ⟨Q6⟩_i, ``q6_std`` = std_i, and percentiles plus min/max
    are taken over the distribution {⟨Q6⟩_i}_i (same for q_tet). Chau–Harding
    *q* is not bounded below by 0, so ⟨q_tet⟩_i and CSV tails (e.g. ``q_tet_min``)
    may be negative.

    Parameters
    ----------
    bead_full_positions_ang
        One (n_atoms, 3) array per bead, Cartesian Å, same atom order each bead.
    box_orth_ang
        (Lx, Ly, Lz) orthorhombic lengths in Å.
    z
        Atomic numbers per atom (e.g. O=8, H=1, M=0 for TIP4P/Flex).
    """
    if len(bead_full_positions_ang) < 1:
        raise ValueError("bead_full_positions_ang must be non-empty")
    box = np.asarray(box_orth_ang, dtype=np.float64).ravel()[:3]
    pos0 = np.asarray(bead_full_positions_ang[0], dtype=np.float64)
    n_atoms = pos0.shape[0]
    for b, p in enumerate(bead_full_positions_ang):
        p = np.asarray(p, dtype=np.float64)
        if p.shape != (n_atoms, 3):
            raise ValueError(f"Bead {b}: expected shape ({n_atoms}, 3), got {p.shape}")

    pos_o0 = oxygen_positions(pos0, z)
    n_o = int(pos_o0.shape[0])
    n_beads = len(bead_full_positions_ang)
    q6_bead = np.full((n_beads, n_o), np.nan, dtype=np.float64)
    qt_bead = np.full((n_beads, n_o), np.nan, dtype=np.float64)

    for b, pos_all in enumerate(bead_full_positions_ang):
        pos_all = np.asarray(pos_all, dtype=np.float64)
        pos_o = oxygen_positions(pos_all, z)
        if pos_o.shape[0] != n_o:
            raise ValueError(f"Bead {b}: oxygen count {pos_o.shape[0]} != {n_o}")
        q6_bead[b] = steinhardt_q6_per_particle(pos_o, box, r_cut=r_cut_q6, r_min=r_min_q6)
        qt_bead[b] = tetrahedral_order_per_particle(pos_o, box, n_neigh=4)

    # Per-oxygen mean over beads: use only finite entries. If an oxygen is NaN on
    # *every* bead (e.g. Q6 undefined — fewer than 4 neighbors in the r_min–r_cut
    # shell on every replica), np.nanmean would warn "Mean of empty slice".
    q6_pi = _mean_over_beads_per_site(q6_bead)
    qt_pi = _mean_over_beads_per_site(qt_bead)

    def _summarize(
        v: np.ndarray,
    ) -> tuple[float, float, float, float, float, float, float, int]:
        m = np.isfinite(v)
        if not np.any(m):
            return (float("nan"),) * 7 + (0,)
        x = v[m]
        return (
            float(np.mean(x)),
            float(np.std(x)),
            float(np.percentile(x, 10)),
            float(np.percentile(x, 50)),
            float(np.percentile(x, 90)),
            float(np.min(x)),
            float(np.max(x)),
            int(x.size),
        )

    q6_m, q6_s, q6_p10, q6_p50, q6_p90, q6_lo, q6_hi, nq6 = _summarize(q6_pi)
    qt_m, qt_s, qt_p10, qt_p50, qt_p90, qt_lo, qt_hi, nqt = _summarize(qt_pi)
    return IceOrderPiResult(
        q6_mean=q6_m,
        q6_std=q6_s,
        q6_p10=q6_p10,
        q6_p50=q6_p50,
        q6_p90=q6_p90,
        q6_min=q6_lo,
        q6_max=q6_hi,
        q_tet_mean=qt_m,
        q_tet_std=qt_s,
        q_tet_p10=qt_p10,
        q_tet_p50=qt_p50,
        q_tet_p90=qt_p90,
        q_tet_min=qt_lo,
        q_tet_max=qt_hi,
        n_oxygen=n_o,
        n_q6_valid=nq6,
        n_q_tet_valid=nqt,
    )


def positions_nm_to_oxygen_angstrom(pos_nm) -> Tuple[np.ndarray, np.ndarray]:
    """OpenMM positions (nm, all atoms O,H,H) → full (n,3) Å + orthorhombic box from state."""
    from openmm import unit

    pos_ang = np.array(
        [[p[i].value_in_unit(unit.nanometer) * 10.0 for i in range(3)] for p in pos_nm],
        dtype=np.float64,
    )
    return pos_ang
