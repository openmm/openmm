"""Utilities for RPMD thermo reporting shared by drivers and tests."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.random import Generator

R_KJ_PER_MOL_K = 0.00831446261815324
# SI constants for velocity draws (nm / ps)
_BOLTZMANN_J_PER_K = 1.380649e-23
_AMU_KG = 1.66053906660e-27
_NM_PER_M = 1.0e9
_PS_PER_S = 1.0e12
_MPS_TO_NM_PS = _NM_PER_M / _PS_PER_S  # 1e-3


def ring_polymer_nm_matrix(nbeads: int) -> np.ndarray:
    """Orthonormal bead↔normal-mode matrix **C** with **q_nm = C @ q_bead** (columns).

    Matches i-PI ``mk_nm_matrix`` / closed-path RPMD convention
    (`ipi.utils.nmtransform.mk_nm_matrix`). Velocities transform the same way:
    **v_bead = C.T @ v_nm** for equal bead masses.

    Parameters
    ----------
    nbeads
        Number of ring-polymer replicas (imaginary-time slices).

    Returns
    -------
    numpy.ndarray
        Shape ``(nbeads, nbeads)``, ``C @ C.T == I``.
    """

    n = int(nbeads)
    if n < 1:
        raise ValueError("nbeads must be >= 1")
    b2nm = np.zeros((n, n), dtype=np.float64)
    b2nm[0, :] = math.sqrt(1.0)
    for j in range(n):
        for i in range(1, n // 2 + 1):
            b2nm[i, j] = math.sqrt(2.0) * math.cos(2.0 * math.pi * j * i / float(n))
        for i in range(n // 2 + 1, n):
            b2nm[i, j] = math.sqrt(2.0) * math.sin(2.0 * math.pi * j * i / float(n))
    if n % 2 == 0:
        b2nm[n // 2, 0:n:2] = 1.0
        b2nm[n // 2, 1:n:2] = -1.0
    return b2nm / math.sqrt(float(n))


def _per_component_velocity_std_nm_ps(mass_da: float, temperature_k: float) -> float:
    """``sqrt(k_B T / m)`` in nm/ps for one Cartesian component (single bead)."""
    if mass_da <= 0.0:
        return 0.0
    m_kg = float(mass_da) * _AMU_KG
    var_m2_s2 = _BOLTZMANN_J_PER_K * float(temperature_k) / m_kg
    sigma_m_s = math.sqrt(var_m2_s2)
    return sigma_m_s * _MPS_TO_NM_PS


def sample_rpmd_velocities_nm_ps(
    rng: Generator,
    masses_da: np.ndarray | list,
    *,
    n_beads: int,
    temperature_k: float,
) -> np.ndarray:
    """Thermal RPMD velocities matching i-PI ``<velocities mode='thermal'>`` NM draw.

    i-PI initializes normal-mode momenta as independent Gaussians with variance
    ``m * nbeads * k_B * T`` per component, then maps to bead momenta. For equal
    bead masses this is equivalent to drawing each NM *velocity* component as
    ``N(0, nbeads * k_B T / m)`` and transforming **v_bead = C.T @ v_nm**
    (`ring_polymer_nm_matrix`).

    Parameters
    ----------
    rng
        NumPy random generator (use ``np.random.default_rng(seed)``).
    masses_da
        Per-atom masses in daltons, shape ``(n_atoms,)``.
    n_beads
        Number of ring-polymer beads.
    temperature_k
        Physical bath temperature in kelvin (same as i-PI ``ensemble.temp``).

    Returns
    -------
    numpy.ndarray
        Velocities with shape ``(n_beads, n_atoms, 3)`` in nm/ps.
    """

    masses = np.asarray(masses_da, dtype=np.float64)
    n_atoms = masses.shape[0]
    n = int(n_beads)
    if n < 1:
        raise ValueError("n_beads must be >= 1")
    c = ring_polymer_nm_matrix(n)
    out = np.zeros((n, n_atoms, 3), dtype=np.float64)
    scale_per_atom = np.array(
        [
            _per_component_velocity_std_nm_ps(float(masses[a]), temperature_k) * math.sqrt(float(n))
            for a in range(n_atoms)
        ],
        dtype=np.float64,
    )
    for a in range(n_atoms):
        sig = scale_per_atom[a]
        if sig == 0.0:
            continue
        for xyz in range(3):
            v_nm = rng.standard_normal(n) * sig
            out[:, a, xyz] = c.T @ v_nm
    return out


def initialize_rpmd_velocities_nm(
    integrator,
    masses_da: np.ndarray | list,
    *,
    n_beads: int,
    temperature_k: float,
    seed: int | None = None,
    rng: Generator | None = None,
) -> None:
    """Set RPMD bead velocities from i-PI-equivalent NM thermal sampling.

    Parameters
    ----------
    integrator
        OpenMM ``RPMDIntegrator`` (or compatible) with ``setVelocities(bead, list)``.
    masses_da
        Same atom order as the integrator / topology, masses in daltons.
    n_beads, temperature_k
        Ring size and bath temperature.
    seed
        Used if ``rng`` is None: ``np.random.default_rng(seed)``.
    rng
        Optional pre-built generator (overrides ``seed``).
    """

    from openmm import Vec3, unit

    if rng is None:
        rng = np.random.default_rng(seed)
    v_arr = sample_rpmd_velocities_nm_ps(
        rng,
        masses_da,
        n_beads=n_beads,
        temperature_k=temperature_k,
    )
    vel_unit = unit.nanometer / unit.picosecond
    for b in range(n_beads):
        vlist = [
            Vec3(
                float(v_arr[b, a, 0]),
                float(v_arr[b, a, 1]),
                float(v_arr[b, a, 2]),
            )
            * vel_unit
            for a in range(v_arr.shape[1])
        ]
        integrator.setVelocities(b, vlist)


def centroid_kinetic_energy_and_temperature(
    all_bead_velocities_nm_ps: np.ndarray | list,
    masses_da: np.ndarray | list,
    *,
    ndof_centroid_ke: int,
) -> tuple[float, float]:
    """Return centroid kinetic energy and derived temperature.

    Parameters
    ----------
    all_bead_velocities_nm_ps
        Bead velocities with shape ``(n_beads, n_atoms, 3)`` in ``nm/ps``.
    masses_da
        Atomic masses in dalton for the corresponding atoms.
    ndof_centroid_ke
        Number of massive centroid degrees of freedom used for equipartition.
    """

    bead_velocities = np.asarray(all_bead_velocities_nm_ps, dtype=np.float64)
    masses = np.asarray(masses_da, dtype=np.float64)
    centroid_velocities = np.mean(bead_velocities, axis=0)
    centroid_ke = float(
        0.5 * np.sum(masses[:, None] * centroid_velocities * centroid_velocities)
    )
    temperature_k = (2.0 * centroid_ke) / (float(ndof_centroid_ke) * R_KJ_PER_MOL_K)
    return centroid_ke, temperature_k
