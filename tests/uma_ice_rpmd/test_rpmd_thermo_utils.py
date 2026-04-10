from __future__ import annotations

import numpy as np
import pytest

from rpmd_thermo_utils import (
    centroid_kinetic_energy_and_temperature,
    initialize_rpmd_velocities_nm,
    ring_polymer_nm_matrix,
    sample_rpmd_velocities_nm_ps,
)


def test_centroid_kinetic_energy_and_temperature_matches_manual_value() -> None:
    velocities_nm_ps = [
        [[2.0, 0.0, 0.0], [0.0, 4.0, 0.0]],
        [[0.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
    ]
    masses_da = [16.0, 1.0]

    kinetic_kj_mol, temperature_k = centroid_kinetic_energy_and_temperature(
        velocities_nm_ps,
        masses_da,
        ndof_centroid_ke=6,
    )

    # Centroid velocities are [1, 0, 0] and [0, 3, 0].
    expected_ke = 0.5 * 16.0 * 1.0**2 + 0.5 * 1.0 * 3.0**2
    expected_temp = (2.0 * expected_ke) / (6.0 * 0.00831446261815324)

    assert kinetic_kj_mol == pytest.approx(expected_ke)
    assert temperature_k == pytest.approx(expected_temp)


def test_ring_polymer_nm_matrix_is_orthogonal() -> None:
    for n in (4, 8, 16, 32):
        c = ring_polymer_nm_matrix(n)
        assert c.shape == (n, n)
        assert np.allclose(c @ c.T, np.eye(n), atol=1e-12, rtol=1e-12)
        assert np.allclose(c.T @ c, np.eye(n), atol=1e-12, rtol=1e-12)


def test_sample_rpmd_velocities_centroid_temperature_mean_near_bath() -> None:
    """NM thermal init: centroid kinetic T should scatter about ``temperature_k``."""
    rng = np.random.default_rng(284759)
    n_beads = 32
    temperature_k = 243.0
    masses_da = np.array([16.0, 1.0, 1.0], dtype=np.float64)
    n_atoms = masses_da.shape[0]
    ndof = 3 * n_atoms

    n_samples = 4000
    temps = []
    for _ in range(n_samples):
        v = sample_rpmd_velocities_nm_ps(
            rng,
            masses_da,
            n_beads=n_beads,
            temperature_k=temperature_k,
        )
        assert v.shape == (n_beads, n_atoms, 3)
        _, t_k = centroid_kinetic_energy_and_temperature(
            v,
            masses_da,
            ndof_centroid_ke=ndof,
        )
        temps.append(t_k)

    mean_t = float(np.mean(temps))
    std_t = float(np.std(temps, ddof=1))
    assert mean_t == pytest.approx(temperature_k, abs=8.0)
    # Instantaneous centroid T has large shot noise for O(10) DOF; only bound loosely.
    assert 15.0 < std_t < 250.0


def test_sample_rpmd_velocities_nonzero_internal_spread() -> None:
    """Internal ring-polymer modes: bead-to-bead spread must not be ~0 (regression)."""
    rng = np.random.default_rng(42)
    v = sample_rpmd_velocities_nm_ps(
        rng,
        np.array([16.0]),
        n_beads=16,
        temperature_k=300.0,
    )
    spread = float(np.mean(np.var(v[:, 0, :], axis=0)))
    assert spread > 1e-4


def test_initialize_rpmd_velocities_nm_populates_every_bead() -> None:
    """``initialize_rpmd_velocities_nm`` drives ``setVelocities`` for each bead index."""
    pytest.importorskip("openmm")

    class _FakeIntegrator:
        def __init__(self) -> None:
            self.calls: list[tuple[int, int]] = []

        def setVelocities(self, bead: int, vlist) -> None:  # noqa: N802 — OpenMM API
            self.calls.append((bead, len(vlist)))

    it = _FakeIntegrator()
    masses = np.array([18.0, 1.0], dtype=np.float64)
    initialize_rpmd_velocities_nm(
        it,
        masses,
        n_beads=6,
        temperature_k=243.0,
        seed=12345,
    )
    assert [c[0] for c in it.calls] == list(range(6))
    assert all(c[1] == 2 for c in it.calls)
