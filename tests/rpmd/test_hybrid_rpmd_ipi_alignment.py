#!/usr/bin/env python3
"""Hybrid RPMD: classical beads follow i-PI frozen-ring / bead-0 force convention."""

from __future__ import annotations

import numpy as np
import pytest

from openmm import Context, HarmonicBondForce, Platform, RPMDIntegrator, System
from openmm import unit

if not hasattr(RPMDIntegrator, "setParticleType"):
    pytest.skip("Extended RPMD API (hybrid particle types) not available", allow_module_level=True)


def _classical_quantum_dimer():
    """Particle 0 classical, particle 1 quantum; single harmonic bond."""
    system = System()
    system.addParticle(12.0)
    system.addParticle(1.0)
    bond = HarmonicBondForce()
    bond.addBond(0, 1, 0.14, 200000.0)
    system.addForce(bond)
    return system


def _make_integrator(num_beads: int):
    integrator = RPMDIntegrator(
        num_beads,
        100.0 * unit.kelvin,
        1.0 / unit.picosecond,
        0.0005 * unit.picosecond,
    )
    integrator.setThermostatType(RPMDIntegrator.NoneThermo)
    integrator.setApplyThermostat(False)
    integrator.setClassicalThermostat(RPMDIntegrator.NoneClassical)
    integrator.setDefaultQuantum(False)
    integrator.setParticleType(0, 0)
    integrator.setParticleType(1, 1)
    integrator.setQuantumParticleTypes({1})
    integrator.setRandomNumberSeed(0)
    return integrator


def _zero_velocities(num_beads: int, n_atoms: int):
    z = unit.Quantity((0.0, 0.0, 0.0), unit.nanometer / unit.picosecond)
    return [z] * n_atoms


def test_hybrid_desynced_classical_quantum_trajectory_matches_bead0_forces():
    """
    Intentionally place the classical atom at different coordinates on a non-zero
    bead than on bead 0. Force evaluation must use bead-0 classical positions for
    every replica (i-PI frozen-ring / HYBRID_RPMD.md), so the quantum atom's
    trajectory matches the fully synced initialization.
    """
    num_beads = 6
    system = _classical_quantum_dimer()
    rng = np.random.RandomState(7)

    def bead_positions():
        """Classical coords identical on all beads; quantum coords vary by bead."""
        out = []
        for bead in range(num_beads):
            classical = (0.10, 0.11, 0.12)
            q = (
                0.14 + 0.008 * bead + 0.002 * rng.randn(),
                0.10 + 0.002 * rng.randn(),
                0.11 + 0.002 * rng.randn(),
            )
            out.append(
                [
                    unit.Quantity(classical, unit.nanometer),
                    unit.Quantity(tuple(q), unit.nanometer),
                ]
            )
        return out

    synced = bead_positions()

    def run(desync_bead: int | None):
        integrator = _make_integrator(num_beads)
        context = Context(system, integrator, Platform.getPlatformByName("Reference"))
        try:
            for b in range(num_beads):
                integrator.setPositions(b, synced[b])
                integrator.setVelocities(b, _zero_velocities(num_beads, 2))
            if desync_bead is not None:
                pos = integrator.getState(desync_bead, getPositions=True).getPositions()
                p0 = pos[0].value_in_unit(unit.nanometer)
                shifted = (p0[0] + 0.08, p0[1] + 0.03, p0[2] - 0.02)
                new_pos = [
                    unit.Quantity(shifted, unit.nanometer),
                    pos[1],
                ]
                integrator.setPositions(desync_bead, new_pos)
            integrator.step(1)
            q_xyz = []
            for b in range(num_beads):
                st = integrator.getState(b, getPositions=True)
                p = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                q_xyz.append(p[1].copy())
            return np.stack(q_xyz, axis=0)
        finally:
            del context, integrator

    ref = run(None)
    trial = run(3)
    np.testing.assert_allclose(trial, ref, rtol=0.0, atol=1e-5)
