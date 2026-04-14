#!/usr/bin/env python3
"""Suzuki–Chin finite-difference correction on RPMD."""

import math

import numpy as np
import pytest

from openmm import Context, HarmonicBondForce, Platform, RPMDIntegrator, System
from openmm import unit

if not hasattr(RPMDIntegrator, "setSuzukiChinEnabled"):
    pytest.skip("Extended RPMD API (Suzuki-Chin) not available", allow_module_level=True)


def _dimer():
    system = System()
    for _ in range(2):
        system.addParticle(1.0)
    bond = HarmonicBondForce()
    bond.addBond(0, 1, 0.15, 150000.0)
    system.addForce(bond)
    return system


def test_suzuki_chin_forces_finite():
    """With SC enabled, short dynamics produces finite energies and no NaNs."""
    system = _dimer()
    integrator = RPMDIntegrator(
        4,
        200.0 * unit.kelvin,
        5.0 / unit.picosecond,
        0.0005 * unit.picosecond,
    )
    integrator.setThermostatType(RPMDIntegrator.Pile)
    integrator.setSuzukiChinEnabled(True)
    integrator.setSuzukiChinEpsilon(0.001 * unit.nanometer)

    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    rng = np.random.RandomState(11)
    for bead in range(4):
        dx = 0.01 * rng.randn(2, 3)
        pos = [
            unit.Quantity(tuple(0.1 + dx[0]), unit.nanometer),
            unit.Quantity(tuple(0.25 + dx[1]), unit.nanometer),
        ]
        integrator.setPositions(bead, pos)
        integrator.setVelocities(
            bead,
            [
                unit.Quantity(tuple(rng.randn(3)), unit.nanometer / unit.picosecond),
                unit.Quantity(tuple(rng.randn(3)), unit.nanometer / unit.picosecond),
            ],
        )

    for _ in range(20):
        integrator.step(1)
        e = integrator.getTotalEnergy()
        assert math.isfinite(e)
        st = integrator.getState(0, getForces=True)
        f = st.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        assert np.all(np.isfinite(f))


def test_suzuki_chin_disabled_matches_toggle():
    """Disabling SC should not raise and keeps finite energy."""
    system = _dimer()
    integrator = RPMDIntegrator(
        4,
        250.0 * unit.kelvin,
        5.0 / unit.picosecond,
        0.0005 * unit.picosecond,
    )
    integrator.setSuzukiChinEnabled(False)

    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    rng = np.random.RandomState(2)
    for bead in range(4):
        pos = [
            unit.Quantity((0.1, 0.0, 0.0), unit.nanometer),
            unit.Quantity((0.26, 0.0, 0.0), unit.nanometer),
        ]
        integrator.setPositions(bead, pos)
    integrator.step(15)
    assert math.isfinite(integrator.getTotalEnergy())
