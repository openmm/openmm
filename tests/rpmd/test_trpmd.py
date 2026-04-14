#!/usr/bin/env python3
"""TRPMD (internal modes thermostatted, centroid NVE) and Fast-Forward Langevin."""

import math

import numpy as np
import pytest

from openmm import Context, HarmonicBondForce, Platform, RPMDIntegrator, System
from openmm import unit

if not hasattr(RPMDIntegrator, "Trpmd"):
    pytest.skip("Extended RPMD API (TRPMD/FFL) not available", allow_module_level=True)


def _dimer_system():
    system = System()
    for _ in range(2):
        system.addParticle(1.0)
    bond = HarmonicBondForce()
    bond.addBond(0, 1, 0.15, 120000.0)
    system.addForce(bond)
    return system


def test_trpmd_enum_and_runs():
    """TRPMD thermostat type runs without error on Reference."""
    assert hasattr(RPMDIntegrator, "Trpmd")
    system = _dimer_system()
    integrator = RPMDIntegrator(
        8,
        250.0 * unit.kelvin,
        4.0 / unit.picosecond,
        0.0005 * unit.picosecond,
    )
    integrator.setThermostatType(RPMDIntegrator.Trpmd)

    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    rng = np.random.RandomState(5)
    for bead in range(8):
        pos = [
            unit.Quantity((0.1 + 0.01 * rng.randn(), 0.0, 0.0), unit.nanometer),
            unit.Quantity((0.25 + 0.01 * rng.randn(), 0.0, 0.0), unit.nanometer),
        ]
        integrator.setPositions(bead, pos)
        integrator.setVelocities(
            bead,
            [
                unit.Quantity(tuple(rng.randn(3)), unit.nanometer / unit.picosecond),
                unit.Quantity(tuple(rng.randn(3)), unit.nanometer / unit.picosecond),
            ],
        )
    integrator.step(200)
    assert math.isfinite(integrator.getTotalEnergy())


def test_fast_forward_langevin_runs():
    """FastForwardLangevin thermostat runs and keeps finite energy."""
    assert hasattr(RPMDIntegrator, "FastForwardLangevin")
    system = _dimer_system()
    integrator = RPMDIntegrator(
        6,
        300.0 * unit.kelvin,
        20.0 / unit.picosecond,
        0.0004 * unit.picosecond,
    )
    integrator.setThermostatType(RPMDIntegrator.FastForwardLangevin)

    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    rng = np.random.RandomState(9)
    for bead in range(6):
        pos = [
            unit.Quantity((0.1, 0.0, 0.0), unit.nanometer),
            unit.Quantity((0.26, 0.0, 0.0), unit.nanometer),
        ]
        integrator.setPositions(bead, pos)
        integrator.setVelocities(
            bead,
            [
                unit.Quantity(tuple(2.0 * rng.randn(3)), unit.nanometer / unit.picosecond),
                unit.Quantity(tuple(2.0 * rng.randn(3)), unit.nanometer / unit.picosecond),
            ],
        )
    integrator.step(100)
    assert math.isfinite(integrator.getTotalEnergy())
