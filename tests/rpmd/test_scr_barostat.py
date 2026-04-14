#!/usr/bin/env python3
"""RPMD stochastic cell rescaling barostat (SCR)."""

import math

import numpy as np
import pytest

from openmm import Context, HarmonicBondForce, NonbondedForce, Platform, RPMDIntegrator, System
from openmm import unit

try:
    from openmm import RPMDStochasticCellRescalingBarostat as RPMDStochasticCellRescalingBarostat
except ImportError:
    RPMDStochasticCellRescalingBarostat = None  # type: ignore

if RPMDStochasticCellRescalingBarostat is None:
    pytest.skip("RPMDStochasticCellRescalingBarostat not available", allow_module_level=True)


def _argon_pair_periodic():
    """Two Lennard-Jones atoms in a periodic box."""
    system = System()
    for _ in range(2):
        system.addParticle(39.9)
    nb = NonbondedForce()
    sigma = 0.34
    eps = 0.5
    nb.addParticle(0.0, sigma, eps)
    nb.addParticle(0.0, sigma, eps)
    nb.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(0.9)
    system.addForce(nb)
    L = 2.2
    system.setDefaultPeriodicBoxVectors(
        unit.Quantity([L, 0, 0], unit.nanometer),
        unit.Quantity([0, L, 0], unit.nanometer),
        unit.Quantity([0, 0, L], unit.nanometer),
    )
    return system


def test_scr_runs_finite_volume():
    """SCR barostat steps run; box volume stays positive and finite."""
    B = RPMDStochasticCellRescalingBarostat
    system = _argon_pair_periodic()
    baro = B(1.0 * unit.bar, 0.5 * unit.picosecond, 2)
    baro.setDefaultPressure(1.0 * unit.bar)
    baro.setRandomNumberSeed(123)
    system.addForce(baro)

    integrator = RPMDIntegrator(
        4,
        120.0 * unit.kelvin,
        5.0 / unit.picosecond,
        0.001 * unit.picosecond,
    )
    integrator.setThermostatType(RPMDIntegrator.Pile)

    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    rng = np.random.RandomState(0)
    for bead in range(4):
        pos = [
            unit.Quantity((0.5 + 0.05 * rng.randn(), 0.5, 0.5), unit.nanometer),
            unit.Quantity((0.9 + 0.05 * rng.randn(), 0.5, 0.5), unit.nanometer),
        ]
        integrator.setPositions(bead, pos)
        integrator.setVelocities(
            bead,
            [
                unit.Quantity(tuple(rng.randn(3)), unit.nanometer / unit.picosecond),
                unit.Quantity(tuple(rng.randn(3)), unit.nanometer / unit.picosecond),
            ],
        )

    st0 = context.getState()
    v0 = st0.getPeriodicBoxVectors()
    vol0 = v0[0][0] * v0[1][1] * v0[2][2]
    integrator.step(400)
    st1 = context.getState()
    v1 = st1.getPeriodicBoxVectors()
    vol1 = v1[0][0] * v1[1][1] * v1[2][2]
    assert vol0 > 0 and vol1 > 0
    assert math.isfinite(vol1)
    assert abs(vol1 - vol0) / vol0 < 0.5


def test_scr_pressure_parameter():
    """Context exposes SCR pressure parameter."""
    B = RPMDStochasticCellRescalingBarostat
    system = _argon_pair_periodic()
    system.addForce(B(2.0 * unit.bar, 1.0 * unit.picosecond, 10))
    integrator = RPMDIntegrator(
        2,
        200.0 * unit.kelvin,
        5.0 / unit.picosecond,
        0.002 * unit.picosecond,
    )
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    for bead in range(2):
        integrator.setPositions(
            bead,
            [
                unit.Quantity((0.5, 0.5, 0.5), unit.nanometer),
                unit.Quantity((0.85, 0.5, 0.5), unit.nanometer),
            ],
        )
    p = context.getParameter(B.Pressure())
    assert math.isfinite(p)
