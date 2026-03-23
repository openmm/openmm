"""
Unit tests: OpenMM RPMD PILE_G thermostat wiring (Bussi centroid + PILE internal).

See ``RPMDIntegrator::ThermostatType::PileG`` in plugins/rpmd and
``test_uma_ice_rpmd.run_simulation`` defaults.
"""
from __future__ import annotations

import pytest

pytest.importorskip("openmm")
from openmm import RPMDIntegrator, unit


def test_rpmd_integrator_pile_g_and_centroid_friction() -> None:
    """PILE_G is available; centroid friction is set without error."""
    friction = 1.0
    gamma_c = 0.5
    integ = RPMDIntegrator(
        4,
        243.0 * unit.kelvin,
        friction / unit.picosecond,
        0.1 * unit.femtoseconds,
    )
    integ.setThermostatType(RPMDIntegrator.PileG)
    integ.setCentroidFriction(gamma_c / unit.picosecond)
    assert integ.getThermostatType() == RPMDIntegrator.PileG
