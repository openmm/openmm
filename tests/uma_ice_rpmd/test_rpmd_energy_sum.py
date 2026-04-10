"""Unit tests for RPMD energy helper in test_uma_ice_rpmd."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from openmm import unit


def test_rpmd_sum_bead_pe_plus_ke_kjmol() -> None:
    from test_uma_ice_rpmd import _rpmd_sum_bead_pe_plus_ke_kjmol

    class _E:
        def value_in_unit(self, u):
            if u == unit.kilojoules_per_mole:
                return 1.0
            raise AssertionError(u)

    class _St:
        def getPotentialEnergy(self):
            return _E()

        def getKineticEnergy(self):
            return _E()

    integ = MagicMock()
    integ.getState.side_effect = lambda i, getEnergy=False: _St()

    s = _rpmd_sum_bead_pe_plus_ke_kjmol(integ, 3)
    assert s == pytest.approx(6.0)  # 3 beads * (1 + 1) kJ/mol
