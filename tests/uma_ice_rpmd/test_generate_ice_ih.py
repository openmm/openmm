"""Tests for generate_ice_ih (verbatim) and ice_ih_openmm bridge."""
from __future__ import annotations

import numpy as np
import pytest

from generate_ice_ih import build_unit_cell, replicate
from ice_ih_openmm import num_molecules_from_supercell, openmm_from_ice_supercell


def test_unit_cell_eight_molecules() -> None:
    waters = build_unit_cell()
    assert len(waters) == 8
    for O, H1, H2 in waters:
        d1 = np.linalg.norm(H1 - O)
        d2 = np.linalg.norm(H2 - O)
        assert d1 == pytest.approx(0.9572, rel=1e-3)
        assert d2 == pytest.approx(0.9572, rel=1e-3)


def test_replicate_1x2x2_is_32_molecules() -> None:
    uc = build_unit_cell()
    waters, supercell = replicate(uc, 1, 2, 2)
    assert len(waters) == 32
    assert num_molecules_from_supercell(1, 2, 2) == 32


def test_openmm_from_supercell_1_2_2() -> None:
    from openmm import unit

    topo, pos, box = openmm_from_ice_supercell(1, 2, 2)
    assert topo.getNumAtoms() == 96
    assert len(pos) == 96
    m = np.array([[box[i][j].value_in_unit(unit.nanometer) for j in range(3)] for i in range(3)])
    assert np.abs(np.linalg.det(m)) > 1e-3


def test_openmm_rejects_invalid_replication() -> None:
    with pytest.raises(ValueError):
        openmm_from_ice_supercell(0, 1, 1)
