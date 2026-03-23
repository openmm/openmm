"""
OpenMM bridge for :mod:`generate_ice_ih` — does not modify ``generate_ice_ih.py``.

Builds ``Topology`` + positions + box from explicit ``nx``, ``ny``, ``nz`` replication.
"""
from __future__ import annotations

from openmm import Vec3, unit
from openmm.app import Element, Topology


def openmm_from_ice_supercell(nx: int, ny: int, nz: int):
    """
    Parameters
    ----------
    nx, ny, nz
        Supercell replication counts (same as ``generate_ice_ih.replicate``).

    Returns
    -------
    topology : Topology
    positions_nm : list[list[float]]
        Cartesian coordinates in nm, O/H/H per residue order.
    box_vectors : tuple
        Orthorhombic periodic box vectors (Vec3 * nanometer).
    """
    if nx < 1 or ny < 1 or nz < 1:
        raise ValueError("nx, ny, nz must be >= 1")

    from generate_ice_ih import build_unit_cell, replicate

    unit_cell = build_unit_cell()
    waters, supercell = replicate(unit_cell, nx, ny, nz)

    topology = Topology()
    positions_nm: list[list[float]] = []
    for O, H1, H2 in waters:
        chain = topology.addChain()
        residue = topology.addResidue("HOH", chain)
        o_atom = topology.addAtom("O", Element.getBySymbol("O"), residue)
        h1_atom = topology.addAtom("H", Element.getBySymbol("H"), residue)
        h2_atom = topology.addAtom("H", Element.getBySymbol("H"), residue)
        topology.addBond(o_atom, h1_atom)
        topology.addBond(o_atom, h2_atom)
        for p in (O, H1, H2):
            positions_nm.append([float(p[0]) / 10.0, float(p[1]) / 10.0, float(p[2]) / 10.0])

    lx_nm = float(supercell[0]) / 10.0
    ly_nm = float(supercell[1]) / 10.0
    lz_nm = float(supercell[2]) / 10.0
    box_vectors = (
        Vec3(lx_nm, 0.0, 0.0) * unit.nanometer,
        Vec3(0.0, ly_nm, 0.0) * unit.nanometer,
        Vec3(0.0, 0.0, lz_nm) * unit.nanometer,
    )
    topology.setPeriodicBoxVectors(box_vectors)
    return topology, positions_nm, box_vectors


def num_molecules_from_supercell(nx: int, ny: int, nz: int) -> int:
    """8 molecules per orthorhombic unit cell."""
    return 8 * nx * ny * nz
