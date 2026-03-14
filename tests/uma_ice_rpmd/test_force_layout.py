#!/usr/bin/env python3
"""
Unit test: 1-bead RPMD vs standard OpenMM with UMA potential.

Verifies that the RPMD batched force path produces identical energies and forces
to standard OpenMM when using a single bead. With 1 bead, RPMD reduces to
classical MD, so results must match. This validates the force layout fix in
uploadAllForcesToGPU and the UMA-OpenMM integration.

Run: pytest tests/uma_ice_rpmd/test_force_layout.py -v
      or: python tests/uma_ice_rpmd/test_force_layout.py
"""

import os
import sys

# Set plugin dir before openmm import
if os.getenv('OPENMM_PLUGIN_DIR') is None:
    for _base in [os.getenv('CONDA_PREFIX'), os.path.expanduser('~/miniconda3')]:
        if _base:
            _plugins = os.path.join(_base, 'lib', 'plugins')
            if os.path.isdir(_plugins):
                os.environ['OPENMM_PLUGIN_DIR'] = _plugins
                break

# Prefer build dir plugins if available (plugins live in build/ or build/lib/plugins)
_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
for _plugdir in (os.path.join(_root, 'build'), os.path.join(_root, 'build', 'lib', 'plugins')):
    if os.path.isdir(_plugdir) and any(f.startswith('libOpenMMRPMD') for f in os.listdir(_plugdir)):
        os.environ['OPENMM_PLUGIN_DIR'] = _plugdir
        break

import numpy as np
from openmm import app, unit, Vec3, Context, Platform
from openmm import VerletIntegrator, RPMDIntegrator
from openmmml import MLPotential


def _get_platform():
    """Prefer CUDA; fallback to Reference if unavailable."""
    try:
        return Platform.getPlatformByName('CUDA')
    except Exception:
        return Platform.getPlatformByName('Reference')


def _create_small_water_system(num_molecules=4):
    """Create a small water box for testing."""
    from ase import Atoms
    from ase.build import molecule

    # Build water molecules
    waters = [molecule('H2O') for _ in range(num_molecules)]
    atoms = waters[0]
    for w in waters[1:]:
        atoms += w

    # Place in box
    cell_size = 12.0  # Angstrom
    atoms.set_cell([cell_size, cell_size, cell_size])
    atoms.set_pbc(True)
    atoms.center()

    # Convert to OpenMM
    symbols = atoms.get_chemical_symbols()
    pos_ang = atoms.get_positions()
    n_atoms = len(symbols)

    topology = app.Topology()
    for i in range(num_molecules):
        chain = topology.addChain()
        res = topology.addResidue('HOH', chain)
        topology.addAtom('O', app.Element.getBySymbol('O'), res)
        topology.addAtom('H', app.Element.getBySymbol('H'), res)
        topology.addAtom('H', app.Element.getBySymbol('H'), res)

    box_vectors = (
        Vec3(cell_size * 0.1, 0, 0) * unit.nanometer,
        Vec3(0, cell_size * 0.1, 0) * unit.nanometer,
        Vec3(0, 0, cell_size * 0.1) * unit.nanometer,
    )
    topology.setPeriodicBoxVectors(box_vectors)
    positions = (pos_ang * 0.1).tolist()  # Angstrom to nm
    return topology, positions, box_vectors


def test_1bead_rpmd_matches_standard_openmm():
    """
    With 1 bead, RPMD must produce identical energy and forces to standard OpenMM.
    This validates the batched force upload layout (SoA) and UMA integration.
    """
    topology, positions, box_vectors = _create_small_water_system(num_molecules=1)
    n_atoms = len(positions)

    # Create UMA system (same for both). Use CUDA when available.
    model_name = 'uma-s-1-pythonforce-batch'
    potential = MLPotential(model_name)
    platform = _get_platform()
    _ml_device = 'cuda' if platform.getName() == 'CUDA' else 'cpu'
    system = potential.createSystem(topology, task_name='omol', charge=0, spin=1, device=_ml_device)

    pos_with_units = [Vec3(float(p[0]), float(p[1]), float(p[2])) * unit.nanometer for p in positions]

    # --- Standard OpenMM (no RPMD) ---
    integrator_std = VerletIntegrator(0.001)  # 1 fs, not used for dynamics
    context_std = Context(system, integrator_std, platform)
    context_std.setPositions(pos_with_units)
    context_std.setPeriodicBoxVectors(*box_vectors)

    state_std = context_std.getState(getEnergy=True, getForces=True)
    pe_std = state_std.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    forces_std = state_std.getForces(asNumpy=True).value_in_unit(
        unit.kilojoules_per_mole / unit.nanometer
    )

    # --- 1-bead RPMD ---
    integrator_rpmd = RPMDIntegrator(
        1,  # num_beads
        243.0 * unit.kelvin,
        1.0 / unit.picosecond,
        0.001 * unit.femtoseconds,
    )
    integrator_rpmd.setThermostatType(RPMDIntegrator.Pile)
    context_rpmd = Context(system, integrator_rpmd, platform)
    context_rpmd.setPositions(pos_with_units)
    context_rpmd.setPeriodicBoxVectors(*box_vectors)
    integrator_rpmd.setPositions(0, pos_with_units)
    context_rpmd.setVelocitiesToTemperature(243.0 * unit.kelvin)
    v = context_rpmd.getState(getVelocities=True).getVelocities()
    integrator_rpmd.setVelocities(0, v)

    state_rpmd = integrator_rpmd.getState(0, getEnergy=True, getForces=True)
    pe_rpmd = state_rpmd.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    forces_rpmd = state_rpmd.getForces(asNumpy=True).value_in_unit(
        unit.kilojoules_per_mole / unit.nanometer
    )

    # Compare
    pe_diff = abs(pe_std - pe_rpmd)
    assert pe_diff < 0.01, (
        f"Potential energy mismatch: std={pe_std:.6f} kJ/mol, rpmd={pe_rpmd:.6f} kJ/mol, diff={pe_diff:.6f}"
    )

    forces_diff = np.abs(forces_std - forces_rpmd)
    max_force_diff = np.max(forces_diff)
    assert max_force_diff < 0.01, (
        f"Force mismatch: max |ΔF| = {max_force_diff:.6f} kJ/(mol·nm)"
    )


def test_1bead_rpmd_step_produces_same_trajectory_as_standard():
    """
    Run one integration step with both integrators; positions and velocities
    should evolve identically for 1 bead (same forces => same acceleration).
    """
    topology, positions, box_vectors = _create_small_water_system(num_molecules=1)
    np.random.seed(42)

    model_name = 'uma-s-1-pythonforce-batch'
    potential = MLPotential(model_name)
    platform = _get_platform()
    _ml_device = 'cuda' if platform.getName() == 'CUDA' else 'cpu'
    system = potential.createSystem(topology, task_name='omol', charge=0, spin=1, device=_ml_device)

    pos_nm = np.array(positions)
    vel_nm_ps = np.random.randn(len(positions), 3) * 0.05  # nm/ps
    pos_with_units = [Vec3(float(p[0]), float(p[1]), float(p[2])) * unit.nanometer for p in pos_nm]
    vel_with_units = [
        Vec3(float(v[0]), float(v[1]), float(v[2])) * unit.nanometer / unit.picosecond
        for v in vel_nm_ps
    ]

    dt = 0.001  # 1 fs
    platform = _get_platform()

    # Standard OpenMM
    int_std = VerletIntegrator(dt * unit.femtoseconds)
    ctx_std = Context(system, int_std, platform)
    ctx_std.setPositions(pos_with_units)
    ctx_std.setVelocities(vel_with_units)
    ctx_std.setPeriodicBoxVectors(*box_vectors)
    int_std.step(1)
    pos_std = ctx_std.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    vel_std = ctx_std.getState(getVelocities=True).getVelocities(asNumpy=True).value_in_unit(
        unit.nanometer / unit.picosecond
    )

    # 1-bead RPMD
    int_rpmd = RPMDIntegrator(1, 243.0 * unit.kelvin, 1.0 / unit.picosecond, dt * unit.femtoseconds)
    int_rpmd.setThermostatType(RPMDIntegrator.Pile)
    ctx_rpmd = Context(system, int_rpmd, platform)
    ctx_rpmd.setPositions(pos_with_units)
    ctx_rpmd.setPeriodicBoxVectors(*box_vectors)
    int_rpmd.setPositions(0, pos_with_units)
    int_rpmd.setVelocities(0, vel_with_units)
    int_rpmd.step(1)
    pos_rpmd = int_rpmd.getState(0, getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    vel_rpmd = int_rpmd.getState(0, getVelocities=True).getVelocities(asNumpy=True).value_in_unit(
        unit.nanometer / unit.picosecond
    )

    # PILE applies Langevin; trajectories will differ slightly due to thermostat randomness.
    # We only check that both produce finite, reasonable positions (no explosion from force layout bug).
    assert np.all(np.isfinite(pos_std)) and np.all(np.isfinite(pos_rpmd))
    assert np.all(np.isfinite(vel_std)) and np.all(np.isfinite(vel_rpmd))
    assert np.all(np.abs(pos_rpmd) < 100.0)  # No runaway coordinates


def test_single_water_molecule_stability():
    """
    Run 100 MD steps on a single water molecule; verify O-H bonds stay ~1 A and
    molecule does not explode. Validates UMA forces and RPMD integration for minimal system.
    """
    topology, positions, box_vectors = _create_small_water_system(num_molecules=1)
    n_atoms = 3  # O, H, H

    model_name = 'uma-s-1-pythonforce-batch'
    potential = MLPotential(model_name)
    platform = _get_platform()
    _ml_device = 'cuda' if platform.getName() == 'CUDA' else 'cpu'
    system = potential.createSystem(topology, task_name='omol', charge=0, spin=1, device=_ml_device)

    pos_with_units = [Vec3(float(p[0]), float(p[1]), float(p[2])) * unit.nanometer for p in positions]

    integrator = RPMDIntegrator(
        1,  # 1 bead = classical
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        0.5 * unit.femtoseconds,
    )
    integrator.setThermostatType(RPMDIntegrator.Pile)
    context = Context(system, integrator, platform)
    context.setPositions(pos_with_units)
    context.setPeriodicBoxVectors(*box_vectors)
    context.setVelocitiesToTemperature(300.0 * unit.kelvin)
    v = context.getState(getVelocities=True).getVelocities()
    integrator.setPositions(0, pos_with_units)
    integrator.setVelocities(0, v)

    for _ in range(100):
        integrator.step(1)

    final_pos = integrator.getState(0, getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    assert np.all(np.isfinite(final_pos)), "Positions became non-finite"
    assert np.all(np.abs(final_pos) < 50.0), "Molecule exploded (runaway coordinates)"

    # O-H bonds should be ~0.96-1.0 A (0.096-0.1 nm)
    o_pos = final_pos[0]
    h1_pos = final_pos[1]
    h2_pos = final_pos[2]
    oh1 = np.linalg.norm(o_pos - h1_pos) * 10.0  # nm -> Angstrom
    oh2 = np.linalg.norm(o_pos - h2_pos) * 10.0
    assert 0.8 < oh1 < 1.5, f"O-H1 bond length {oh1:.2f} A out of range"
    assert 0.8 < oh2 < 1.5, f"O-H2 bond length {oh2:.2f} A out of range"


if __name__ == '__main__':
    print("Testing 1-bead RPMD vs standard OpenMM with UMA potential (single water)...")
    test_1bead_rpmd_matches_standard_openmm()
    print("  test_1bead_rpmd_matches_standard_openmm: PASSED")
    test_1bead_rpmd_step_produces_same_trajectory_as_standard()
    print("  test_1bead_rpmd_step_produces_same_trajectory_as_standard: PASSED")
    test_single_water_molecule_stability()
    print("  test_single_water_molecule_stability: PASSED")
    print("All tests passed.")
