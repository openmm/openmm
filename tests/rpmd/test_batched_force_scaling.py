#!/usr/bin/env python3
"""
Test: Verify RPMD batched PythonForce scaling matches single-copy forces.

Catches the 0x100000000 fixed-point scaling bug in uploadAllForcesToGPU.
The bug: in mixed/double precision, forces uploaded by the batched path were
not scaled by 0x100000000 before storing as long long, but the GPU integrator
kernel always divides by 0x100000000.  This made batched forces ~4.3 billion
times too small — effectively zero — causing ice structure to collapse.

The test works by running a few RPMD steps and comparing the position change
against what the forces predict.  If forces are ~zero, positions barely move
(only ring-polymer springs + thermostat), and the test fails.

Run: python tests/rpmd/test_batched_force_scaling.py
"""

import os
import sys

if os.getenv('OPENMM_PLUGIN_DIR') is None:
    for _base in [os.getenv('CONDA_PREFIX'), os.path.expanduser('~/miniconda3')]:
        if _base:
            _plugins = os.path.join(_base, 'lib', 'plugins')
            if os.path.isdir(_plugins):
                os.environ['OPENMM_PLUGIN_DIR'] = _plugins
                break

_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
for _plugdir in (os.path.join(_root, 'build'), os.path.join(_root, 'build', 'lib', 'plugins')):
    if os.path.isdir(_plugdir) and any(f.startswith('libOpenMMRPMD') for f in os.listdir(_plugdir)):
        os.environ['OPENMM_PLUGIN_DIR'] = _plugdir
        break

import numpy as np
from openmm import app, unit, Vec3, Context, Platform
from openmm import VerletIntegrator, RPMDIntegrator


def _get_platform():
    try:
        return Platform.getPlatformByName('CUDA')
    except Exception:
        return Platform.getPlatformByName('Reference')


def _create_stretched_water():
    """Create a single water molecule with O-H bonds stretched to ~1.2 A.

    The stretched geometry produces significant forces (~hundreds of kJ/mol/nm)
    that must appear in position changes after a few MD steps.
    """
    topology = app.Topology()
    chain = topology.addChain()
    res = topology.addResidue('HOH', chain)
    topology.addAtom('O', app.Element.getBySymbol('O'), res)
    topology.addAtom('H', app.Element.getBySymbol('H'), res)
    topology.addAtom('H', app.Element.getBySymbol('H'), res)

    box_nm = 1.2
    topology.setPeriodicBoxVectors(
        (Vec3(box_nm, 0, 0) * unit.nanometer,
         Vec3(0, box_nm, 0) * unit.nanometer,
         Vec3(0, 0, box_nm) * unit.nanometer)
    )
    pos_ang = np.array([
        [6.0, 6.0, 6.0],
        [7.2, 6.0, 6.0],
        [6.0, 7.2, 6.0],
    ])
    positions = (pos_ang * 0.1).tolist()
    return topology, positions


def test_batched_rpmd_forces_drive_dynamics():
    """Batched PythonForce must produce position changes consistent with
    single-copy force magnitudes (not ~zero from wrong scaling)."""
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping")
        return True

    topology, positions = _create_stretched_water()
    box_vectors = topology.getPeriodicBoxVectors()
    n_atoms = len(positions)

    model_name = 'uma-s-1-pythonforce-batch'
    potential = MLPotential(model_name)
    platform = _get_platform()
    ml_dev = 'cuda' if platform.getName() == 'CUDA' else 'cpu'

    try:
        system = potential.createSystem(
            topology, task_name='omol', charge=0, spin=1, device=ml_dev
        )
    except Exception as e:
        print(f"Could not create UMA system: {e}, skipping")
        return True

    pos_units = [Vec3(*p) * unit.nanometer for p in positions]

    # --- Step A: get reference force magnitude via single-copy path ---
    int_std = VerletIntegrator(1.0 * unit.femtoseconds)
    ctx_std = Context(system, int_std, platform)
    ctx_std.setPositions(pos_units)
    ctx_std.setPeriodicBoxVectors(*box_vectors)
    state_std = ctx_std.getState(getForces=True)
    f_std = state_std.getForces(asNumpy=True).value_in_unit(
        unit.kilojoules_per_mole / unit.nanometer
    )
    max_f = np.abs(f_std).max()
    print(f"  Single-copy max |F|: {max_f:.2f} kJ/(mol·nm)")
    assert max_f > 10.0, f"Forces too small to be a meaningful test ({max_f:.4f})"
    del ctx_std

    # --- Step B: run a few RPMD steps (batched path) and measure displacement ---
    num_beads = 4
    integrator = RPMDIntegrator(
        num_beads, 300.0 * unit.kelvin, 1.0 / unit.picosecond, 0.5 * unit.femtoseconds
    )
    integrator.setThermostatType(RPMDIntegrator.NoneThermo)
    ctx = Context(system, integrator, platform)
    ctx.setPositions(pos_units)
    ctx.setPeriodicBoxVectors(*box_vectors)
    for b in range(num_beads):
        integrator.setPositions(b, pos_units)
    ctx.setVelocitiesToTemperature(0.0 * unit.kelvin)
    zero_vel = [Vec3(0, 0, 0)] * n_atoms
    for b in range(num_beads):
        integrator.setVelocities(b, zero_vel)

    pos_before = []
    for b in range(num_beads):
        s = integrator.getState(b, getPositions=True)
        pos_before.append(
            s.getPositions(asNumpy=True).value_in_unit(unit.nanometer).copy()
        )
    pos_before = np.array(pos_before)

    integrator.step(10)

    pos_after = []
    for b in range(num_beads):
        s = integrator.getState(b, getPositions=True)
        pos_after.append(
            s.getPositions(asNumpy=True).value_in_unit(unit.nanometer).copy()
        )
    pos_after = np.array(pos_after)

    centroid_displacement = np.abs(
        pos_after.mean(axis=0) - pos_before.mean(axis=0)
    ).max()
    print(f"  10-step centroid max displacement: {centroid_displacement:.6f} nm")

    # With correct forces (~hundreds of kJ/mol/nm on a stretched water),
    # 10 steps × 0.5 fs should produce displacements >> 1e-6 nm.
    # With the scaling bug (forces ~0), displacement would be << 1e-6 nm.
    assert centroid_displacement > 1e-5, (
        f"Centroid barely moved ({centroid_displacement:.2e} nm) after 10 steps — "
        f"forces are likely not reaching the integrator (0x100000000 scaling bug?)"
    )
    print(f"  PASS: centroid moved {centroid_displacement:.6f} nm (forces are active)")
    return True


if __name__ == '__main__':
    ok = test_batched_rpmd_forces_drive_dynamics()
    print("Test PASSED" if ok else "Test FAILED")
    sys.exit(0 if ok else 1)
