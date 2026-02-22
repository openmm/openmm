#!/usr/bin/env python3
"""
Test: Hybrid RPMD with cavity particle excluded from RPMD.

Verifies:
1. Cavity particle beads stay identical (classical behavior)
2. Molecular particles develop bead spread (quantum delocalization)
3. Simulation runs without errors on CUDA
4. Energy is physically reasonable
"""

import sys
import numpy as np

try:
    from openmm import openmm, unit
except ImportError:
    import openmm
    from openmm import unit

# Constants
BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5


def create_small_dimer_system(num_molecules=5, box_size_nm=1.5):
    """Create a small dimer system for testing."""
    np.random.seed(42)
    system = openmm.System()
    positions = []

    mass_O = 16.0  # amu
    k_OO_au = 0.73204
    r0_OO_au = 2.281655158
    k_OO = k_OO_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
    r0_OO = r0_OO_au * BOHR_TO_NM

    charge_magnitude = 0.3
    sigma_O = 6.230426584 * BOHR_TO_NM
    epsilon_O = 0.00016685201 * HARTREE_TO_KJMOL

    bond_force = openmm.HarmonicBondForce()
    nonbonded_force = openmm.NonbondedForce()
    nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nonbonded_force.setCutoffDistance(0.7)

    side = int(np.ceil(num_molecules**(1/3)))
    spacing = box_size_nm / side

    mol_idx = 0
    for i in range(side):
        for j in range(side):
            for kk in range(side):
                if mol_idx >= num_molecules:
                    break
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing
                z = (kk + 0.5) * spacing
                theta = np.random.rand() * 2 * np.pi
                phi = np.arccos(2 * np.random.rand() - 1)
                d = np.array([np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)])
                r1 = np.array([x, y, z]) - 0.5 * r0_OO * d
                r2 = np.array([x, y, z]) + 0.5 * r0_OO * d
                idx1 = system.addParticle(mass_O)
                idx2 = system.addParticle(mass_O)
                positions.append(openmm.Vec3(*r1) * unit.nanometer)
                positions.append(openmm.Vec3(*r2) * unit.nanometer)
                bond_force.addBond(idx1, idx2, r0_OO, k_OO)
                nonbonded_force.addParticle(-charge_magnitude, sigma_O, epsilon_O)
                nonbonded_force.addParticle(+charge_magnitude, sigma_O, epsilon_O)
                nonbonded_force.addException(idx1, idx2, 0.0, 1.0, 0.0)
                mol_idx += 1
            if mol_idx >= num_molecules:
                break
        if mol_idx >= num_molecules:
            break

    system.addForce(bond_force)
    system.addForce(nonbonded_force)
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_size_nm, 0, 0),
        openmm.Vec3(0, box_size_nm, 0),
        openmm.Vec3(0, 0, box_size_nm)
    )
    return system, positions


def test_hybrid_rpmd_cavity():
    """Test that cavity particle is treated classically in hybrid RPMD."""
    print("=" * 60)
    print("Test: Hybrid RPMD with Classical Cavity Particle")
    print("=" * 60)

    num_molecules = 5
    num_beads = 4
    temperature_K = 100.0
    dt_ps = 0.001  # 1 fs
    num_steps = 500

    # Create system
    system, positions = create_small_dimer_system(num_molecules)
    num_molecular_particles = system.getNumParticles()
    print(f"  Molecular particles: {num_molecular_particles}")

    # Add cavity particle
    photon_mass = 1.0 / 1822.888  # amu (electron mass)
    cavity_index = system.addParticle(photon_mass)
    positions.append(openmm.Vec3(0.01, 0, 0) * unit.nanometer)

    # Add cavity to nonbonded (no interactions)
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.addParticle(0.0, 0.1, 0.0)

    print(f"  Cavity particle index: {cavity_index}")
    print(f"  Total particles: {system.getNumParticles()}")

    # Create RPMD integrator
    integrator = openmm.RPMDIntegrator(
        num_beads,
        temperature_K * unit.kelvin,
        1.0 / unit.picosecond,
        dt_ps * unit.picosecond
    )
    integrator.setThermostatType(openmm.RPMDIntegrator.PileG)
    integrator.setCentroidFriction(0.5 / unit.picosecond)

    # Mark cavity as classical (type 1 -- not in quantum types)
    integrator.setParticleType(cavity_index, 1)
    integrator.setClassicalThermostat(openmm.RPMDIntegrator.LangevinClassical)
    print(f"  Cavity marked as classical (type 1)")

    # Select platform
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        print(f"  Platform: CUDA")
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
        print(f"  Platform: Reference (CUDA not available)")

    # Create context
    context = openmm.Context(system, integrator, platform)

    # Set positions for all beads
    for bead in range(num_beads):
        bead_positions = []
        for i, pos in enumerate(positions):
            pos_nm = pos.value_in_unit(unit.nanometer)
            if i == cavity_index:
                # Same position for cavity on all beads
                bead_positions.append(openmm.Vec3(0.01, 0, 0) * unit.nanometer)
            else:
                dx = np.random.randn() * 0.005
                dy = np.random.randn() * 0.005
                dz = np.random.randn() * 0.005
                bead_positions.append(
                    openmm.Vec3(pos_nm[0]+dx, pos_nm[1]+dy, pos_nm[2]+dz) * unit.nanometer
                )
        integrator.setPositions(bead, bead_positions)

    context.setPositions(positions)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)

    # Set velocities for all beads
    state = context.getState(getVelocities=True)
    base_vel = state.getVelocities()
    for bead in range(num_beads):
        bead_vel = []
        for v in base_vel:
            vx = v[0].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01*np.random.randn()
            vy = v[1].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01*np.random.randn()
            vz = v[2].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01*np.random.randn()
            bead_vel.append(openmm.Vec3(vx, vy, vz) * unit.nanometers/unit.picoseconds)
        integrator.setVelocities(bead, bead_vel)

    print(f"\n--- Running {num_steps} steps ---")
    integrator.step(num_steps)
    print(f"  Completed {num_steps} steps successfully")

    # Verify: get positions from all beads
    print(f"\n--- Verifying Bead Behavior ---")

    all_cavity_pos = []
    all_mol_pos = {i: [] for i in range(min(4, num_molecular_particles))}

    for bead in range(num_beads):
        state = integrator.getState(bead, getPositions=True)
        pos = state.getPositions(asNumpy=True)
        pos_nm = pos.value_in_unit(unit.nanometer)

        # Cavity position
        all_cavity_pos.append(pos_nm[cavity_index])

        # First few molecular particles
        for i in all_mol_pos:
            all_mol_pos[i].append(pos_nm[i])

    all_cavity_pos = np.array(all_cavity_pos)
    
    # Check cavity beads are identical
    cavity_spread = np.std(all_cavity_pos, axis=0)
    max_cavity_spread = np.max(cavity_spread)
    print(f"  Cavity bead spread (std): {cavity_spread}")
    print(f"  Max cavity bead spread: {max_cavity_spread:.2e} nm")

    # Check molecular particles have bead spread
    mol_spreads = []
    for i, positions_list in all_mol_pos.items():
        pos_arr = np.array(positions_list)
        spread = np.std(pos_arr, axis=0)
        mol_spreads.append(np.max(spread))

    avg_mol_spread = np.mean(mol_spreads)
    print(f"  Average molecular bead spread: {avg_mol_spread:.6f} nm")

    # Assertions
    passed = True

    # Cavity beads should be identical (spread < tiny tolerance)
    if max_cavity_spread > 1e-6:
        print(f"  FAIL: Cavity bead spread too large: {max_cavity_spread:.2e}")
        passed = False
    else:
        print(f"  PASS: Cavity beads are identical (spread < 1e-6)")

    # Molecular beads should have some spread (quantum delocalization)
    if avg_mol_spread < 1e-6:
        print(f"  FAIL: Molecular beads have no spread (no quantum delocalization)")
        passed = False
    else:
        print(f"  PASS: Molecular beads show quantum delocalization")

    # Check energy is finite
    energy = integrator.getTotalEnergy()
    if hasattr(energy, 'value_in_unit'):
        energy_val = energy.value_in_unit(unit.kilojoule_per_mole)
    else:
        energy_val = float(energy)

    if not np.isfinite(energy_val):
        print(f"  FAIL: Energy is not finite: {energy_val}")
        passed = False
    else:
        print(f"  PASS: Energy is finite: {energy_val:.2f} kJ/mol")

    print(f"\n{'='*60}")
    if passed:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print(f"{'='*60}")
    return passed


if __name__ == "__main__":
    success = test_hybrid_rpmd_cavity()
    sys.exit(0 if success else 1)
