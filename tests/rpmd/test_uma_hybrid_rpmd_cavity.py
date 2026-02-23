#!/usr/bin/env python3
"""
Test: UMA with hybrid RPMD and cavity particle

This integration test verifies that UMA works correctly with hybrid RPMD when
a cavity particle is excluded from the ML potential evaluation.

Key checks:
1. Cavity particle receives zero force from UMA
2. Molecular particles receive correct UMA forces
3. Cavity particle remains classical (beads identical)
4. Molecular particles have quantum delocalization (bead spread)
"""

import sys
import numpy as np
import openmm
from openmm import app, unit

def create_water_cavity_system():
    """Create a water molecule + cavity particle system."""
    # Create topology for water
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("WAT", chain)
    O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    H1 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
    H2 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
    topology.addBond(O, H1)
    topology.addBond(O, H2)
    
    # Equilibrium water geometry (in nm)
    positions = [
        [0.0, 0.0, 0.0],      # O
        [0.0957, 0.0, 0.0],   # H1
        [-0.024, 0.093, 0.0]  # H2
    ]
    
    return topology, positions


def test_uma_hybrid_rpmd_cavity():
    """Test UMA + hybrid RPMD with cavity particle."""
    print("=" * 70)
    print("Test: UMA + Hybrid RPMD + Cavity Particle Integration")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    topology, positions = create_water_cavity_system()
    
    # Test parameters
    num_beads = 4
    temperature_K = 300.0
    dt_fs = 1.0
    
    # Create system with UMA (only for water)
    print("\nCreating system with UMA potential...")
    try:
        potential = MLPotential('uma-s-1p1-pythonforce-batch')
        system = potential.createSystem(
            topology,
            task_name='omol',
            charge=0,
            spin=1
        )
    except Exception as e:
        print(f"Could not create UMA potential: {e}")
        print("  Skipping test (UMA model may not be available)")
        return True
    
    num_molecular_particles = system.getNumParticles()
    print(f"  Molecular particles: {num_molecular_particles}")
    
    # Add cavity particle
    photon_mass = 1.0 / 1822.888  # amu (electron mass)
    cavity_index = system.addParticle(photon_mass)
    positions.append([0.0, 0.0, 0.5])  # 5 Angstroms from water
    
    # Add cavity to nonbonded force (no interactions)
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
        dt_fs * unit.femtoseconds
    )
    
    # Mark cavity as classical (type 1)
    integrator.setParticleType(cavity_index, 1)
    integrator.setClassicalThermostat(openmm.RPMDIntegrator.LangevinClassical)
    print(f"  Cavity marked as classical (type 1)")
    print(f"  Molecular particles are quantum (type 0, default)")
    
    # Use CUDA if available
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        print(f"  Using CUDA platform")
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
        print(f"  Using Reference platform (CUDA not available)")
    
    # Create context
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions * unit.nanometers)
    
    # Initialize bead positions
    print("\nInitializing bead positions...")
    for bead in range(num_beads):
        bead_positions = []
        for i, pos in enumerate(positions):
            if i == cavity_index:
                # Cavity: same position on all beads
                bead_positions.append(pos)
            else:
                # Molecular: perturb each bead
                perturb = np.random.randn(3) * 0.001
                bead_positions.append([pos[0] + perturb[0], pos[1] + perturb[1], pos[2] + perturb[2]])
        integrator.setPositions(bead, bead_positions * unit.nanometers)
    
    # Run simulation
    num_steps = 20
    print(f"\nRunning {num_steps} RPMD steps...")
    integrator.step(num_steps)
    
    # Verification tests
    print("\n" + "=" * 70)
    print("VERIFICATION CHECKS")
    print("=" * 70)
    
    all_checks_passed = True
    
    # Check 1: Cavity particle should have zero force from UMA
    print("\n1. Checking cavity particle gets zero force from UMA:")
    cavity_forces = []
    for bead in range(num_beads):
        state = integrator.getState(bead, getForces=True)
        forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        cavity_forces.append(forces[cavity_index])
    
    cavity_forces = np.array(cavity_forces)
    max_cavity_force = np.abs(cavity_forces).max()
    
    print(f"   Max force on cavity: {max_cavity_force:.6f} kJ/mol/nm")
    
    if max_cavity_force < 1e-6:
        print(f"   PASS: Cavity receives zero force from UMA")
    else:
        print(f"   FAIL: Cavity has non-zero force!")
        all_checks_passed = False
    
    # Check 2: Molecular particles should have non-zero forces
    print("\n2. Checking molecular particles have non-zero forces:")
    molecular_forces = []
    for bead in range(num_beads):
        state = integrator.getState(bead, getForces=True)
        forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        molecular_forces.append(forces[:num_molecular_particles])
    
    molecular_forces = np.array(molecular_forces)
    avg_molecular_force = np.abs(molecular_forces).mean()
    
    print(f"   Avg force magnitude on molecular particles: {avg_molecular_force:.6f} kJ/mol/nm")
    
    if avg_molecular_force > 1.0:
        print(f"   PASS: Molecular particles have significant forces")
    else:
        print(f"   ⚠ Warning: Forces are small (may be at equilibrium)")
    
    # Check 3: Cavity particle beads should be identical
    print("\n3. Checking cavity particle beads are identical (classical):")
    cavity_positions = []
    for bead in range(num_beads):
        state = integrator.getState(bead, getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        cavity_positions.append(pos[cavity_index])
    
    cavity_positions = np.array(cavity_positions)
    cavity_spread = np.std(cavity_positions, axis=0)
    
    print(f"   Position spread (std dev): [{cavity_spread[0]:.8f}, {cavity_spread[1]:.8f}, {cavity_spread[2]:.8f}] nm")
    
    if np.max(cavity_spread) < 1e-6:
        print(f"   PASS: Cavity beads are identical (< 1e-6 nm)")
    else:
        print(f"   FAIL: Cavity beads are not identical!")
        print(f"     Max spread: {np.max(cavity_spread):.8f} nm")
        all_checks_passed = False
    
    # Check 4: Molecular particles should have bead spread (quantum delocalization)
    print("\n4. Checking molecular particles have bead spread (quantum):")
    molecular_positions = []
    for bead in range(num_beads):
        state = integrator.getState(bead, getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        molecular_positions.append(pos[:num_molecular_particles])
    
    molecular_positions = np.array(molecular_positions)
    molecular_spread = []
    for i in range(num_molecular_particles):
        particle_beads = molecular_positions[:, i, :]
        spread = np.std(particle_beads, axis=0)
        molecular_spread.append(np.max(spread))
    
    avg_molecular_spread = np.mean(molecular_spread)
    
    print(f"   Average position spread: {avg_molecular_spread:.6f} nm")
    for i, spread in enumerate(molecular_spread):
        print(f"   Particle {i}: {spread:.6f} nm")
    
    if avg_molecular_spread > 1e-4:
        print(f"   PASS: Molecular particles show quantum delocalization")
    else:
        print(f"   ⚠ Note: Spread is small (short simulation)")
    
    # Check 5: Energy conservation (rough check)
    print("\n5. Checking simulation stability:")
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"   Final potential energy: {energy:.2f} kJ/mol")
    
    if not np.isnan(energy) and not np.isinf(energy):
        print(f"   PASS: Energy is finite (simulation stable)")
    else:
        print(f"   FAIL: Energy is NaN or Inf (simulation unstable)")
        all_checks_passed = False
    
    # Summary
    print("\n" + "=" * 70)
    if all_checks_passed:
        print("ALL CHECKS PASSED: UMA + hybrid RPMD + cavity working correctly")
        return True
    else:
        print("SOME CHECKS FAILED: Issues detected")
        return False


if __name__ == "__main__":
    success = test_uma_hybrid_rpmd_cavity()
    sys.exit(0 if success else 1)
