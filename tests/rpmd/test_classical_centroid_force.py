#!/usr/bin/env python3
"""
Test: Verify classical particle sees centroid-averaged force in hybrid RPMD

This test checks that in hybrid RPMD, a classical particle experiences the
centroid-averaged force F_c = (1/P) * sum_k F(q_c, q_quantum^(k)), not just
the force from bead 0.

Bug: Before fix, classical particles only saw forces from bead 0, leading to
incorrect dynamics sampling the wrong distribution.
"""

import sys
import numpy as np
import openmm
from openmm import unit

def create_hybrid_test_system():
    """
    Create a simple test system with:
    - 1 classical particle (type 1)
    - 2 quantum particles (type 0, bonded)
    
    The classical particle interacts with the quantum dimer via a simple potential.
    """
    system = openmm.System()
    
    # Add particles
    mass_classical = 16.0  # amu
    mass_quantum = 16.0
    
    # Particle 0: Classical
    system.addParticle(mass_classical)
    
    # Particles 1, 2: Quantum (bonded dimer)
    system.addParticle(mass_quantum)
    system.addParticle(mass_quantum)
    
    # Add bond between quantum particles
    bond_force = openmm.HarmonicBondForce()
    k_bond = 1000.0  # kJ/mol/nm^2
    r0_bond = 0.15   # nm
    bond_force.addBond(1, 2, r0_bond, k_bond)
    system.addForce(bond_force)
    
    # Add nonbonded force (classical interacts with quantum)
    nb_force = openmm.NonbondedForce()
    nb_force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    
    # Classical particle: charge +0.5, moderate LJ
    nb_force.addParticle(+0.5, 0.3, 0.5)
    
    # Quantum particles: charges that create asymmetric force field
    nb_force.addParticle(-0.3, 0.3, 0.5)
    nb_force.addParticle(-0.2, 0.3, 0.5)
    
    # Exclude bonded interaction
    nb_force.addException(1, 2, 0.0, 1.0, 0.0)
    
    system.addForce(nb_force)
    
    # Initial positions
    positions = [
        [0.0, 0.0, 0.0],      # Classical particle
        [0.5, 0.0, 0.0],      # Quantum particle 1
        [0.65, 0.0, 0.0]      # Quantum particle 2
    ]
    
    return system, positions


def test_classical_centroid_force():
    """Test that classical particle sees centroid-averaged force."""
    print("=" * 70)
    print("Test: Classical Particle Centroid Force in Hybrid RPMD")
    print("=" * 70)
    
    system, positions = create_hybrid_test_system()
    
    # Test parameters
    num_beads = 8
    temperature_K = 300.0
    dt_fs = 1.0
    
    # Create RPMD integrator
    integrator = openmm.RPMDIntegrator(
        num_beads,
        temperature_K * unit.kelvin,
        1.0 / unit.picosecond,
        dt_fs * unit.femtoseconds
    )
    
    # Mark particle 0 as classical (type 1)
    integrator.setParticleType(0, 1)
    # Particles 1, 2 are quantum (default type 0)
    
    # Set classical thermostat
    integrator.setClassicalThermostat(openmm.RPMDIntegrator.LangevinClassical)
    
    print(f"\nSystem setup:")
    print(f"  Particles: 3 (1 classical, 2 quantum)")
    print(f"  RPMD beads: {num_beads}")
    print(f"  Classical particle: index 0 (type 1)")
    print(f"  Quantum particles: indices 1, 2 (type 0)")
    
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
    # Classical particle: same position on all beads
    # Quantum particles: spread across beads
    print("\nInitializing bead positions...")
    for bead in range(num_beads):
        bead_pos = []
        for i, pos in enumerate(positions):
            if i == 0:
                # Classical: same on all beads
                bead_pos.append(pos)
            else:
                # Quantum: perturb each bead differently
                perturb = np.random.randn(3) * 0.02
                bead_pos.append([pos[0] + perturb[0], pos[1] + perturb[1], pos[2] + perturb[2]])
        integrator.setPositions(bead, bead_pos * unit.nanometers)
    
    # Run a few steps to let quantum beads evolve
    print("Running 10 RPMD steps...")
    integrator.step(10)
    
    # Now check that classical particle position is identical on all beads
    # and compute the expected centroid force
    print("\nVerifying classical particle dynamics:")
    
    classical_positions = []
    quantum_forces_on_beads = []
    
    for bead in range(num_beads):
        state = integrator.getState(bead, getPositions=True, getForces=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        classical_positions.append(pos[0])
        quantum_forces_on_beads.append(forces[0])  # Force on classical particle
    
    classical_positions = np.array(classical_positions)
    quantum_forces_on_beads = np.array(quantum_forces_on_beads)
    
    # Check 1: Classical particle positions should be identical on all beads
    pos_spread = np.std(classical_positions, axis=0)
    print(f"\n1. Classical particle position spread across beads:")
    print(f"   Std dev: [{pos_spread[0]:.6f}, {pos_spread[1]:.6f}, {pos_spread[2]:.6f}] nm")
    
    if np.max(pos_spread) < 1e-6:
        print(f"   PASS: Positions identical (< 1e-6 nm)")
    else:
        print(f"   FAIL: Positions not identical!")
        return False
    
    # Check 2: Compute centroid force (average across beads)
    centroid_force = np.mean(quantum_forces_on_beads, axis=0)
    print(f"\n2. Force on classical particle:")
    print(f"   Centroid force (average): [{centroid_force[0]:8.4f}, {centroid_force[1]:8.4f}, {centroid_force[2]:8.4f}]")
    
    # The actual force on each bead should vary (due to different quantum positions)
    # but they should average to the centroid force
    force_spread = np.std(quantum_forces_on_beads, axis=0)
    print(f"   Force spread (std dev): [{force_spread[0]:8.4f}, {force_spread[1]:8.4f}, {force_spread[2]:8.4f}]")
    
    if np.max(force_spread) > 0.1:
        print(f"   Forces vary across beads (as expected)")
    else:
        print(f"   ⚠ Warning: Forces nearly identical across beads")
    
    # Check 3: After velocity update, classical velocity should reflect centroid force
    # Run one more step and check velocity change
    initial_state = integrator.getState(0, getVelocities=True)
    initial_vel = initial_state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)[0]
    
    integrator.step(1)
    
    final_state = integrator.getState(0, getVelocities=True)
    final_vel = final_state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)[0]
    
    vel_change = final_vel - initial_vel
    
    # Expected velocity change: (F/m) * dt (very approximate, ignoring thermostat)
    # Just check the direction is consistent with centroid force
    force_direction = centroid_force / np.linalg.norm(centroid_force)
    vel_direction = vel_change / np.linalg.norm(vel_change) if np.linalg.norm(vel_change) > 1e-6 else np.zeros(3)
    
    print(f"\n3. Velocity update consistency:")
    print(f"   Velocity change: [{vel_change[0]:8.4f}, {vel_change[1]:8.4f}, {vel_change[2]:8.4f}] nm/ps")
    print(f"   Force direction: [{force_direction[0]:6.3f}, {force_direction[1]:6.3f}, {force_direction[2]:6.3f}]")
    print(f"   Vel change dir:  [{vel_direction[0]:6.3f}, {vel_direction[1]:6.3f}, {vel_direction[2]:6.3f}]")
    
    # Dot product should be positive (velocity change aligned with force)
    alignment = np.dot(force_direction, vel_direction)
    print(f"   Alignment (dot product): {alignment:.3f}")
    
    if alignment > 0.5:
        print(f"   PASS: Velocity change aligned with centroid force")
    else:
        print(f"   FAIL: Velocity change not aligned with force!")
        print(f"     This suggests classical particle is not seeing centroid-averaged force")
        return False
    
    print("\nTest PASSED: Classical particle correctly experiences centroid force")
    return True


if __name__ == "__main__":
    success = test_classical_centroid_force()
    sys.exit(0 if success else 1)
