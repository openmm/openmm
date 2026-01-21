#!/usr/bin/env python3
"""
Diagnostic: Check if cavity displacement matches expected value from dipole
"""
import sys
import numpy as np
from openmm import openmm, unit

BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5

def create_simple_system():
    """Create a simple system with a few charged particles."""
    system = openmm.System()
    positions = []
    
    # Add 10 O-O dimers
    mass_O = 16.0
    charge_O = -0.3
    bond_length = 0.12  # nm
    
    # Arrange in a line for simplicity
    for i in range(10):
        x = i * 0.3  # Space them out
        
        system.addParticle(mass_O)
        system.addParticle(mass_O)
        positions.append(openmm.Vec3(x, 0.0, 0.0) * unit.nanometer)
        positions.append(openmm.Vec3(x + bond_length, 0.0, 0.0) * unit.nanometer)
    
    # Add nonbonded force
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(1.0 * unit.nanometer)
    
    for i in range(20):
        nb.addParticle(charge_O, 0.1, 0.0)
        if i % 2 == 1:  # Exception for bonded pairs
            nb.addException(i-1, i, 0.0, 1.0, 0.0)
    
    system.addForce(nb)
    
    # Periodic box
    box_size = 5.0
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_size, 0, 0),
        openmm.Vec3(0, box_size, 0),
        openmm.Vec3(0, 0, box_size)
    )
    
    return system, positions

def calculate_dipole(context, cavity_idx):
    """Calculate molecular dipole moment (excluding cavity)."""
    state = context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    
    # Get charges from system
    system = context.getSystem()
    charges = []
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for i in range(system.getNumParticles()):
                if i == cavity_idx:
                    charges.append(0.0)
                else:
                    q, _, _ = force.getParticleParameters(i)
                    charges.append(q._value)  # Get underlying value
            break
    
    # Calculate dipole (excluding cavity)
    dipole = np.zeros(3)
    for i, (p, q) in enumerate(zip(pos, charges)):
        if i != cavity_idx:
            dipole[0] += q * p[0]
            dipole[1] += q * p[1]
            dipole[2] += q * p[2]
    
    return dipole

def main():
    print("=" * 70)
    print("Cavity Displacement Diagnostic")
    print("Checking if displacement magnitude matches dipole expectation")
    print("=" * 70)
    
    # Create system
    system, positions = create_simple_system()
    num_particles = system.getNumParticles()
    print(f"\nSystem: {num_particles} particles")
    
    # Add cavity particle
    cavity_idx = system.addParticle(1.0)  # 1 amu
    positions.append(openmm.Vec3(0.0, 0.0, 0.0) * unit.nanometer)
    
    # Add cavity to nonbonded (charge=0)
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.addParticle(0.0, 0.1, 0.0)
    
    print(f"Cavity particle index: {cavity_idx}")
    
    # Cavity parameters
    omegac_au = 0.00913  # Hartree
    lambda_coupling = 0.001
    photon_mass = 1.0
    
    print(f"\nCavity parameters:")
    print(f"  omega_c: {omegac_au} a.u.")
    print(f"  lambda: {lambda_coupling}")
    print(f"  photon mass: {photon_mass} amu")
    
    # Add cavity force
    cavity_force = openmm.CavityForce(cavity_idx, omegac_au, lambda_coupling, photon_mass)
    system.addForce(cavity_force)
    
    # Create simulation
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName('CUDA')
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(100.0 * unit.kelvin)
    
    print("\n--- Running Short Simulation ---")
    print("Steps | Cavity Position (nm) | Dipole (e*nm) | Expected q_eq (nm)")
    print("-" * 70)
    
    # Run for a few steps and track
    for step in range(0, 501, 100):
        if step > 0:
            integrator.step(100)
        
        # Get current state
        state = context.getState(getPositions=True, getEnergy=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Cavity position
        cavity_pos = pos[cavity_idx]
        
        # Calculate molecular dipole
        dipole = calculate_dipole(context, cavity_idx)
        dipole_mag = np.linalg.norm(dipole[:2])  # Only x,y matter
        
        # Expected equilibrium position: q_eq = -(lambda/omega) * d
        # In SI-like units for OpenMM
        factor = -lambda_coupling / omegac_au
        expected_q = factor * dipole
        expected_q_mag = np.linalg.norm(expected_q[:2])
        
        # Actual cavity displacement magnitude
        actual_q_mag = np.linalg.norm(cavity_pos[:2])
        
        print(f"{step:5d} | ({cavity_pos[0]:8.4f}, {cavity_pos[1]:8.4f}, {cavity_pos[2]:8.4f}) | " 
              f"({dipole[0]:7.4f}, {dipole[1]:7.4f}, {dipole[2]:7.4f}) | "
              f"({expected_q[0]:8.4f}, {expected_q[1]:8.4f})")
        
        if step > 0:
            ratio = actual_q_mag / expected_q_mag if expected_q_mag > 1e-10 else 0
            print(f"      | Actual |q|: {actual_q_mag:8.4f} nm | Expected |q_eq|: {expected_q_mag:8.4f} nm | Ratio: {ratio:.2f}x")
    
    # Final analysis
    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)
    state = context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    cavity_pos = pos[cavity_idx]
    dipole = calculate_dipole(context, cavity_idx)
    
    factor = -lambda_coupling / omegac_au
    expected_q = factor * dipole
    
    actual_mag = np.linalg.norm(cavity_pos[:2])
    expected_mag = np.linalg.norm(expected_q[:2])
    
    print(f"\nFinal cavity position: ({cavity_pos[0]:.4f}, {cavity_pos[1]:.4f}, {cavity_pos[2]:.4f}) nm")
    print(f"Final dipole moment: ({dipole[0]:.4f}, {dipole[1]:.4f}, {dipole[2]:.4f}) e*nm")
    print(f"\nExpected q_eq = -(lambda/omega) * d:")
    print(f"  q_eq = -({lambda_coupling}/{omegac_au}) * d")
    print(f"  q_eq = {factor:.4f} * d")
    print(f"  q_eq = ({expected_q[0]:.4f}, {expected_q[1]:.4f}) nm")
    
    print(f"\nMagnitude comparison:")
    print(f"  |q_actual| = {actual_mag:.4f} nm")
    print(f"  |q_expected| = {expected_mag:.4f} nm")
    print(f"  Ratio: {actual_mag/expected_mag if expected_mag > 1e-10 else float('inf'):.2f}x")
    
    if actual_mag / expected_mag > 10:
        print("\n⚠️  WARNING: Actual displacement is >10x larger than expected!")
        print("   This suggests a bug in force calculation or unit conversion.")
    elif actual_mag / expected_mag > 2:
        print("\n⚠️  Actual displacement is larger than expected")
        print("   May indicate force imbalance or incorrect dynamics")
    else:
        print("\n✓ Displacement magnitude is reasonable")
    
    print("=" * 70)

if __name__ == '__main__':
    try:
        import openmm
        print("✓ OpenMM loaded successfully\n")
    except ImportError as e:
        print(f"✗ Failed to import OpenMM: {e}")
        sys.exit(1)
    
    main()
