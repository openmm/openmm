#!/usr/bin/env python3
"""
Test: Verify UMA + other forces are correctly accumulated in batched RPMD path

This test checks that when UMA is used with other forces (like HarmonicBondForce),
the total force is correctly computed as the sum of all force contributions.

Bug: Before fix, the batched RPMD path would overwrite UMA forces when copying
non-batched forces, resulting in incorrect total forces.
"""

import sys
import numpy as np
import openmm
from openmm import app, unit

def create_simple_molecule_system():
    """Create a simple 2-atom molecule with a harmonic bond."""
    # Create topology
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("MOL", chain)
    atom1 = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    atom2 = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    topology.addBond(atom1, atom2)
    
    # Initial positions (equilibrium bond length ~1.5 Angstrom)
    positions = [
        [0.0, 0.0, 0.0],
        [0.15, 0.0, 0.0]
    ]
    
    return topology, positions


def test_uma_force_accumulation():
    """Test that UMA forces and bond forces are correctly accumulated."""
    print("=" * 70)
    print("Test: UMA + HarmonicBond Force Accumulation in Batched RPMD")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    topology, positions = create_simple_molecule_system()
    
    # Test parameters
    num_beads = 4
    temperature_K = 300.0
    dt_fs = 1.0
    
    # Create system with UMA
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
    
    # Add a harmonic bond force (in addition to UMA)
    bond_force = openmm.HarmonicBondForce()
    k_bond = 1000.0  # kJ/mol/nm^2
    r0_bond = 0.15   # nm
    bond_force.addBond(0, 1, r0_bond, k_bond)
    bond_force.setForceGroup(1)  # Different force group than UMA (group 0)
    system.addForce(bond_force)
    
    print(f"  System has {system.getNumForces()} forces")
    print(f"  Force 0 (UMA): group {system.getForce(0).getForceGroup()}")
    print(f"  Force 1 (Bond): group {system.getForce(1).getForceGroup()}")
    
    # Create RPMD integrator
    integrator = openmm.RPMDIntegrator(
        num_beads,
        temperature_K * unit.kelvin,
        1.0 / unit.picosecond,
        dt_fs * unit.femtoseconds
    )
    
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
    
    # Initialize bead positions (slightly perturbed)
    for bead in range(num_beads):
        bead_pos = []
        for pos in positions:
            # Add small random perturbation to each bead
            perturb = np.random.randn(3) * 0.001
            bead_pos.append([pos[0] + perturb[0], pos[1] + perturb[1], pos[2] + perturb[2]])
        integrator.setPositions(bead, bead_pos * unit.nanometers)
    
    # Run a single step to compute forces
    print("\nRunning single RPMD step to test force computation...")
    integrator.step(1)
    
    # Get state for each bead and check forces
    print("\nVerifying force accumulation on each bead:")
    all_forces_ok = True
    
    for bead in range(num_beads):
        # Get total force
        state = integrator.getState(bead, getForces=True)
        total_forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        # Get UMA force alone (force group 0)
        state_uma = integrator.getState(bead, getForces=True, groups={0})
        uma_forces = state_uma.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        # Get bond force alone (force group 1)
        state_bond = integrator.getState(bead, getForces=True, groups={1})
        bond_forces = state_bond.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        # Compute expected total (sum of components)
        expected_total = uma_forces + bond_forces
        
        # Check if total matches expected
        force_diff = np.abs(total_forces - expected_total).max()
        force_magnitude = np.abs(total_forces).max()
        relative_error = force_diff / force_magnitude if force_magnitude > 1e-6 else 0.0
        
        print(f"  Bead {bead}:")
        print(f"    UMA force atom 0: [{uma_forces[0][0]:8.4f}, {uma_forces[0][1]:8.4f}, {uma_forces[0][2]:8.4f}]")
        print(f"    Bond force atom 0: [{bond_forces[0][0]:8.4f}, {bond_forces[0][1]:8.4f}, {bond_forces[0][2]:8.4f}]")
        print(f"    Total force atom 0: [{total_forces[0][0]:8.4f}, {total_forces[0][1]:8.4f}, {total_forces[0][2]:8.4f}]")
        print(f"    Expected total: [{expected_total[0][0]:8.4f}, {expected_total[0][1]:8.4f}, {expected_total[0][2]:8.4f}]")
        print(f"    Max absolute error: {force_diff:.6f} kJ/mol/nm")
        print(f"    Relative error: {relative_error*100:.4f}%")
        
        # Check accumulation is correct (tolerance: 1% relative error or 0.1 kJ/mol/nm absolute)
        if relative_error > 0.01 and force_diff > 0.1:
            print(f"    FAIL: Forces not correctly accumulated!")
            all_forces_ok = False
        else:
            print(f"    PASS")
    
    if all_forces_ok:
        print("\nTest PASSED: All forces correctly accumulated")
        return True
    else:
        print("\nTest FAILED: Force accumulation bug detected")
        print("  This indicates the batched RPMD path is overwriting forces")
        print("  instead of adding them.")
        return False


if __name__ == "__main__":
    success = test_uma_force_accumulation()
    sys.exit(0 if success else 1)
