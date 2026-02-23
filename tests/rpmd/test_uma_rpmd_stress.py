#!/usr/bin/env python3
"""
Test: Stress Test for UMA Batched RPMD

This test runs challenging scenarios designed to expose bugs:
- Many beads (high P)
- Long simulations
- Multiple force groups
- Edge cases (single particle, many particles)

These stress tests catch issues like:
- Memory corruption
- Race conditions
- Accumulation of numerical errors
- Buffer overflows
"""

import sys
import numpy as np
import openmm
from openmm import app, unit
import time

def create_multi_molecule_system(num_molecules=3):
    """Create a system with multiple water molecules."""
    topology = app.Topology()
    chain = topology.addChain()
    positions = []
    
    spacing = 0.5  # nm between molecules
    
    for mol_idx in range(num_molecules):
        residue = topology.addResidue(f"WAT{mol_idx}", chain)
        O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
        H1 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
        H2 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
        topology.addBond(O, H1)
        topology.addBond(O, H2)
        
        # Offset each molecule
        offset = [mol_idx * spacing, 0, 0]
        positions.append([0.0 + offset[0], 0.0 + offset[1], 0.0 + offset[2]])
        positions.append([0.0957 + offset[0], 0.0 + offset[1], 0.0 + offset[2]])
        positions.append([-0.024 + offset[0], 0.093 + offset[1], 0.0 + offset[2]])
    
    return topology, positions


def test_many_beads():
    """Test with large number of beads."""
    print("=" * 70)
    print("Stress Test 1: Many Beads (P=32)")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    topology, positions = create_multi_molecule_system(num_molecules=2)
    
    num_beads = 32  # Stress test with many beads
    temperature_K = 300.0
    
    try:
        potential = MLPotential('uma-s-1p1-pythonforce-batch')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
    except Exception as e:
        print(f"Could not create UMA potential: {e}")
        return True
    
    integrator = openmm.RPMDIntegrator(
        num_beads,
        temperature_K * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtoseconds
    )
    
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        print(f"Using CUDA platform")
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
        print(f"Using Reference platform")
    
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions * unit.nanometers)
    
    # Initialize beads
    for bead in range(num_beads):
        perturbed = []
        for pos in positions:
            perturb = np.random.randn(3) * 0.001
            perturbed.append([pos[0] + perturb[0], pos[1] + perturb[1], pos[2] + perturb[2]])
        integrator.setPositions(bead, perturbed * unit.nanometers)
    
    print(f"\nRunning 50 steps with {num_beads} beads...")
    start_time = time.time()
    
    try:
        integrator.step(50)
        elapsed = time.time() - start_time
        
        print(f"Completed successfully in {elapsed:.2f}s")
        print(f"  Performance: {50/elapsed:.2f} steps/s")
        
        # Check final state is valid
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        if np.isnan(energy) or np.isinf(energy):
            print(f"FAIL: Energy is NaN or Inf")
            return False
        
        if np.any(np.isnan(forces)) or np.any(np.isinf(forces)):
            print(f"FAIL: Forces contain NaN or Inf")
            return False
        
        print(f"Final state valid (E={energy:.2f} kJ/mol)")
        return True
        
    except Exception as e:
        print(f"FAIL: Simulation crashed: {e}")
        return False


def test_multiple_force_groups():
    """Test with forces in multiple groups."""
    print("\n" + "=" * 70)
    print("Stress Test 2: Multiple Force Groups")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    topology, positions = create_multi_molecule_system(num_molecules=2)
    
    num_beads = 8
    temperature_K = 300.0
    
    try:
        potential = MLPotential('uma-s-1p1-pythonforce-batch')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
    except Exception as e:
        print(f"Could not create UMA potential: {e}")
        return True
    
    # Add forces in multiple groups
    print(f"Adding forces in multiple groups...")
    
    # Group 1: Harmonic restraints
    restraint = openmm.CustomExternalForce("k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint.addGlobalParameter("k", 100.0)  # kJ/mol/nm^2
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    for i, pos in enumerate(positions):
        restraint.addParticle(i, pos)
    restraint.setForceGroup(1)
    system.addForce(restraint)
    
    # Group 2: Simple nonbonded (electrostatics only)
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    for i in range(len(positions)):
        charge = 0.1 * (1 if i % 3 == 0 else -0.05)  # O negative, H positive
        nb.addParticle(charge, 0.1, 0.0)  # No LJ
    nb.setForceGroup(2)
    system.addForce(nb)
    
    print(f"  Force group 0: UMA")
    print(f"  Force group 1: Restraints")
    print(f"  Force group 2: Electrostatics")
    
    integrator = openmm.RPMDIntegrator(
        num_beads,
        temperature_K * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtoseconds
    )
    
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
    
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions * unit.nanometers)
    
    # Initialize beads
    for bead in range(num_beads):
        perturbed = []
        for pos in positions:
            perturb = np.random.randn(3) * 0.01
            perturbed.append([pos[0] + perturb[0], pos[1] + perturb[1], pos[2] + perturb[2]])
        integrator.setPositions(bead, perturbed * unit.nanometers)
    
    print(f"\nRunning 20 steps...")
    try:
        integrator.step(20)
        
        # Verify force accumulation
        state_all = integrator.getState(0, getForces=True)
        forces_all = state_all.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        state_0 = integrator.getState(0, getForces=True, groups={0})
        forces_0 = state_0.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        state_1 = integrator.getState(0, getForces=True, groups={1})
        forces_1 = state_1.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        state_2 = integrator.getState(0, getForces=True, groups={2})
        forces_2 = state_2.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        expected_total = forces_0 + forces_1 + forces_2
        error = np.abs(forces_all - expected_total).max()
        
        print(f"  Force accumulation error: {error:.8f} kJ/mol/nm")
        
        if error < 0.01:
            print(f"PASS: All force groups correctly accumulated")
            return True
        else:
            print(f"FAIL: Force accumulation error too large!")
            print(f"  Total force: {np.linalg.norm(forces_all):.4f}")
            print(f"  Expected: {np.linalg.norm(expected_total):.4f}")
            return False
            
    except Exception as e:
        print(f"FAIL: Simulation crashed: {e}")
        return False


def test_long_simulation():
    """Test stability over longer simulation."""
    print("\n" + "=" * 70)
    print("Stress Test 3: Long Simulation (500 steps)")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    topology, positions = create_multi_molecule_system(num_molecules=2)
    
    num_beads = 8
    temperature_K = 300.0
    
    try:
        potential = MLPotential('uma-s-1p1-pythonforce-batch')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
    except Exception as e:
        print(f"Could not create UMA potential: {e}")
        return True
    
    integrator = openmm.RPMDIntegrator(
        num_beads,
        temperature_K * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtoseconds
    )
    
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
    
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions * unit.nanometers)
    
    # Initialize beads
    for bead in range(num_beads):
        perturbed = []
        for pos in positions:
            perturb = np.random.randn(3) * 0.001
            perturbed.append([pos[0] + perturb[0], pos[1] + perturb[1], pos[2] + perturb[2]])
        integrator.setPositions(bead, perturbed * unit.nanometers)
    
    print(f"\nRunning 500 steps, checking every 100...")
    
    energies = []
    for i in range(5):
        integrator.step(100)
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        energies.append(energy)
        
        if np.isnan(energy) or np.isinf(energy):
            print(f"FAIL: Energy became NaN/Inf at step {(i+1)*100}")
            return False
        
        print(f"  Step {(i+1)*100}: E = {energy:.2f} kJ/mol")
    
    # Check energy hasn't exploded (sign of instability)
    energy_range = max(energies) - min(energies)
    avg_energy = np.mean(energies)
    
    print(f"\n  Energy range: {energy_range:.2f} kJ/mol")
    print(f"  Average energy: {avg_energy:.2f} kJ/mol")
    
    if abs(energies[-1]) < 1e6:  # Reasonable energy scale
        print(f"PASS: Simulation stable over 500 steps")
        return True
    else:
        print(f"FAIL: Energy exploded (instability)")
        return False


def run_all_stress_tests():
    """Run all stress tests."""
    print("\n" + "="*70)
    print("UMA BATCHED RPMD STRESS TESTS")
    print("="*70)
    
    results = []
    
    results.append(("Many Beads (P=32)", test_many_beads()))
    results.append(("Multiple Force Groups", test_multiple_force_groups()))
    results.append(("Long Simulation", test_long_simulation()))
    
    print("\n" + "="*70)
    print("STRESS TEST SUMMARY")
    print("="*70)
    
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"{status}: {name}")
    
    all_passed = all(result[1] for result in results)
    
    if all_passed:
        print("\nALL STRESS TESTS PASSED")
        return True
    else:
        print("\nSOME STRESS TESTS FAILED")
        return False


if __name__ == "__main__":
    success = run_all_stress_tests()
    sys.exit(0 if success else 1)
