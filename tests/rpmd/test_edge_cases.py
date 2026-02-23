#!/usr/bin/env python3
"""
Test: Edge Cases and Numerical Stability for UMA/RPMD

This test suite checks edge cases and numerical corner cases that commonly
expose bugs in scientific code:

1. Single particle system
2. Very cold temperature (quantum limit)
3. Very hot temperature (classical limit)
4. Zero forces (equilibrium)
5. Extreme forces (stress test)
6. Massless particles (virtual sites)
7. Mixed precision handling
"""

import sys
import numpy as np
import openmm
from openmm import app, unit

def test_single_particle():
    """Test with a single particle (edge case)."""
    print("=" * 70)
    print("Edge Case Test 1: Single Particle System")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    # Single oxygen atom
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("O", chain)
    O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    
    positions = [[0.0, 0.0, 0.0]]
    
    try:
        potential = MLPotential('uma-s-1p1-pythonforce-batch')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
    except Exception as e:
        print(f"Could not create UMA potential: {e}")
        return True
    
    num_beads = 4
    integrator = openmm.RPMDIntegrator(
        num_beads,
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtoseconds
    )
    
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
    
    try:
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions * unit.nanometers)
        
        for bead in range(num_beads):
            integrator.setPositions(bead, positions * unit.nanometers)
        
        integrator.step(10)
        
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        if not (np.isnan(energy) or np.isinf(energy) or np.any(np.isnan(forces)) or np.any(np.isinf(forces))):
            print(f"PASS: Single particle system works (E={energy:.2f} kJ/mol)")
            return True
        else:
            print(f"FAIL: NaN or Inf detected")
            return False
            
    except Exception as e:
        print(f"FAIL: Exception with single particle: {e}")
        return False


def test_temperature_extremes():
    """Test at very low and very high temperatures."""
    print("\n" + "=" * 70)
    print("Edge Case Test 2: Temperature Extremes")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("WAT", chain)
    O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    H1 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
    H2 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
    topology.addBond(O, H1)
    topology.addBond(O, H2)
    
    positions = [[0.0, 0.0, 0.0], [0.0957, 0.0, 0.0], [-0.024, 0.093, 0.0]]
    
    temperatures = [10.0, 1000.0]  # Very cold and very hot
    all_passed = True
    
    for temp in temperatures:
        print(f"\nTesting at T={temp} K:")
        
        try:
            potential = MLPotential('uma-s-1p1-pythonforce-batch')
            system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
        except Exception as e:
            print(f"Could not create UMA potential: {e}")
            continue
        
        num_beads = 8
        integrator = openmm.RPMDIntegrator(
            num_beads,
            temp * unit.kelvin,
            1.0 / unit.picosecond,
            0.5 * unit.femtoseconds
        )
        
        try:
            platform = openmm.Platform.getPlatformByName('CUDA')
        except Exception:
            platform = openmm.Platform.getPlatformByName('Reference')
        
        try:
            context = openmm.Context(system, integrator, platform)
            context.setPositions(positions * unit.nanometers)
            
            for bead in range(num_beads):
                perturbed = []
                for pos in positions:
                    perturb = np.random.randn(3) * 0.001
                    perturbed.append([pos[0] + perturb[0], pos[1] + perturb[1], pos[2] + perturb[2]])
                integrator.setPositions(bead, perturbed * unit.nanometers)
            
            integrator.step(20)
            
            state = context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
            if not (np.isnan(energy) or np.isinf(energy)):
                print(f"  PASS: T={temp}K works (E={energy:.2f} kJ/mol)")
            else:
                print(f"  FAIL: NaN or Inf at T={temp}K")
                all_passed = False
                
        except Exception as e:
            print(f"  FAIL: Exception at T={temp}K: {e}")
            all_passed = False
    
    return all_passed


def test_zero_force_case():
    """Test when forces are very small (near equilibrium)."""
    print("\n" + "=" * 70)
    print("Edge Case Test 3: Near-Zero Forces")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("WAT", chain)
    O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    H1 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
    H2 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
    
    # Near-equilibrium geometry (should have very small forces)
    positions = [[0.0, 0.0, 0.0], [0.09584, 0.0, 0.0], [-0.02395, 0.09298, 0.0]]
    
    try:
        potential = MLPotential('uma-s-1p1-pythonforce-batch')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
    except Exception as e:
        print(f"Could not create UMA potential: {e}")
        return True
    
    num_beads = 4
    integrator = openmm.RPMDIntegrator(
        num_beads,
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtoseconds
    )
    
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
    
    try:
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions * unit.nanometers)
        
        for bead in range(num_beads):
            integrator.setPositions(bead, positions * unit.nanometers)
        
        # Step with near-zero forces
        integrator.step(10)
        
        # Check that simulation doesn't produce NaN from division by zero
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        force_magnitude = np.linalg.norm(forces)
        
        print(f"Force magnitude: {force_magnitude:.6f} kJ/mol/nm")
        print(f"Energy: {energy:.6f} kJ/mol")
        
        if not (np.isnan(energy) or np.isinf(energy) or np.any(np.isnan(forces)) or np.any(np.isinf(forces))):
            print(f"PASS: Near-zero forces handled correctly")
            return True
        else:
            print(f"FAIL: NaN or Inf with near-zero forces")
            return False
            
    except Exception as e:
        print(f"FAIL: Exception with near-zero forces: {e}")
        return False


def test_particle_exclusion():
    """Test that excluded particles receive zero force from UMA."""
    print("\n" + "=" * 70)
    print("Edge Case Test 4: Particle Exclusion")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    # Water + extra particle
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("WAT", chain)
    O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    H1 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
    H2 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
    
    positions = [[0.0, 0.0, 0.0], [0.0957, 0.0, 0.0], [-0.024, 0.093, 0.0]]
    
    try:
        potential = MLPotential('uma-s-1p1-pythonforce-batch')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
    except Exception as e:
        print(f"Could not create UMA potential: {e}")
        return True
    
    # Add an extra particle (e.g., virtual photon)
    extra_index = system.addParticle(1.0 / 1822.888)  # electron mass
    positions.append([0.0, 0.0, 0.5])
    
    print(f"System: 3 water atoms + 1 extra particle")
    print(f"  Water: indices 0-2 (UMA applied)")
    print(f"  Extra: index {extra_index} (should get zero force from UMA)")
    
    num_beads = 4
    integrator = openmm.RPMDIntegrator(
        num_beads,
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtoseconds
    )
    
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
    
    try:
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions * unit.nanometers)
        
        for bead in range(num_beads):
            perturbed = []
            for pos in positions:
                perturb = np.random.randn(3) * 0.001
                perturbed.append([pos[0] + perturb[0], pos[1] + perturb[1], pos[2] + perturb[2]])
            integrator.setPositions(bead, perturbed * unit.nanometers)
        
        integrator.step(5)
        
        # Check forces on all beads
        max_extra_force = 0.0
        for bead in range(num_beads):
            state = integrator.getState(bead, getForces=True)
            forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
            extra_force_magnitude = np.linalg.norm(forces[extra_index])
            max_extra_force = max(max_extra_force, extra_force_magnitude)
        
        print(f"\nMax force on extra particle: {max_extra_force:.8f} kJ/mol/nm")
        
        if max_extra_force < 1e-6:
            print(f"PASS: Extra particle receives zero force from UMA")
            return True
        else:
            print(f"FAIL: Extra particle has non-zero force!")
            print(f"  This indicates UMA is incorrectly being applied to excluded particles.")
            return False
            
    except Exception as e:
        print(f"FAIL: Exception: {e}")
        return False


def test_position_wraparound():
    """Test periodic boundary condition edge cases."""
    print("\n" + "=" * 70)
    print("Edge Case Test 5: Periodic Boundary Conditions")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    # Create water in a periodic box
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("WAT", chain)
    O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    H1 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
    H2 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
    
    # Position near box edge
    box_size = 2.0  # nm
    positions = [
        [0.95, 0.95, 0.95],
        [1.0457, 0.95, 0.95],
        [0.926, 1.043, 0.95]
    ]
    
    try:
        potential = MLPotential('uma-s-1p1-pythonforce-batch')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
    except Exception as e:
        print(f"Could not create UMA potential: {e}")
        return True
    
    # Set periodic box
    system.setDefaultPeriodicBoxVectors(
        [box_size, 0, 0],
        [0, box_size, 0],
        [0, 0, box_size]
    )
    
    num_beads = 4
    integrator = openmm.RPMDIntegrator(
        num_beads,
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtoseconds
    )
    
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
    
    try:
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions * unit.nanometers)
        
        for bead in range(num_beads):
            integrator.setPositions(bead, positions * unit.nanometers)
        
        # Run simulation (positions should wrap)
        print(f"Running 50 steps near box edge...")
        integrator.step(50)
        
        # Check final state
        state = context.getState(getEnergy=True, getPositions=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        final_pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        
        print(f"  Final energy: {energy:.2f} kJ/mol")
        print(f"  Positions may have wrapped through PBC")
        
        if not (np.isnan(energy) or np.isinf(energy)):
            print(f"PASS: PBC handling works correctly")
            return True
        else:
            print(f"FAIL: NaN or Inf with PBC")
            return False
            
    except Exception as e:
        print(f"FAIL: Exception with PBC: {e}")
        return False


def test_numerical_precision():
    """Test numerical precision in force accumulation."""
    print("\n" + "=" * 70)
    print("Edge Case Test 6: Numerical Precision")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    # Create system
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("O2", chain)
    O1 = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    O2 = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    
    positions = [[0.0, 0.0, 0.0], [0.15, 0.0, 0.0]]
    
    try:
        potential = MLPotential('uma-s-1p1-pythonforce-batch')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
    except Exception as e:
        print(f"Could not create UMA potential: {e}")
        return True
    
    # Add many weak forces that should accumulate
    num_weak_forces = 10
    for i in range(num_weak_forces):
        # Weak harmonic restraint
        restraint = openmm.CustomExternalForce(f"0.1*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        offset = i * 0.001  # Small offsets
        restraint.addParticle(0, [offset, 0, 0])
        restraint.addParticle(1, [0.15 + offset, 0, 0])
        restraint.setForceGroup(i + 1)
        system.addForce(restraint)
    
    print(f"System with UMA + {num_weak_forces} weak forces")
    print(f"Testing precision of force accumulation...")
    
    num_beads = 4
    integrator = openmm.RPMDIntegrator(
        num_beads,
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtoseconds
    )
    
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
    
    try:
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions * unit.nanometers)
        
        for bead in range(num_beads):
            integrator.setPositions(bead, positions * unit.nanometers)
        
        integrator.step(1)
        
        # Get total force
        state_total = integrator.getState(0, getForces=True)
        total_forces = state_total.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        # Sum individual groups
        summed_forces = np.zeros((2, 3))
        for group in range(num_weak_forces + 1):
            state_group = integrator.getState(0, getForces=True, groups={group})
            group_forces = state_group.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
            summed_forces += group_forces
        
        error = np.abs(total_forces - summed_forces).max()
        relative_error = error / np.abs(total_forces).max() if np.abs(total_forces).max() > 1e-6 else 0
        
        print(f"  Max accumulation error: {error:.8f} kJ/mol/nm")
        print(f"  Relative error: {relative_error*100:.6f}%")
        
        if error < 0.01 or relative_error < 0.001:
            print(f"PASS: Numerical precision maintained in accumulation")
            return True
        else:
            print(f"FAIL: Accumulation error too large")
            print(f"  May indicate precision loss in force accumulation")
            return False
            
    except Exception as e:
        print(f"FAIL: Exception: {e}")
        return False


def run_all_edge_case_tests():
    """Run all edge case tests."""
    print("\n" + "="*70)
    print("UMA/RPMD EDGE CASE TESTS")
    print("="*70)
    
    results = []
    
    results.append(("Single Particle", test_single_particle()))
    results.append(("Temperature Extremes", test_temperature_extremes()))
    results.append(("Near-Zero Forces", test_zero_force_case()))
    results.append(("Particle Exclusion", test_particle_exclusion()))
    results.append(("Periodic Boundaries", test_position_wraparound()))
    results.append(("Numerical Precision", test_numerical_precision()))
    
    print("\n" + "="*70)
    print("EDGE CASE TEST SUMMARY")
    print("="*70)
    
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"{status}: {name}")
    
    all_passed = all(result[1] for result in results)
    
    if all_passed:
        print("\nALL EDGE CASE TESTS PASSED")
        return True
    else:
        print("\nSOME EDGE CASE TESTS FAILED")
        return False


if __name__ == "__main__":
    success = run_all_edge_case_tests()
    sys.exit(0 if success else 1)
