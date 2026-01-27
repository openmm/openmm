#!/usr/bin/env python3
"""
Test script for RPMD with UMA potentials.
Tests UMA model compatibility with RPMD integrator including GPU execution.
"""
import sys
import os
import openmm
from openmm import unit
import openmm.app as app
import numpy as np

def test_rpmd_uma_basic():
    """Test basic RPMD with UMA potential"""
    print("=" * 60)
    print("Testing RPMD with UMA Potential")
    print("=" * 60)
    
    # Import UMA after imports are done
    try:
        from openmmml import MLPotential
    except ImportError as e:
        print(f"Failed to import openmmml: {e}")
        return 1
    
    # Load a small test system
    test_dir = os.path.join(os.path.dirname(__file__), "fairchem", "openmm-ml", "test", "data")
    if not os.path.exists(test_dir):
        # Try alternative path
        test_dir = os.path.join(os.path.dirname(__file__), "tests", "dimer_system")
        if not os.path.exists(test_dir):
            print(f"Test data directory not found. Tried: {test_dir}")
            print("Creating a simple test system instead...")
            # Create a simple water molecule
            from openmm.app import Topology, Element
            topology = Topology()
            chain = topology.addChain()
            residue = topology.addResidue("WAT", chain)
            O = topology.addAtom("O", Element.getBySymbol("O"), residue)
            H1 = topology.addAtom("H1", Element.getBySymbol("H"), residue)
            H2 = topology.addAtom("H2", Element.getBySymbol("H"), residue)
            topology.addBond(O, H1)
            topology.addBond(O, H2)
            
            positions = np.array([
                [0.0, 0.0, 0.0],
                [0.0957, 0.0, 0.0],
                [-0.024, 0.093, 0.0]
            ]) * unit.nanometers
    else:
        pdb_file = os.path.join(test_dir, "toluene", "toluene.pdb")
        if os.path.exists(pdb_file):
            pdb = app.PDBFile(pdb_file)
            topology = pdb.topology
            positions = pdb.getPositions(asNumpy=True)
        else:
            print(f"Test PDB file not found: {pdb_file}")
            return 1
    
    num_particles = topology.getNumAtoms()
    num_beads = 4  # Small number for testing
    print(f"\nSystem: {num_particles} particles, {num_beads} beads")
    
    # Test 1: Create UMA system
    print(f"\n[TEST 1] Creating UMA system...")
    try:
        potential = MLPotential('uma-s-1p1')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
        print(f"  ✓ System created with {system.getNumForces()} forces")
    except Exception as e:
        print(f"  ✗ Failed to create UMA system: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Test 2: Create RPMD integrator
    print(f"\n[TEST 2] Creating RPMD integrator...")
    temperature = 300.0
    friction = 1.0
    dt = 0.5
    
    print(f"  Temperature: {temperature} K")
    print(f"  Friction: {friction} ps^-1")
    print(f"  Time step: {dt} fs")
    print(f"  Number of beads: {num_beads}")
    
    try:
        integrator = openmm.RPMDIntegrator(
            num_beads,
            temperature * unit.kelvin,
            friction / unit.picosecond,
            dt * unit.femtosecond
        )
        print("  ✓ RPMD integrator created")
    except Exception as e:
        print(f"  ✗ Failed to create RPMD integrator: {e}")
        return 1
    
    # Test 3: Create context (try CUDA first, fall back to CPU)
    print(f"\n[TEST 3] Creating context...")
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        print("  Using CUDA platform")
    except:
        try:
            platform = openmm.Platform.getPlatformByName('OpenCL')
            print("  Using OpenCL platform")
        except:
            platform = openmm.Platform.getPlatformByName('CPU')
            print("  Using CPU platform")
    
    try:
        context = openmm.Context(system, integrator, platform)
        print("  ✓ Context created successfully")
    except Exception as e:
        print(f"  ✗ Failed to create context: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Test 4: Set initial positions for all beads
    print(f"\n[TEST 4] Setting initial positions...")
    try:
        # Set positions for the context
        context.setPositions(positions)
        
        # Set positions for each bead with small perturbations
        for bead in range(num_beads):
            bead_positions = []
            for p in positions:
                px = p[0].value_in_unit(unit.nanometers) + 0.001 * np.random.randn()
                py = p[1].value_in_unit(unit.nanometers) + 0.001 * np.random.randn()
                pz = p[2].value_in_unit(unit.nanometers) + 0.001 * np.random.randn()
                bead_positions.append(openmm.Vec3(px, py, pz) * unit.nanometers)
            integrator.setPositions(bead, bead_positions)
        print("  ✓ Positions set for all beads")
    except Exception as e:
        print(f"  ✗ Failed to set positions: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Test 5: Initialize velocities
    print(f"\n[TEST 5] Initializing velocities...")
    try:
        context.setVelocitiesToTemperature(temperature * unit.kelvin)
        
        # Get velocities from context
        state = context.getState(getVelocities=True)
        base_velocities = state.getVelocities()
        
        # Set slightly perturbed velocities for each bead
        for bead in range(num_beads):
            bead_velocities = []
            for v in base_velocities:
                vx = v[0].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01 * np.random.randn()
                vy = v[1].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01 * np.random.randn()
                vz = v[2].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01 * np.random.randn()
                bead_velocities.append(openmm.Vec3(vx, vy, vz) * unit.nanometers/unit.picoseconds)
            integrator.setVelocities(bead, bead_velocities)
        print("  ✓ Velocities initialized for all beads")
    except Exception as e:
        print(f"  ✗ Failed to initialize velocities: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Test 6: Run simulation
    print(f"\n[TEST 6] Running RPMD simulation (50 steps)...")
    try:
        initial_energy = integrator.getTotalEnergy()
        print(f"  Initial total energy: {initial_energy} kJ/mol")
        
        # Get initial state for one bead
        state0_initial = integrator.getState(0, getEnergy=True)
        pe_initial = state0_initial.getPotentialEnergy()
        print(f"  Initial potential energy (bead 0): {pe_initial}")
        
        # Run simulation
        integrator.step(50)
        
        final_energy = integrator.getTotalEnergy()
        print(f"  Final total energy: {final_energy} kJ/mol")
        
        # Get final state
        state0_final = integrator.getState(0, getEnergy=True)
        pe_final = state0_final.getPotentialEnergy()
        print(f"  Final potential energy (bead 0): {pe_final}")
        
        print("  ✓ Simulation completed successfully")
    except Exception as e:
        print(f"  ✗ Simulation failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Test 7: Check energies are finite
    print(f"\n[TEST 7] Checking energies are finite...")
    try:
        assert np.isfinite(final_energy), f"Total energy is not finite: {final_energy}"
        assert np.isfinite(pe_final.value_in_unit(unit.kilojoules_per_mole)), \
            f"Potential energy is not finite: {pe_final}"
        print("  ✓ All energies are finite")
    except AssertionError as e:
        print(f"  ✗ Energy check failed: {e}")
        return 1
    
    print("\n" + "=" * 60)
    print("ALL TESTS PASSED!")
    print("=" * 60)
    return 0


def test_rpmd_uma_pile_g():
    """Test RPMD with UMA potential and PILE_G thermostat"""
    print("\n" + "=" * 60)
    print("Testing RPMD with UMA and PILE_G Thermostat")
    print("=" * 60)
    
    try:
        from openmmml import MLPotential
        from openmm.app import Topology, Element
    except ImportError as e:
        print(f"Failed to import required modules: {e}")
        return 1
    
    # Create a simple test system
    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("WAT", chain)
    O = topology.addAtom("O", Element.getBySymbol("O"), residue)
    H1 = topology.addAtom("H1", Element.getBySymbol("H"), residue)
    H2 = topology.addAtom("H2", Element.getBySymbol("H"), residue)
    topology.addBond(O, H1)
    topology.addBond(O, H2)
    
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.0957, 0.0, 0.0],
        [-0.024, 0.093, 0.0]
    ]) * unit.nanometers
    
    print("\nSystem: water molecule, 4 beads")
    
    # Create UMA system
    print("\n[TEST 1] Creating UMA system with PILE_G thermostat...")
    try:
        potential = MLPotential('uma-s-1p1')
        system = potential.createSystem(topology, task_name='omol')
        
        integrator = openmm.RPMDIntegrator(
            4,  # num_beads
            300.0 * unit.kelvin,
            1.0 / unit.picosecond,
            0.5 * unit.femtosecond
        )
        
        # Set to PILE_G thermostat
        integrator.setThermostatType(openmm.RPMDIntegrator.PileG)
        assert integrator.getThermostatType() == openmm.RPMDIntegrator.PileG
        print("  ✓ PILE_G thermostat set")
        
        # Set centroid friction
        integrator.setCentroidFriction(0.5 / unit.picoseconds)
        print("  ✓ Centroid friction configured")
        
        # Create context
        platform = openmm.Platform.getPlatform(0)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions)
        
        # Set positions for all beads
        for bead in range(4):
            bead_pos = [openmm.Vec3(p[0].value_in_unit(unit.nanometers) + 0.001*np.random.randn(),
                                    p[1].value_in_unit(unit.nanometers) + 0.001*np.random.randn(),
                                    p[2].value_in_unit(unit.nanometers) + 0.001*np.random.randn())
                       * unit.nanometers for p in positions]
            integrator.setPositions(bead, bead_pos)
        
        context.setVelocitiesToTemperature(300 * unit.kelvin)
        
        # Run short simulation
        integrator.step(20)
        
        energy = integrator.getTotalEnergy()
        assert np.isfinite(energy), f"Energy not finite: {energy}"
        
        print("  ✓ PILE_G simulation completed successfully")
        print(f"  Final energy: {energy} kJ/mol")
        
    except Exception as e:
        print(f"  ✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    print("\n" + "=" * 60)
    print("PILE_G TEST PASSED!")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    # Run basic test
    result1 = test_rpmd_uma_basic()
    
    # Run PILE_G test
    result2 = test_rpmd_uma_pile_g()
    
    sys.exit(max(result1, result2))
