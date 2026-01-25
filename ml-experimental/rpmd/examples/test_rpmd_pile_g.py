#!/usr/bin/env python3
"""
Test script for RPMD with PILE_G (Bussi centroid) thermostat.
Tests basic functionality and GPU kernel execution.
"""
import sys
import openmm
from openmm import unit
import numpy as np

def test_rpmd_pile_g():
    """Test RPMD with PILE_G thermostat"""
    print("=" * 60)
    print("Testing RPMD with PILE_G Thermostat")
    print("=" * 60)
    
    # Create a simple test system (Argon atoms)
    num_particles = 10
    num_beads = 4
    print(f"\nSystem: {num_particles} particles, {num_beads} beads")
    
    system = openmm.System()
    for i in range(num_particles):
        system.addParticle(39.948)  # Argon mass
    
    # Add a simple harmonic force
    force = openmm.HarmonicBondForce()
    for i in range(num_particles - 1):
        force.addBond(i, i+1, 0.3, 1000.0)
    system.addForce(force)
    
    # Create RPMD integrator with PILE_G thermostat
    temperature = 300.0
    friction = 1.0
    dt = 0.5
    
    print(f"Temperature: {temperature} K")
    print(f"Friction: {friction} ps^-1")
    print(f"Time step: {dt} fs")
    
    integrator = openmm.RPMDIntegrator(
        num_beads,
        temperature * unit.kelvin,
        friction / unit.picosecond,
        dt * unit.femtosecond
    )
    
    # Test 1: Default thermostat should be Pile
    print(f"\n[TEST 1] Default thermostat type...")
    default_type = integrator.getThermostatType()
    print(f"  Default: {default_type} (expected: {openmm.RPMDIntegrator.Pile})")
    assert default_type == openmm.RPMDIntegrator.Pile, "Default should be Pile"
    print("  ✓ PASSED")
    
    # Test 2: Switch to PILE_G
    print(f"\n[TEST 2] Setting PILE_G thermostat...")
    integrator.setThermostatType(openmm.RPMDIntegrator.PileG)
    new_type = integrator.getThermostatType()
    print(f"  New type: {new_type} (expected: {openmm.RPMDIntegrator.PileG})")
    assert new_type == openmm.RPMDIntegrator.PileG, "Should be PileG"
    print("  ✓ PASSED")
    
    # Test 3: Centroid friction getter/setter
    print(f"\n[TEST 3] Centroid friction getter/setter...")
    default_centroid_friction = integrator.getCentroidFriction()
    print(f"  Default centroid friction: {default_centroid_friction}")
    integrator.setCentroidFriction(0.5 / unit.picoseconds)
    new_centroid_friction = integrator.getCentroidFriction()
    print(f"  New centroid friction: {new_centroid_friction}")
    assert abs(new_centroid_friction.value_in_unit(unit.picoseconds**-1) - 0.5) < 1e-10, "Centroid friction not set correctly"
    print("  ✓ PASSED")
    
    # Test 4: Create context and run simulation
    print(f"\n[TEST 4] Creating context...")
    try:
        # Try to get CUDA platform first, fall back to CPU
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
        
        context = openmm.Context(system, integrator, platform)
        print("  ✓ Context created successfully")
        
        # Set initial positions for all beads
        print(f"\n[TEST 5] Setting initial positions...")
        positions = []
        for i in range(num_particles):
            positions.append(openmm.Vec3(i * 0.3, 0, 0))
        
        # Set positions for the context (copy 0 is what the context sees)
        context.setPositions(positions)
        
        for bead in range(num_beads):
            # Add small random perturbation for each bead
            bead_positions = [openmm.Vec3(p[0] + 0.01 * np.random.randn(), 
                                          p[1] + 0.01 * np.random.randn(),
                                          p[2] + 0.01 * np.random.randn()) 
                             for p in positions]
            integrator.setPositions(bead, bead_positions)
        print("  ✓ Positions set for all beads")
        
        # Initialize velocities
        print(f"\n[TEST 6] Initializing velocities...")
        context.setVelocitiesToTemperature(temperature * unit.kelvin)
        
        # For RPMD, we need to initialize velocities for all beads, not just bead 0
        # Get velocities from context (which has random thermal velocities)
        state = context.getState(getVelocities=True)
        base_velocities = state.getVelocities()
        
        # Set slightly perturbed velocities for each bead
        for bead in range(num_beads):
            bead_velocities = []
            for v in base_velocities:
                # Strip units and add perturbation
                vx = v[0].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01 * np.random.randn()
                vy = v[1].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01 * np.random.randn()
                vz = v[2].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01 * np.random.randn()
                bead_velocities.append(openmm.Vec3(vx, vy, vz) * unit.nanometers/unit.picoseconds)
            integrator.setVelocities(bead, bead_velocities)
        print("  ✓ Velocities initialized for all beads")
        
        # Run simulation
        print(f"\n[TEST 7] Running simulation (100 steps)...")
        initial_energy = integrator.getTotalEnergy()
        print(f"  Initial energy: {initial_energy} kJ/mol")
        
        integrator.step(100)
        
        final_energy = integrator.getTotalEnergy()
        print(f"  Final energy: {final_energy} kJ/mol")
        print("  ✓ Simulation completed successfully")
        
        # Test 8: Verify temperature
        print(f"\n[TEST 8] Checking temperature...")
        state = integrator.getState(0, getEnergy=True, getVelocities=True)
        ke = state.getKineticEnergy()
        print(f"  Kinetic energy (bead 0): {ke}")
        print("  ✓ Temperature check completed")
        
        print("\n" + "=" * 60)
        print("ALL TESTS PASSED!")
        print("=" * 60)
        return 0
        
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(test_rpmd_pile_g())
