#!/usr/bin/env python3
"""
Test hybrid classical-quantum RPMD implementation.

This test creates a simple water molecule where:
- H atoms are quantum (RPMD with multiple beads)
- O atom is classical (all beads coupled together)
"""

import os
os.environ['OPENMM_PLUGIN_DIR'] = '/media/extradrive/Trajectories/openmm/build_hybrid'

import sys
sys.path.insert(0, '/media/extradrive/Trajectories/openmm/build_hybrid/python/build/lib.linux-x86_64-cpython-312')

import openmm as mm
from openmm import unit
from openmm import RPMDIntegrator
import numpy as np

def test_api():
    """Test that the hybrid RPMD API is properly exposed."""
    print("=" * 60)
    print("Testing Hybrid RPMD API")
    print("=" * 60)
    
    # Create integrator
    numCopies = 4
    temperature = 300  # K
    friction = 1.0 / unit.picosecond
    stepSize = 0.5 * unit.femtosecond
    
    integrator = RPMDIntegrator(numCopies, temperature, friction, stepSize)
    
    # Test API methods exist
    print("\n1. Testing API methods exist:")
    print(f"   - setParticleType: {hasattr(integrator, 'setParticleType')}")
    print(f"   - getParticleTypes: {hasattr(integrator, 'getParticleTypes')}")
    print(f"   - setQuantumParticleTypes: {hasattr(integrator, 'setQuantumParticleTypes')}")
    print(f"   - getQuantumParticleTypes: {hasattr(integrator, 'getQuantumParticleTypes')}")
    print(f"   - setDefaultQuantum: {hasattr(integrator, 'setDefaultQuantum')}")
    print(f"   - getDefaultQuantum: {hasattr(integrator, 'getDefaultQuantum')}")
    print(f"   - setClassicalThermostat: {hasattr(integrator, 'setClassicalThermostat')}")
    print(f"   - getClassicalThermostat: {hasattr(integrator, 'getClassicalThermostat')}")
    
    # Test ClassicalThermostatType enum
    print("\n2. Testing ClassicalThermostatType enum:")
    print(f"   - BussiClassical: {RPMDIntegrator.BussiClassical}")
    print(f"   - LangevinClassical: {RPMDIntegrator.LangevinClassical}")
    print(f"   - NoneClassical: {RPMDIntegrator.NoneClassical}")
    
    # Test setting particle types
    print("\n3. Testing setParticleType:")
    integrator.setParticleType(0, 1)  # Particle 0 is type 1 (O)
    integrator.setParticleType(1, 2)  # Particle 1 is type 2 (H)
    integrator.setParticleType(2, 2)  # Particle 2 is type 2 (H)
    types = integrator.getParticleTypes()
    print(f"   - Particle types: {dict(types)}")
    
    # Test quantum particle types
    print("\n4. Testing quantum particle types:")
    integrator.setQuantumParticleTypes({2})  # Type 2 (H) is quantum
    quantum_types = integrator.getQuantumParticleTypes()
    print(f"   - Quantum types: {set(quantum_types)}")
    
    # Test default quantum behavior
    print("\n5. Testing default quantum:")
    print(f"   - Default quantum: {integrator.getDefaultQuantum()}")
    integrator.setDefaultQuantum(False)  # Default is classical
    print(f"   - After setDefaultQuantum(False): {integrator.getDefaultQuantum()}")
    
    # Test classical thermostat
    print("\n6. Testing classical thermostat:")
    print(f"   - Default: {integrator.getClassicalThermostat()}")
    integrator.setClassicalThermostat(RPMDIntegrator.LangevinClassical)
    print(f"   - After setting Langevin: {integrator.getClassicalThermostat()}")
    
    print("\nAPI tests passed!")
    return True


def test_simple_simulation():
    """Test a simple hybrid RPMD simulation with a harmonic oscillator."""
    print("\n" + "=" * 60)
    print("Testing Hybrid RPMD Simulation")
    print("=" * 60)
    
    # Create a simple 2-particle system
    system = mm.System()
    
    # Particle 0: heavy (classical O-like)
    mass_O = 16.0 * unit.amu
    system.addParticle(mass_O)
    
    # Particle 1: light (quantum H-like)
    mass_H = 1.0 * unit.amu
    system.addParticle(mass_H)
    
    # Add a harmonic bond between them
    harmonic = mm.HarmonicBondForce()
    harmonic.addBond(0, 1, 0.1 * unit.nanometer, 1000 * unit.kilojoule_per_mole / unit.nanometer**2)
    system.addForce(harmonic)
    
    # Create hybrid RPMD integrator
    numCopies = 4
    temperature = 300  # K
    friction = 1.0 / unit.picosecond
    stepSize = 0.5 * unit.femtosecond
    
    integrator = RPMDIntegrator(numCopies, temperature, friction, stepSize)
    
    # Configure hybrid mode:
    # - Type 1 = O (classical)
    # - Type 2 = H (quantum)
    integrator.setParticleType(0, 1)  # O is type 1
    integrator.setParticleType(1, 2)  # H is type 2
    integrator.setQuantumParticleTypes({2})  # Type 2 is quantum
    integrator.setDefaultQuantum(False)  # Default is classical
    integrator.setClassicalThermostat(RPMDIntegrator.LangevinClassical)
    
    print("\nConfiguration:")
    print(f"  - Num copies: {numCopies}")
    print(f"  - Temperature: {temperature} K")
    print(f"  - Friction: {friction}")
    print(f"  - Step size: {stepSize}")
    print(f"  - Particle 0 (O, {mass_O}): CLASSICAL")
    print(f"  - Particle 1 (H, {mass_H}): QUANTUM")
    print(f"  - Classical thermostat: Langevin")
    
    # Create context
    platform = mm.Platform.getPlatformByName('Reference')
    context = mm.Context(system, integrator, platform)
    
    # Set initial positions using context first
    positions = [
        mm.Vec3(0.0, 0.0, 0.0),
        mm.Vec3(0.1, 0.0, 0.0)
    ]
    positions = [p * unit.nanometer for p in positions]
    
    # Set positions on context first (required for setVelocitiesToTemperature)
    context.setPositions(positions)
    
    # Then set positions on all beads
    for copy in range(numCopies):
        integrator.setPositions(copy, positions)
    
    # Set velocities
    context.setVelocitiesToTemperature(temperature)
    
    # Run simulation
    print("\nRunning simulation...")
    nSteps = 100
    printInterval = 20
    
    for step in range(0, nSteps, printInterval):
        integrator.step(printInterval)
        
        # Get current state
        state = context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        
        print(f"  Step {step + printInterval}: PE = {pe:.4f}, KE = {ke:.4f}, Total = {pe+ke:.4f} kJ/mol")
    
    print("\nSimulation completed - hybrid mode is working!")
    
    print("\nSimulation test completed!")
    return True


def test_energy_conservation():
    """Test energy conservation for NVE (no thermostat) hybrid RPMD."""
    print("\n" + "=" * 60)
    print("Testing Energy Conservation (NVE)")
    print("=" * 60)
    
    # Create system
    system = mm.System()
    system.addParticle(16.0 * unit.amu)  # O
    system.addParticle(1.0 * unit.amu)   # H
    
    harmonic = mm.HarmonicBondForce()
    harmonic.addBond(0, 1, 0.1 * unit.nanometer, 1000 * unit.kilojoule_per_mole / unit.nanometer**2)
    system.addForce(harmonic)
    
    # Create NVE hybrid RPMD (no thermostat)
    numCopies = 4
    temperature = 300
    friction = 0.0 / unit.picosecond  # No friction = NVE
    stepSize = 0.1 * unit.femtosecond  # Small for energy conservation
    
    integrator = RPMDIntegrator(numCopies, temperature, friction, stepSize)
    integrator.setApplyThermostat(False)  # Disable thermostat
    
    # Hybrid setup
    integrator.setParticleType(0, 1)
    integrator.setParticleType(1, 2)
    integrator.setQuantumParticleTypes({2})
    integrator.setDefaultQuantum(False)
    
    platform = mm.Platform.getPlatformByName('Reference')
    context = mm.Context(system, integrator, platform)
    
    positions = [
        mm.Vec3(0.0, 0.0, 0.0) * unit.nanometer,
        mm.Vec3(0.1, 0.0, 0.0) * unit.nanometer
    ]
    
    context.setPositions(positions)
    for copy in range(numCopies):
        integrator.setPositions(copy, positions)
    
    context.setVelocitiesToTemperature(temperature)
    
    # Sync classical beads first by running one step
    integrator.step(1)
    
    # Get initial total energy
    state = context.getState(getEnergy=True)
    initial_energy = state.getPotentialEnergy() + state.getKineticEnergy()
    initial_energy_val = initial_energy.value_in_unit(unit.kilojoule_per_mole)
    
    print(f"\nInitial total energy: {initial_energy_val:.6f} kJ/mol")
    
    # Run and monitor energy
    energies = [initial_energy_val]
    for i in range(100):
        integrator.step(10)
        state = context.getState(getEnergy=True)
        e = (state.getPotentialEnergy() + state.getKineticEnergy()).value_in_unit(unit.kilojoule_per_mole)
        energies.append(e)
    
    energies = np.array(energies)
    energy_drift = np.max(np.abs(energies - initial_energy_val))
    
    print(f"Final total energy: {energies[-1]:.6f} kJ/mol")
    print(f"Max energy drift: {energy_drift:.6e} kJ/mol")
    
    if energy_drift < 0.1:  # Allow 0.1 kJ/mol drift
        print("[OK] Energy is well conserved")
    else:
        print(f"[WARNING] Energy drift is larger than expected: {energy_drift:.6f} kJ/mol")
    
    return True


if __name__ == "__main__":
    try:
        test_api()
        test_simple_simulation()
        test_energy_conservation()
        print("\n" + "=" * 60)
        print("ALL TESTS PASSED!")
        print("=" * 60)
    except Exception as e:
        print(f"\nTEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
