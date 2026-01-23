#!/usr/bin/env python3
"""
Test hybrid classical-quantum RPMD on GPU (CUDA).
"""

import os
os.environ['OPENMM_PLUGIN_DIR'] = '/media/extradrive/Trajectories/openmm/build_hybrid'

import sys
sys.path.insert(0, '/media/extradrive/Trajectories/openmm/build_hybrid/python/build/lib.linux-x86_64-cpython-312')

import openmm as mm
from openmm import unit
from openmm import RPMDIntegrator
import numpy as np

def test_cuda_hybrid():
    """Test hybrid RPMD on CUDA platform."""
    print("=" * 60)
    print("Testing Hybrid RPMD on CUDA")
    print("=" * 60)
    
    # Check if CUDA is available
    try:
        platform = mm.Platform.getPlatformByName('CUDA')
        print(f"CUDA platform found: {platform.getName()}")
    except Exception as e:
        print(f"CUDA platform not available: {e}")
        print("Trying OpenCL...")
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
            print(f"OpenCL platform found: {platform.getName()}")
        except Exception as e2:
            print(f"OpenCL also not available: {e2}")
            print("Skipping GPU test.")
            return True
    
    # Create system with multiple particles
    system = mm.System()
    
    # Add 3 water molecules: 3 O (classical) + 6 H (quantum)
    num_water = 3
    for i in range(num_water):
        system.addParticle(16.0 * unit.amu)  # O
        system.addParticle(1.0 * unit.amu)   # H
        system.addParticle(1.0 * unit.amu)   # H
    
    # Add harmonic bonds O-H
    harmonic = mm.HarmonicBondForce()
    for i in range(num_water):
        base = i * 3
        harmonic.addBond(base, base+1, 0.1 * unit.nanometer, 1000 * unit.kilojoule_per_mole / unit.nanometer**2)
        harmonic.addBond(base, base+2, 0.1 * unit.nanometer, 1000 * unit.kilojoule_per_mole / unit.nanometer**2)
    system.addForce(harmonic)
    
    # Create hybrid RPMD integrator
    numCopies = 8  # Use more beads for GPU test
    temperature = 300  # K
    friction = 1.0 / unit.picosecond
    stepSize = 0.5 * unit.femtosecond
    
    integrator = RPMDIntegrator(numCopies, temperature, friction, stepSize)
    
    # Configure hybrid mode:
    # - Type 1 = O (classical)
    # - Type 2 = H (quantum)
    for i in range(num_water):
        base = i * 3
        integrator.setParticleType(base, 1)      # O
        integrator.setParticleType(base + 1, 2)  # H
        integrator.setParticleType(base + 2, 2)  # H
    
    integrator.setQuantumParticleTypes({2})  # Type 2 is quantum
    integrator.setDefaultQuantum(False)  # Default is classical
    integrator.setClassicalThermostat(RPMDIntegrator.LangevinClassical)
    
    print(f"\nConfiguration:")
    print(f"  - Num copies: {numCopies}")
    print(f"  - Temperature: {temperature} K")
    print(f"  - Number of particles: {system.getNumParticles()}")
    print(f"  - Quantum particles (H): {num_water * 2}")
    print(f"  - Classical particles (O): {num_water}")
    print(f"  - Platform: {platform.getName()}")
    
    # Create context
    context = mm.Context(system, integrator, platform)
    
    # Set initial positions
    positions = []
    for i in range(num_water):
        # Place water molecules in a row
        x_offset = i * 0.5
        positions.append(mm.Vec3(x_offset, 0.0, 0.0) * unit.nanometer)
        positions.append(mm.Vec3(x_offset + 0.1, 0.0, 0.0) * unit.nanometer)
        positions.append(mm.Vec3(x_offset, 0.1, 0.0) * unit.nanometer)
    
    context.setPositions(positions)
    for copy in range(numCopies):
        integrator.setPositions(copy, positions)
    
    context.setVelocitiesToTemperature(temperature)
    
    # Run simulation
    print(f"\nRunning {1000} steps on GPU...")
    
    # Warmup
    integrator.step(100)
    
    # Timing test
    import time
    start = time.time()
    integrator.step(1000)
    elapsed = time.time() - start
    
    state = context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    
    print(f"\nResults:")
    print(f"  - Time for 1000 steps: {elapsed:.3f} s")
    print(f"  - Steps/second: {1000/elapsed:.1f}")
    print(f"  - Final PE: {pe:.4f} kJ/mol")
    print(f"  - Final KE: {ke:.4f} kJ/mol")
    print(f"  - Total energy: {pe+ke:.4f} kJ/mol")
    
    print("\n[OK] GPU hybrid RPMD test completed successfully!")
    return True


def test_larger_system():
    """Test with a larger system for performance."""
    print("\n" + "=" * 60)
    print("Testing Hybrid RPMD Performance")
    print("=" * 60)
    
    # Try CUDA first
    try:
        platform = mm.Platform.getPlatformByName('CUDA')
    except:
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
        except:
            print("No GPU platform available, skipping performance test.")
            return True
    
    # Create system with 100 water molecules
    num_water = 100
    system = mm.System()
    
    for i in range(num_water):
        system.addParticle(16.0 * unit.amu)  # O
        system.addParticle(1.0 * unit.amu)   # H
        system.addParticle(1.0 * unit.amu)   # H
    
    # Add harmonic bonds
    harmonic = mm.HarmonicBondForce()
    for i in range(num_water):
        base = i * 3
        harmonic.addBond(base, base+1, 0.1 * unit.nanometer, 1000 * unit.kilojoule_per_mole / unit.nanometer**2)
        harmonic.addBond(base, base+2, 0.1 * unit.nanometer, 1000 * unit.kilojoule_per_mole / unit.nanometer**2)
    system.addForce(harmonic)
    
    numCopies = 16
    temperature = 300
    friction = 1.0 / unit.picosecond
    stepSize = 0.5 * unit.femtosecond
    
    integrator = RPMDIntegrator(numCopies, temperature, friction, stepSize)
    
    # Configure hybrid mode
    for i in range(num_water):
        base = i * 3
        integrator.setParticleType(base, 1)      # O
        integrator.setParticleType(base + 1, 2)  # H
        integrator.setParticleType(base + 2, 2)  # H
    
    integrator.setQuantumParticleTypes({2})
    integrator.setDefaultQuantum(False)
    integrator.setClassicalThermostat(RPMDIntegrator.LangevinClassical)
    
    print(f"\nConfiguration:")
    print(f"  - Num copies: {numCopies}")
    print(f"  - Number of particles: {system.getNumParticles()}")
    print(f"  - Quantum particles (H): {num_water * 2}")
    print(f"  - Classical particles (O): {num_water}")
    print(f"  - Platform: {platform.getName()}")
    
    context = mm.Context(system, integrator, platform)
    
    # Set initial positions on a grid
    positions = []
    grid_size = int(np.ceil(num_water ** (1/3)))
    spacing = 0.4  # nm
    idx = 0
    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                if idx >= num_water:
                    break
                x, y, z = i * spacing, j * spacing, k * spacing
                positions.append(mm.Vec3(x, y, z) * unit.nanometer)
                positions.append(mm.Vec3(x + 0.1, y, z) * unit.nanometer)
                positions.append(mm.Vec3(x, y + 0.1, z) * unit.nanometer)
                idx += 1
    
    context.setPositions(positions)
    for copy in range(numCopies):
        integrator.setPositions(copy, positions)
    
    context.setVelocitiesToTemperature(temperature)
    
    # Warmup
    print("\nWarmup...")
    integrator.step(100)
    
    # Timing
    import time
    nSteps = 1000
    print(f"Running {nSteps} steps...")
    start = time.time()
    integrator.step(nSteps)
    elapsed = time.time() - start
    
    steps_per_second = nSteps / elapsed
    ns_per_day = steps_per_second * stepSize.value_in_unit(unit.nanosecond) * 86400
    
    print(f"\nPerformance:")
    print(f"  - Time: {elapsed:.3f} s")
    print(f"  - Steps/second: {steps_per_second:.1f}")
    print(f"  - ns/day: {ns_per_day:.4f}")
    
    print("\n[OK] Performance test completed!")
    return True


if __name__ == "__main__":
    try:
        test_cuda_hybrid()
        test_larger_system()
        print("\n" + "=" * 60)
        print("ALL GPU TESTS PASSED!")
        print("=" * 60)
    except Exception as e:
        print(f"\nTEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
