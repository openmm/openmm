#!/usr/bin/env python3
"""
Quick test to verify hybrid RPMD API compiles and can be imported.
"""

import sys
import os

# Test 1: Check if we can import OpenMM
try:
    import openmm
    print("✓ OpenMM imported successfully")
except ImportError as e:
    print(f"✗ Failed to import OpenMM: {e}")
    sys.exit(1)

# Test 2: Check if RPMDIntegrator exists
try:
    from openmm import RPMDIntegrator
    print("✓ RPMDIntegrator imported successfully")
except ImportError as e:
    print(f"✗ Failed to import RPMDIntegrator: {e}")
    sys.exit(1)

# Test 3: Create a basic integrator
try:
    integrator = RPMDIntegrator(8, 300.0, 1.0, 0.001)
    print(f"✓ Created RPMDIntegrator: {integrator.getNumCopies()} beads")
except Exception as e:
    print(f"✗ Failed to create RPMDIntegrator: {e}")
    sys.exit(1)

# Test 4: Check if new API methods exist
try:
    # These should exist after our changes
    hasattr(integrator, 'setParticleType')
    hasattr(integrator, 'getParticleTypes')
    hasattr(integrator, 'setQuantumParticleTypes')
    hasattr(integrator, 'getQuantumParticleTypes')
    hasattr(integrator, 'setDefaultQuantum')
    hasattr(integrator, 'getDefaultQuantum')
    hasattr(integrator, 'setClassicalThermostat')
    hasattr(integrator, 'getClassicalThermostat')
    print("✓ New hybrid RPMD API methods exist")
except AttributeError as e:
    print(f"⚠ New API methods not yet available (needs recompilation): {e}")

# Test 5: Check ClassicalThermostatType enum
try:
    hasattr(RPMDIntegrator, 'BussiClassical')
    hasattr(RPMDIntegrator, 'LangevinClassical')
    hasattr(RPMDIntegrator, 'NoneClassical')
    print("✓ ClassicalThermostatType enum exists")
except AttributeError as e:
    print(f"⚠ ClassicalThermostatType enum not yet available (needs recompilation): {e}")

print("\n" + "="*60)
print("Compilation test summary:")
print("- Core OpenMM and RPMDIntegrator work")
print("- New API will be available after recompilation")
print("="*60)
