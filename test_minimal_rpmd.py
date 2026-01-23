#!/usr/bin/env python3
"""
Minimal RPMD test to isolate double-free issue.
"""
import sys
import openmm
from openmm import unit

print("Creating system...")
system = openmm.System()
for i in range(2):
    system.addParticle(1.0)

print("Creating integrator...")
integrator = openmm.RPMDIntegrator(2, 300*unit.kelvin, 1.0/unit.picosecond, 0.5*unit.femtosecond)

print("Creating context...")
try:
    context = openmm.Context(system, integrator)
    print("✓ Context created")
    del context
    print("✓ Context deleted")
except Exception as e:
    print(f"✗ Error: {e}")
    sys.exit(1)

print("✓ Test completed successfully")
sys.exit(0)
