#!/usr/bin/env python3
"""
Quick test to verify eSEN model works with OpenMM RPMD batching.
"""

import sys
import numpy as np
from openmm import app, unit
from openmm import RPMDIntegrator, Context, Platform
from openmmml import MLPotential

print("="*80)
print("eSEN RPMD Batching Test")
print("="*80)

# Create a small water system
print("\nCreating 10 water molecules...")
from ase.build import molecule
from ase import Atoms

positions = []
symbols = []
for i in range(10):
    base_x = i * 0.3  # nm spacing
    positions.extend([
        [base_x, 0, 0],           # O
        [base_x + 0.09572, 0, 0], # H1
        [base_x - 0.024, 0.093, 0] # H2
    ])
    symbols.extend(['O', 'H', 'H'])

# Create topology
topology = app.Topology()
chain = topology.addChain()
for i in range(10):
    residue = topology.addResidue('HOH', chain)
    o_atom = topology.addAtom('O', app.Element.getBySymbol('O'), residue)
    h1_atom = topology.addAtom('H', app.Element.getBySymbol('H'), residue)
    h2_atom = topology.addAtom('H', app.Element.getBySymbol('H'), residue)
    topology.addBond(o_atom, h1_atom)
    topology.addBond(o_atom, h2_atom)

box_vectors = ([3.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0])
topology.setPeriodicBoxVectors([v*unit.nanometer for v in box_vectors])

# Test with eSEN model
print("\nLoading eSEN model...")
potential = MLPotential('esen-sm-conserving-all-omol-pythonforce-batch')

print("Creating system...")
system = potential.createSystem(
    topology,
    task_name='omol',
    charge=0,
    spin=1
)

print(f"System PBC: {system.usesPeriodicBoundaryConditions()}")

# Create RPMD integrator with 4 beads
num_beads = 4
print(f"\nCreating RPMD integrator with {num_beads} beads...")
integrator = RPMDIntegrator(
    num_beads,
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    1.0 * unit.femtoseconds
)

# Create context
print("Creating context...")
platform = Platform.getPlatformByName('CUDA')
properties = {'Precision': 'mixed'}
context = Context(system, integrator, platform, properties)

# Initialize bead positions
print(f"Initializing {num_beads} beads...")
positions_nm = np.array(positions)
for i in range(num_beads):
    perturbed = positions_nm + np.random.randn(*positions_nm.shape) * 0.001
    integrator.setPositions(i, perturbed * unit.nanometer)
    integrator.setVelocities(i, np.zeros_like(perturbed) * unit.nanometer / unit.picosecond)

# Get initial energies (this will trigger batched evaluation!)
print("\nGetting initial bead energies (triggers batched evaluation)...")
try:
    for i in range(num_beads):
        state = integrator.getState(i, getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"  Bead {i}: E = {energy:.2f} kJ/mol")
    
    print("\nSUCCESS! eSEN model works with batched RPMD!")
    print("\nRunning a few MD steps...")
    integrator.step(5)
    print("MD steps completed successfully!")
    
    sys.exit(0)
    
except Exception as e:
    print(f"\nFAILED: {type(e).__name__}: {str(e)}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
