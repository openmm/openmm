#!/usr/bin/env python3
"""
Test UMA with RPMD integrator using PythonForce.
"""

import sys
import openmm
from openmm import app, unit
import numpy as np

print("="*70)
print("UMA + PythonForce + RPMD Test")
print("="*70)

try:
    from openmmml import MLPotential
    from openmm import RPMDIntegrator
    
    # Create water molecule
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("WAT", chain)
    O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    H1 = topology.addAtom("H1", app.Element.getBySymbol("H"), residue)
    H2 = topology.addAtom("H2", app.Element.getBySymbol("H"), residue)
    topology.addBond(O, H1)
    topology.addBond(O, H2)
    
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.0957, 0.0, 0.0],
        [-0.024, 0.093, 0.0]
    ])
    
    # Create UMA potential
    print("\nLoading UMA model with PythonForce...")
    potential = MLPotential('uma-s-1p1-pythonforce')
    
    print("Creating system...")
    system = potential.createSystem(
        topology,
        task_name='omol',
        charge=0,
        spin=1
    )
    
    # Create RPMD integrator
    print("Creating RPMD integrator with 4 beads...")
    integrator = RPMDIntegrator(
        4,  # numCopies
        300*unit.kelvin,  # temperature
        1.0/unit.picosecond,  # frictionCoeff
        0.5*unit.femtoseconds  # stepSize
    )
    # Default is PILE thermostat, which is fine
    
    # Create context
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        print(f"Using CUDA platform")
    except:
        platform = openmm.Platform.getPlatformByName('CPU')
        print(f"Using CPU platform")
    
    context = openmm.Context(system, integrator, platform)
    
    # Initialize all beads with same positions
    print("Initializing RPMD beads...")
    for i in range(integrator.getNumCopies()):
        integrator.setPositions(i, positions * unit.nanometers)
        # Add small perturbation to break symmetry
        pert = np.random.randn(3, 3) * 0.001
        integrator.setPositions(i, (positions + pert) * unit.nanometers)
    
    # Check initial energies
    print("\nInitial bead energies:")
    for i in range(integrator.getNumCopies()):
        state = integrator.getState(i, getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"  Bead {i}: {energy:.2f} kJ/mol")
        assert np.isfinite(energy), f"Non-finite energy in bead {i}"
    
    # Run RPMD simulation
    print("\nRunning 50 RPMD steps...")
    integrator.step(50)
    
    # Check final energies
    print("\nFinal bead energies:")
    energies = []
    for i in range(integrator.getNumCopies()):
        state = integrator.getState(i, getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        energies.append(energy)
        print(f"  Bead {i}: {energy:.2f} kJ/mol")
        assert np.isfinite(energy), f"Non-finite energy in bead {i}"
    
    # Statistics
    mean_energy = np.mean(energies)
    std_energy = np.std(energies)
    print(f"\nEnergy statistics:")
    print(f"  Mean: {mean_energy:.2f} kJ/mol")
    print(f"  Std:  {std_energy:.2f} kJ/mol")
    
    print("\n" + "="*70)
    print("✓ RPMD TEST PASSED!")
    print("="*70)
    print("\nUMA + PythonForce + RPMD is working!")
    print("You now have GPU-accelerated RPMD with UMA models in OpenMM.")
    
    sys.exit(0)
    
except Exception as e:
    print(f"\n✗ Test FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
