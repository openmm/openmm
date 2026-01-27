#!/usr/bin/env python3
"""
Test UMA+LES model in OpenMM-ML for molecular dynamics.

This script demonstrates how to use a trained UMA+LES model
in OpenMM for MD simulations with improved long-range electrostatics.

Usage:
    python test_uma_les_openmm.py --checkpoint /path/to/checkpoint.pt --molecules water
"""

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import openmm
    from openmm import app, unit
    from openmmml import MLPotential
    from ase import Atoms
except ImportError as e:
    print(f"Error: Missing dependency: {e}")
    print("Install with: pip install openmm openmmtorch openmmml ase")
    sys.exit(1)


def create_water_molecule():
    """Create a simple water molecule topology and positions."""
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("WAT", chain)
    
    O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    H1 = topology.addAtom("H1", app.Element.getBySymbol("H"), residue)
    H2 = topology.addAtom("H2", app.Element.getBySymbol("H"), residue)
    
    topology.addBond(O, H1)
    topology.addBond(O, H2)
    
    # Water geometry (Ångströms -> nm)
    positions = np.array([
        [0.0, 0.0, 0.0],      # O
        [0.0957, 0.0, 0.0],   # H1
        [-0.024, 0.093, 0.0]  # H2
    ]) * 0.1  # Å to nm
    
    return topology, positions * unit.nanometers


def create_methane_molecule():
    """Create a simple methane molecule topology and positions."""
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("MET", chain)
    
    C = topology.addAtom("C", app.Element.getBySymbol("C"), residue)
    H1 = topology.addAtom("H1", app.Element.getBySymbol("H"), residue)
    H2 = topology.addAtom("H2", app.Element.getBySymbol("H"), residue)
    H3 = topology.addAtom("H3", app.Element.getBySymbol("H"), residue)
    H4 = topology.addAtom("H4", app.Element.getBySymbol("H"), residue)
    
    topology.addBond(C, H1)
    topology.addBond(C, H2)
    topology.addBond(C, H3)
    topology.addBond(C, H4)
    
    # Tetrahedral geometry
    positions = np.array([
        [0.0, 0.0, 0.0],      # C
        [0.109, 0.109, 0.109],   # H1
        [-0.109, -0.109, 0.109], # H2
        [-0.109, 0.109, -0.109], # H3
        [0.109, -0.109, -0.109]  # H4
    ]) * 0.1  # Å to nm
    
    return topology, positions * unit.nanometers


def test_openmm_integration(checkpoint_path: str, molecule: str = "water", 
                           steps: int = 1000, temperature: float = 300.0):
    """
    Test trained UMA+LES model in OpenMM.
    
    Args:
        checkpoint_path: Path to trained model checkpoint
        molecule: Molecule to simulate ('water' or 'methane')
        steps: Number of MD steps
        temperature: Temperature in Kelvin
    """
    print("=" * 60)
    print("UMA+LES OpenMM Integration Test")
    print("=" * 60)
    print(f"Checkpoint: {checkpoint_path}")
    print(f"Molecule: {molecule}")
    print(f"Steps: {steps}")
    print(f"Temperature: {temperature} K")
    print("=" * 60)
    
    # Create molecule
    print(f"\n1. Creating {molecule} molecule...")
    if molecule.lower() == "water":
        topology, positions = create_water_molecule()
    elif molecule.lower() == "methane":
        topology, positions = create_methane_molecule()
    else:
        print(f"Error: Unknown molecule '{molecule}'")
        print("Supported: water, methane")
        return False
    
    print(f"   ✓ Topology created: {topology.getNumAtoms()} atoms")
    
    # Load UMA+LES model
    print("\n2. Loading UMA+LES model...")
    try:
        potential = MLPotential(checkpoint_path)
        print("   ✓ Model loaded")
    except Exception as e:
        print(f"   ✗ Failed to load model: {e}")
        return False
    
    # Create OpenMM system
    print("\n3. Creating OpenMM system...")
    try:
        system = potential.createSystem(
            topology,
            task_name='omol',
            charge=0,
            spin=1
        )
        print(f"   ✓ System created: {system.getNumParticles()} particles")
    except Exception as e:
        print(f"   ✗ Failed to create system: {e}")
        return False
    
    # Set up integrator
    print(f"\n4. Setting up Langevin integrator (T={temperature}K)...")
    integrator = openmm.LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtoseconds
    )
    
    # Create simulation
    print("\n5. Creating simulation context...")
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed'}
        print("   ✓ Using CUDA platform")
    except Exception:
        try:
            platform = openmm.Platform.getPlatformByName('CPU')
            properties = {}
            print("   ✓ Using CPU platform (CUDA not available)")
        except Exception as e:
            print(f"   ✗ Failed to get platform: {e}")
            return False
    
    simulation = app.Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
    
    # Get initial state
    print("\n6. Computing initial state...")
    state = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
    initial_energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    initial_forces = state.getForces(asNumpy=True)
    max_force = np.max(np.linalg.norm(initial_forces, axis=1))
    
    print(f"   Initial energy: {initial_energy:.4f} kcal/mol")
    print(f"   Max force: {max_force:.4f} kJ/(mol·nm)")
    
    # Run MD
    print(f"\n7. Running {steps} MD steps...")
    try:
        integrator.step(steps)
        print(f"   ✓ MD completed")
    except Exception as e:
        print(f"   ✗ MD failed: {e}")
        return False
    
    # Get final state
    print("\n8. Computing final state...")
    state = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
    final_energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    final_forces = state.getForces(asNumpy=True)
    max_force = np.max(np.linalg.norm(final_forces, axis=1))
    final_positions = state.getPositions(asNumpy=True)
    
    print(f"   Final energy: {final_energy:.4f} kcal/mol")
    print(f"   Max force: {max_force:.4f} kJ/(mol·nm)")
    print(f"   Energy change: {final_energy - initial_energy:.4f} kcal/mol")
    
    # Check RMSD
    initial_pos = positions.value_in_unit(unit.nanometers)
    final_pos = final_positions.value_in_unit(unit.nanometers)
    rmsd = np.sqrt(np.mean((initial_pos - final_pos)**2)) * 10  # nm to Å
    print(f"   RMSD: {rmsd:.4f} Å")
    
    print("\n" + "=" * 60)
    print("✓ OpenMM integration test completed successfully!")
    print("=" * 60)
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Test UMA+LES model in OpenMM-ML",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--checkpoint", "-c",
        type=str,
        required=True,
        help="Path to trained model checkpoint (.pt file)"
    )
    
    parser.add_argument(
        "--molecule", "-m",
        type=str,
        default="water",
        choices=["water", "methane"],
        help="Molecule to simulate"
    )
    
    parser.add_argument(
        "--steps", "-s",
        type=int,
        default=1000,
        help="Number of MD steps"
    )
    
    parser.add_argument(
        "--temperature", "-T",
        type=float,
        default=300.0,
        help="Temperature in Kelvin"
    )
    
    args = parser.parse_args()
    
    # Validate checkpoint
    if not Path(args.checkpoint).exists():
        print(f"Error: Checkpoint file not found: {args.checkpoint}")
        sys.exit(1)
    
    # Run test
    success = test_openmm_integration(
        checkpoint_path=args.checkpoint,
        molecule=args.molecule,
        steps=args.steps,
        temperature=args.temperature
    )
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
