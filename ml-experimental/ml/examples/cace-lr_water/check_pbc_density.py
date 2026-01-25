#!/usr/bin/env python3
"""
Check if periodic boundary conditions are working correctly
and verify the density of the water system.
"""

import sys
import numpy as np
from pathlib import Path
sys.path.append(str(Path(__file__).parent))

import openmm
from openmm import unit
from openmmml import MLPotential
from run_water_cace_lr import create_water_topology

def check_system(num_molecules=32):
    """Check PBC and density of the system."""
    print("=" * 60)
    print("Checking Periodic Boundary Conditions and Density")
    print("=" * 60)
    
    # Create topology
    topology, positions, box_size_nm = create_water_topology(num_molecules)
    
    # Set periodic box vectors on topology
    topology.setPeriodicBoxVectors((
        [box_size_nm, 0, 0] * unit.nanometer,
        [0, box_size_nm, 0] * unit.nanometer,
        [0, 0, box_size_nm] * unit.nanometer
    ))
    
    print(f"\n--- System Setup ---")
    print(f"  Molecules: {num_molecules}")
    print(f"  Box size: {box_size_nm:.4f} nm")
    
    # Calculate volume and density
    volume_nm3 = box_size_nm ** 3
    volume_cm3 = volume_nm3 * 1e-21  # nm³ to cm³
    mass_g = num_molecules * 18.01528 / 6.022e23  # g per molecule
    density_gcm3 = mass_g / volume_cm3
    
    print(f"  Volume: {volume_nm3:.6f} nm³ = {volume_cm3*1e21:.6f} Å³")
    print(f"  Mass: {mass_g*1e3:.6f} g")
    print(f"  Density: {density_gcm3:.4f} g/cm³")
    print(f"  Expected (liquid water): ~1.0 g/cm³")
    
    if density_gcm3 < 0.5:
        print(f"\n⚠ WARNING: Density is too low! This will behave like a gas.")
        print(f"  Consider reducing box size or increasing number of molecules.")
    elif density_gcm3 > 1.2:
        print(f"\n⚠ WARNING: Density is too high! This may cause issues.")
    else:
        print(f"\n✓ Density is reasonable for liquid water.")
    
    # Create system and check PBC
    model_path = "/media/extradrive/Trajectories/openmm/LES-BEC/water/fit/fit_version_1/best_model.pth"
    potential = MLPotential('cace-lr', model_path=model_path)
    system = potential.createSystem(topology)
    
    # Set box vectors
    system.setDefaultPeriodicBoxVectors(
        [box_size_nm, 0, 0] * unit.nanometer,
        [0, box_size_nm, 0] * unit.nanometer,
        [0, 0, box_size_nm] * unit.nanometer
    )
    
    # Check if system uses PBC
    uses_pbc = system.usesPeriodicBoundaryConditions()
    print(f"\n--- Periodic Boundary Conditions ---")
    print(f"  System uses PBC: {uses_pbc}")
    
    # Get box vectors from system
    a, b, c = system.getDefaultPeriodicBoxVectors()
    print(f"  Box vector a: {a}")
    print(f"  Box vector b: {b}")
    print(f"  Box vector c: {c}")
    
    # Check topology
    top_box = topology.getPeriodicBoxVectors()
    if top_box is not None:
        print(f"  Topology box vectors: {top_box}")
    else:
        print(f"  Topology box vectors: None (this is OK if set on system)")
    
    # Check forces
    print(f"\n--- Forces in System ---")
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force_pbc = force.usesPeriodicBoundaryConditions()
        print(f"  Force {i} ({type(force).__name__}): uses PBC = {force_pbc}")
    
    # Create a context to check actual state
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.0005*unit.picosecond)
    platform = openmm.Platform.getPlatformByName('Reference')
    context = openmm.Context(system, integrator, platform)
    
    # Set positions
    pos_with_units = [openmm.Vec3(p[0], p[1], p[2]) * unit.nanometer for p in positions]
    context.setPositions(pos_with_units)
    
    # Get state and check box vectors
    state = context.getState()
    state_box = state.getPeriodicBoxVectors(asNumpy=True)
    print(f"\n--- Context State ---")
    print(f"  State box vectors:")
    print(f"    a: {state_box[0]}")
    print(f"    b: {state_box[1]}")
    print(f"    c: {state_box[2]}")
    
    # Calculate minimum image distance
    # For a cubic box, the cutoff should be < box_size/2
    cutoff = 5.5  # Angstrom, from CACE model
    cutoff_nm = cutoff / 10.0
    print(f"\n--- Neighbor List Check ---")
    print(f"  CACE cutoff: {cutoff} Å = {cutoff_nm:.4f} nm")
    print(f"  Box size: {box_size_nm:.4f} nm")
    print(f"  Box size / 2: {box_size_nm/2:.4f} nm")
    
    if cutoff_nm > box_size_nm / 2:
        print(f"  ⚠ WARNING: Cutoff ({cutoff_nm:.4f} nm) > box_size/2 ({box_size_nm/2:.4f} nm)")
        print(f"     This means each atom will interact with multiple periodic images!")
        print(f"     This is OK but may cause issues with neighbor lists.")
    else:
        print(f"  ✓ Cutoff is less than box_size/2, so each atom sees at most one image of each neighbor")
    
    # Check actual distances between molecules
    print(f"\n--- Intermolecular Distances ---")
    n_atoms = len(positions)
    n_mols = num_molecules
    atoms_per_mol = n_atoms // n_mols
    
    # Get center of mass of each molecule
    mol_coms = []
    for i in range(n_mols):
        start_idx = i * atoms_per_mol
        mol_pos = np.array(positions[start_idx:start_idx+atoms_per_mol])
        com = mol_pos.mean(axis=0)
        mol_coms.append(com)
    
    mol_coms = np.array(mol_coms)
    
    # Calculate minimum image distances
    min_distances = []
    for i in range(min(10, n_mols)):  # Check first 10 molecules
        for j in range(i+1, min(10, n_mols)):
            dr = mol_coms[j] - mol_coms[i]
            # Apply minimum image convention
            dr = dr - box_size_nm * np.round(dr / box_size_nm)
            dist = np.linalg.norm(dr)
            min_distances.append(dist)
    
    if min_distances:
        avg_dist = np.mean(min_distances)
        min_dist = np.min(min_distances)
        print(f"  Average minimum image distance (first 10 molecules): {avg_dist:.4f} nm")
        print(f"  Minimum distance: {min_dist:.4f} nm")
        print(f"  Typical O-O distance in liquid water: ~0.28 nm")
        
        if avg_dist > 0.5:
            print(f"  ⚠ WARNING: Molecules are too far apart! This looks like a gas.")
        elif avg_dist < 0.2:
            print(f"  ⚠ WARNING: Molecules are too close! May cause overlap.")
        else:
            print(f"  ✓ Intermolecular distances look reasonable for liquid water.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--molecules', type=int, default=32)
    args = parser.parse_args()
    
    check_system(args.molecules)
