#!/usr/bin/env python3
"""
Baseline Water MD Simulation using UMA ML Potential
=====================================================

System: Water molecules using FAIRChem UMA potential
Purpose: Validate IR spectrum with ML potential before cavity experiments
Protocol: NO cavity - pure water baseline for spectroscopic validation

UMA (Universal Machine Learning Architecture) captures all interactions
including bonds, angles, and nonbonded forces through neural networks
trained on quantum mechanical data.

Expected IR Peaks:
- OH stretch: ~3200-3700 cm⁻¹ (broad, strong)
- HOH bending: ~1645 cm⁻¹
- Librational modes: 400-900 cm⁻¹
"""

import sys
import numpy as np
import time
from pathlib import Path

try:
    import openmm
    from openmm import unit
    from openmm.app import Topology, Element
    from openmmml import MLPotential
    import torch
    print("OpenMM and OpenMM-ML loaded successfully")
    print(f"  PyTorch version: {torch.__version__}")
    print(f"  CUDA available: {torch.cuda.is_available()}")
except ImportError as e:
    print(f"Error importing required packages: {e}")
    print("Make sure OpenMM, openmmml, and fairchem are installed.")
    sys.exit(1)


def create_water_topology(num_molecules=10):
    """
    Create a water box topology for UMA simulation.
    
    Parameters
    ----------
    num_molecules : int
        Number of water molecules
        
    Returns
    -------
    topology : openmm.app.Topology
        OpenMM topology with water molecules
    positions : list of lists
        Atomic positions in nm
    """
    print(f"\n--- Creating Water Topology ---")
    print(f"  Model: UMA ML Potential (no explicit force field)")
    print(f"  Molecules: {num_molecules}")
    
    # Water geometry
    oh_bond = 0.09572  # nm
    hoh_angle = 104.52 * np.pi / 180.0  # radians
    
    # Calculate box size (density ~1 g/cm³)
    # 1 water = 18 g/mol, 1 g/cm³ -> 1 mol/18 cm³ -> 18 cm³/mol
    # 18 cm³/mol / (6.022e23 mol⁻¹) = 2.99e-23 cm³ = 29.9 Å³ per molecule
    vol_per_mol_nm3 = 0.0299  # nm³ per water molecule
    total_vol = num_molecules * vol_per_mol_nm3
    box_size_nm = total_vol ** (1/3)
    # Add some extra space
    box_size_nm = max(box_size_nm * 1.2, 1.0)
    
    print(f"  Box size: {box_size_nm:.2f} nm")
    
    # Create topology
    topology = Topology()
    positions = []
    
    # Arrange molecules in a grid
    side = int(np.ceil(num_molecules**(1/3)))
    spacing = box_size_nm / side
    
    mol_count = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if mol_count >= num_molecules:
                    break
                    
                # Add residue
                chain = topology.addChain()
                residue = topology.addResidue('HOH', chain)
                
                # Add atoms
                o_atom = topology.addAtom('O', Element.getBySymbol('O'), residue)
                h1_atom = topology.addAtom('H', Element.getBySymbol('H'), residue)
                h2_atom = topology.addAtom('H', Element.getBySymbol('H'), residue)
                
                # Center position
                cx = (i + 0.5) * spacing
                cy = (j + 0.5) * spacing
                cz = (k + 0.5) * spacing
                
                # Random rotation for each molecule
                np.random.seed(42 + mol_count)
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                psi = np.random.rand() * 2 * np.pi
                
                # Water molecule geometry
                h1_local = np.array([oh_bond * np.sin(hoh_angle/2), 0, oh_bond * np.cos(hoh_angle/2)])
                h2_local = np.array([-oh_bond * np.sin(hoh_angle/2), 0, oh_bond * np.cos(hoh_angle/2)])
                o_local = np.array([0, 0, 0])
                
                # Rotation matrices
                Rz = np.array([
                    [np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]
                ])
                Ry = np.array([
                    [np.cos(phi), 0, np.sin(phi)],
                    [0, 1, 0],
                    [-np.sin(phi), 0, np.cos(phi)]
                ])
                Rx = np.array([
                    [1, 0, 0],
                    [0, np.cos(psi), -np.sin(psi)],
                    [0, np.sin(psi), np.cos(psi)]
                ])
                R = Rz @ Ry @ Rx
                
                # Rotate and translate
                o_pos = R @ o_local + np.array([cx, cy, cz])
                h1_pos = R @ h1_local + np.array([cx, cy, cz])
                h2_pos = R @ h2_local + np.array([cx, cy, cz])
                
                positions.append(o_pos.tolist())
                positions.append(h1_pos.tolist())
                positions.append(h2_pos.tolist())
                
                mol_count += 1
                
            if mol_count >= num_molecules:
                break
        if mol_count >= num_molecules:
            break
    
    print(f"  Created {mol_count} water molecules ({len(positions)} atoms)")
    
    return topology, positions, box_size_nm


def compute_dipole(positions, charges):
    """Compute total molecular dipole moment."""
    dipole = np.zeros(3)
    for i, pos in enumerate(positions):
        dipole += charges[i] * np.array(pos)
    return dipole


def run_baseline_simulation(num_molecules=10, temperature_K=300.0,
                           dt_ps=0.0005, equilibration_ps=10.0, production_ps=100.0,
                           output_prefix="uma_water_baseline"):
    """
    Run baseline water MD with UMA ML potential - NO CAVITY.
    
    This validates that:
    1. UMA potential works for water
    2. IR spectrum shows expected peaks
    3. System is numerically stable
    
    Parameters
    ----------
    num_molecules : int
        Number of water molecules
    temperature_K : float
        Temperature in K
    dt_ps : float
        Timestep in ps (MUST be small for flexible bonds, recommend 0.5 fs)
    equilibration_ps : float
        Equilibration time in ps
    production_ps : float
        Production time in ps
    output_prefix : str
        Prefix for output files
    """
    print("=" * 80)
    print("Baseline Water MD Simulation - UMA ML Potential")
    print("=" * 80)
    print("\nPurpose: Validate IR spectrum with ML potential")
    print("Expected: OH stretch peak at 3200-3700 cm⁻¹")
    
    print(f"\nSimulation Parameters:")
    print(f"  Water molecules: {num_molecules}")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt_ps} ps ({dt_ps*1000:.2f} fs)")
    print(f"  Equilibration: {equilibration_ps} ps")
    print(f"  Production: {production_ps} ps")
    
    # Calculate steps
    equilibration_steps = int(equilibration_ps / dt_ps)
    production_steps = int(production_ps / dt_ps)
    total_steps = equilibration_steps + production_steps
    
    # Create topology
    np.random.seed(42)
    topology, positions, box_size_nm = create_water_topology(num_molecules)
    
    # Create UMA potential
    print(f"\n--- Creating UMA Potential ---")
    print(f"  Loading uma-s-1p1 model...")
    potential = MLPotential('uma-s-1p1-pythonforce')
    
    # Create system
    print(f"\n--- Creating OpenMM System ---")
    system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
    print(f"  System created with {system.getNumParticles()} particles")
    print(f"  Forces: {[type(f).__name__ for f in system.getForces()]}")
    
    # Charges for dipole calculation (standard water charges)
    # O: -0.8476, H: +0.4238 (TIP3P-like for dipole calculation)
    charges = []
    for i, atom in enumerate(topology.atoms()):
        if atom.element.symbol == 'O':
            charges.append(-0.8476)
        else:
            charges.append(0.4238)
    
    # Create integrator
    print(f"\n--- Creating Integrator ---")
    friction = 1.0  # ps⁻¹
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt_ps * unit.picosecond
    )
    print(f"  LangevinMiddleIntegrator: friction = {friction} ps⁻¹, dt = {dt_ps} ps")
    
    # Create context
    print(f"\n--- Creating Simulation Context ---")
    platform = openmm.Platform.getPlatformByName('CUDA')
    print(f"  Using CUDA platform (GPU acceleration)")
    
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions * unit.nanometer)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    
    # Get initial state
    state = context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"  Initial PE: {pe:.1f} kJ/mol")
    
    # Storage for trajectory
    dipole_times = []
    dipole_values = []
    
    output_file = f"{output_prefix}_lambda0.0000.npz"
    
    print(f"\n--- Running Simulation ---")
    print(f"  Total steps: {total_steps:,}")
    print(f"  Equilibration: {equilibration_steps:,} steps ({equilibration_ps} ps)")
    print(f"  Production: {production_steps:,} steps ({production_ps} ps)")
    print(f"  Output file: {output_file}")
    print(f"  Progress updates every 100 steps")
    
    start_time = time.time()
    report_interval = 100  # Report every 100 steps for smaller systems
    
    for step in range(total_steps):
        integrator.step(1)
        
        # Record data every step during production
        if step >= equilibration_steps:
            state = context.getState(getPositions=True)
            pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            dipole = compute_dipole(pos, charges)
            
            dipole_times.append((step - equilibration_steps) * dt_ps)
            dipole_values.append(dipole.copy())
        
        # Progress report
        if (step + 1) % report_interval == 0 or (step + 1) == total_steps:
            phase = "EQUIL" if step < equilibration_steps else "PROD"
            pct = 100 * (step + 1) / total_steps
            elapsed = time.time() - start_time
            rate = (step + 1) / max(elapsed, 1e-8)
            eta = (total_steps - step - 1) / max(rate, 1e-8)
            
            # Calculate simulation speed (ns/day)
            ns_per_day = (rate * dt_ps * 86400) / 1000
            
            state = context.getState(getEnergy=True, getPositions=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            dipole_now = compute_dipole(pos, charges)
            
            # Save NPZ file during production
            if phase == "PROD" and len(dipole_times) > 0:
                np.savez(output_file,
                        time_ps=np.array(dipole_times),
                        dipole_nm=np.array(dipole_values),
                        metadata={
                            'lambda_coupling': 0.0,  # No cavity
                            'dt_ps': dt_ps,
                            'num_molecules': num_molecules,
                            'box_size_nm': box_size_nm,
                            'temperature_K': temperature_K,
                            'status': 'running',
                            'step': step + 1,
                            'total_steps': total_steps,
                            'model': 'uma-s-1p1',
                            'flexible_water': True
                        })
            
            print(f"[{pct:5.1f}%] {phase} | PE: {pe:10.1f} kJ/mol | "
                  f"d_xy: ({dipole_now[0]:6.2f}, {dipole_now[1]:6.2f}) e·nm | "
                  f"Speed: {ns_per_day:5.3f} ns/day | ETA: {eta/60:5.1f}m", flush=True)
    
    # Final save
    elapsed = time.time() - start_time
    print(f"\n--- Simulation Complete ---")
    print(f"  Total time: {elapsed/60:.1f} min ({elapsed/3600:.2f} hr)")
    print(f"  Performance: {(step+1)/elapsed:.1f} steps/s")
    
    np.savez(output_file,
             time_ps=np.array(dipole_times),
             dipole_nm=np.array(dipole_values),
             metadata={
                 'lambda_coupling': 0.0,
                 'dt_ps': dt_ps,
                 'num_molecules': num_molecules,
                 'box_size_nm': box_size_nm,
                 'temperature_K': temperature_K,
                 'status': 'complete',
                 'step': total_steps,
                 'total_steps': total_steps,
                 'elapsed_time_s': elapsed,
                 'model': 'uma-s-1p1',
                 'flexible_water': True
             })
    
    print(f"  Saved: {output_file}")
    print(f"  Data points: {len(dipole_times):,}")
    print(f"  Production time: {dipole_times[-1]:.1f} ps")
    print(f"\nNext step: Analyze IR spectrum")
    print(f"  python analyze_spectrum.py {output_file} --xlim 0 4500")


def main():
    """Main entry point with command-line argument parsing."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Baseline Water MD with UMA ML Potential for IR Validation'
    )
    parser.add_argument('--molecules', type=int, default=10,
                       help='Number of water molecules (default: 10)')
    parser.add_argument('--temp', type=float, default=300.0,
                       help='Temperature in K (default: 300)')
    parser.add_argument('--dt', type=float, default=0.0005,
                       help='Timestep in ps (default: 0.0005 = 0.5 fs)')
    parser.add_argument('--equil', type=float, default=10.0,
                       help='Equilibration time in ps (default: 10)')
    parser.add_argument('--prod', type=float, default=100.0,
                       help='Production time in ps (default: 100)')
    parser.add_argument('--output', type=str, default='uma_water_baseline',
                       help='Output prefix (default: uma_water_baseline)')
    parser.add_argument('--test', action='store_true',
                       help='Quick test mode (3 molecules, 1 ps)')
    
    args = parser.parse_args()
    
    # Override for test mode
    if args.test:
        print("\n*** QUICK TEST MODE ***\n")
        args.molecules = 3
        args.equil = 0.5
        args.prod = 1.0
    
    run_baseline_simulation(
        num_molecules=args.molecules,
        temperature_K=args.temp,
        dt_ps=args.dt,
        equilibration_ps=args.equil,
        production_ps=args.prod,
        output_prefix=args.output
    )


if __name__ == "__main__":
    main()
