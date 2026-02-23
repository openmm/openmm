#!/usr/bin/env python3
"""
Molecular Dynamics Simulation of Water with CACE-LR
==================================================

System: Water molecules using CACE-LR potential (with LES charges)
Goal: Run MD and record dipole moments for IR spectra calculation.
"""

import sys
import numpy as np
import time
import os
from pathlib import Path
from ase import Atoms
from ase.build import molecule

# Add current directory to path for imports
sys.path.append(str(Path(__file__).parent))

try:
    import openmm
    from openmm import unit
    from openmm.app import Topology, Element, Simulation, StateDataReporter, PDBFile
    from openmmml import MLPotential
    import torch
    from cace_dipole_calculator import CACEDipoleCalculator
    print("OpenMM, OpenMM-ML, and CACE loaded successfully")
except ImportError as e:
    print(f"Error importing required packages: {e}")
    sys.exit(1)

# Constants
HARTREE_TO_KJMOL = 2625.5
AMU_TO_AU = 1822.888
DEBYE_TO_E_A = 1.0 / 4.803

def create_water_topology(num_molecules=32):
    """Create a water box topology."""
    print(f"\n--- Creating Water Topology ---")
    print(f"  Molecules: {num_molecules}")
    
    # Simple TIP3P-like geometry for initialization
    oh_bond = 0.09572  # nm
    hoh_angle = 104.52 * np.pi / 180.0
    
    # Calculate box size for liquid water density (~1.0 g/cm³ at 300K)
    # Molar volume of water at 300K: ~18.07 cm³/mol = 0.0301 nm³/molecule
    vol_per_mol_nm3 = 0.0301  # More accurate for 300K
    total_vol = num_molecules * vol_per_mol_nm3
    box_size_nm = (total_vol ** (1/3))  # Start at target density, barostat will adjust
    
    print(f"  Box size: {box_size_nm:.2f} nm")
    
    topology = Topology()
    positions = []
    
    side = int(np.ceil(num_molecules**(1/3)))
    spacing = box_size_nm / side
    
    mol_count = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if mol_count >= num_molecules:
                    break
                    
                chain = topology.addChain()
                residue = topology.addResidue('HOH', chain)
                topology.addAtom('O', Element.getBySymbol('O'), residue)
                topology.addAtom('H', Element.getBySymbol('H'), residue)
                topology.addAtom('H', Element.getBySymbol('H'), residue)
                
                cx = (i + 0.5) * spacing
                cy = (j + 0.5) * spacing
                cz = (k + 0.5) * spacing
                
                # Random rotation
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                psi = np.random.rand() * 2 * np.pi
                
                h1_local = np.array([oh_bond * np.sin(hoh_angle/2), 0, oh_bond * np.cos(hoh_angle/2)])
                h2_local = np.array([-oh_bond * np.sin(hoh_angle/2), 0, oh_bond * np.cos(hoh_angle/2)])
                o_local = np.array([0, 0, 0])
                
                # Rotation matrix
                Rz = np.array([[np.cos(theta), -np.sin(theta), 0],
                               [np.sin(theta), np.cos(theta), 0],
                               [0, 0, 1]])
                Ry = np.array([[np.cos(phi), 0, np.sin(phi)],
                               [0, 1, 0],
                               [-np.sin(phi), 0, np.cos(phi)]])
                Rx = np.array([[1, 0, 0],
                               [0, np.cos(psi), -np.sin(psi)],
                               [0, np.sin(psi), np.cos(psi)]])
                R = Rz @ Ry @ Rx
                
                positions.append((R @ o_local + np.array([cx, cy, cz])).tolist())
                positions.append((R @ h1_local + np.array([cx, cy, cz])).tolist())
                positions.append((R @ h2_local + np.array([cx, cy, cz])).tolist())
                
                mol_count += 1
    
    return topology, positions, box_size_nm

def run_simulation(num_molecules=32, temperature_K=300.0, 
                  dt_ps=0.004, equilibration_ps=10.0, production_ps=50.0,
                  output_prefix="water_cace_lr"):
    """
    Run CACE-LR water MD simulation.
    """
    print("=" * 80)
    print("CACE-LR Water MD Simulation")
    print("=" * 80)
    
    # Model path (openmm root is 2 levels up from tests/cace-lr_water/)
    openmm_root = Path(__file__).resolve().parents[2]
    model_path = str(openmm_root / "cace" / "water" / "fit" / "fit_version_1" / "best_model.pth")
    if not os.path.exists(model_path):
        print(f"Error: Model not found at {model_path}")
        sys.exit(1)
        
    # Setup topology
    np.random.seed(42)
    topology, positions, box_size_nm = create_water_topology(num_molecules)
    n_atoms = len(positions)
    
    # Set periodic box vectors on topology BEFORE creating system
    # PBC detected from box
    # Topology.setPeriodicBoxVectors() takes a single tuple of 3 vectors
    topology.setPeriodicBoxVectors((
        [box_size_nm, 0, 0] * unit.nanometer,
        [0, box_size_nm, 0] * unit.nanometer,
        [0, 0, box_size_nm] * unit.nanometer
    ))
    
    # Create CACE potential
    print(f"\n--- Creating CACE-LR Potential ---")
    potential = MLPotential('cace-lr', model_path=model_path)
    
    # Create system
    system = potential.createSystem(topology)
    # Also set on system (redundant check)
    system.setDefaultPeriodicBoxVectors(
        [box_size_nm, 0, 0] * unit.nanometer,
        [0, box_size_nm, 0] * unit.nanometer,
        [0, 0, box_size_nm] * unit.nanometer
    )
    
    # Add NPT barostat for constant pressure simulation
    pressure_bar = 1.0  # 1 bar = 1 atm
    barostat = openmm.MonteCarloBarostat(
        pressure_bar * unit.bar,
        temperature_K * unit.kelvin
    )
    barostat.setFrequency(25)  # Attempt volume change every 25 steps
    system.addForce(barostat)
    print(f"  Added MonteCarloBarostat: {pressure_bar} bar, {temperature_K} K")
    
    # Verify PBC is enabled
    uses_pbc = system.usesPeriodicBoundaryConditions()
    print(f"  System uses PBC: {uses_pbc}")
    
    # Create dipole calculator for recording
    # We use CPU for recording to not interfere with GPU simulation if any
    dipole_calc = CACEDipoleCalculator(model_path, device='cpu')
    
    # Create integrator
    friction = 1.0  # ps⁻¹
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt_ps * unit.picosecond
    )
    
    # Create context
    platform = openmm.Platform.getPlatformByName('CUDA') if torch.cuda.is_available() else openmm.Platform.getPlatformByName('Reference')
    print(f"  Using {platform.getName()} platform")
    
    context = openmm.Context(system, integrator, platform)
    
    # Set positions
    pos_with_units = [openmm.Vec3(p[0], p[1], p[2]) * unit.nanometer for p in positions]
    context.setPositions(pos_with_units)
    
    # Energy minimization to remove bad contacts
    print(f"\n--- Running Energy Minimization ---", flush=True)
    try:
        initial_state = context.getState(getEnergy=True)
        initial_pe = initial_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        print(f"  Initial PE: {initial_pe:.2f} kJ/mol", flush=True)
        
        print(f"  Running minimizer (max 1000 iterations)...", flush=True)
        openmm.LocalEnergyMinimizer.minimize(context, tolerance=10.0, maxIterations=1000)
        
        minimized_state = context.getState(getEnergy=True)
        minimized_pe = minimized_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        print(f"  Minimized PE: {minimized_pe:.2f} kJ/mol", flush=True)
        print(f"  Energy change: {minimized_pe - initial_pe:.2f} kJ/mol", flush=True)
    except Exception as e:
        print(f"  Energy minimization failed: {e}", flush=True)
        print(f"  Continuing without minimization...", flush=True)
        import traceback
        traceback.print_exc()
    
    # Now set velocities
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    
    # Equilibration
    equil_steps = int(equilibration_ps / dt_ps)
    print(f"\n--- Running Equilibration ({equilibration_ps} ps, {equil_steps} steps) ---")
    
    equil_report = max(1, equil_steps // 10)
    for i in range(equil_steps):
        integrator.step(1)
        if (i + 1) % equil_report == 0:
            state = context.getState(getEnergy=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            print(f"  Equil Step {i+1}/{equil_steps} | PE: {pe:.2f} kJ/mol", flush=True)
    
    # Production
    prod_steps = int(production_ps / dt_ps)
    output_file = f"{output_prefix}.npz"
    pdb_file = f"{output_prefix}.pdb"
    print(f"\n--- Running Production ({production_ps} ps, {prod_steps} steps) ---")
    print(f"  Saving to: {output_file}")
    print(f"  PDB trajectory: {pdb_file}")
    print(f"  NPT simulation: {pressure_bar} bar, {temperature_K} K")
    print(f"  Timestep: {dt_ps*1000:.1f} fs | Sampling rate: {1.0/dt_ps:.0f} THz")
    if dt_ps > 0.001:
        print(f"  ⚠ Warning: Timestep > 1 fs may be too coarse for IR spectroscopy.")
        print(f"    Recommended: 0.5-1.0 fs (0.0005-0.001 ps) for high-frequency vibrations")
    
    # Helper function to write PDB frame with box vectors from state
    def write_pdb_frame(topology, state, file_handle, model_num=None):
        """Write a PDB frame with box vectors from state."""
        from openmm.app.internal.unitcell import computeLengthsAndAngles
        from datetime import date
        from openmm import Platform
        import math
        
        if model_num is not None:
            file_handle.write(f"MODEL     {model_num:4d}\n")
        
        # Write header
        print("REMARK   1 CREATED WITH OPENMM %s, %s" % (Platform.getOpenMMVersion(), str(date.today())), file=file_handle)
        
        # Get box vectors from state and write CRYST1 record
        box_vectors = state.getPeriodicBoxVectors()
        if box_vectors is not None:
            # computeLengthsAndAngles returns unitless values: lengths in nm, angles in radians
            a_nm, b_nm, c_nm, alpha_rad, beta_rad, gamma_rad = computeLengthsAndAngles(box_vectors)
            
            # Convert angles to degrees and lengths to Angstroms for PDB format
            alpha_deg = alpha_rad * 180.0 / math.pi
            beta_deg = beta_rad * 180.0 / math.pi
            gamma_deg = gamma_rad * 180.0 / math.pi
            
            # Write CRYST1 record (lengths in Angstroms)
            print("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 " % (
                a_nm*10, b_nm*10, c_nm*10, alpha_deg, beta_deg, gamma_deg), file=file_handle)
        
        # Write positions
        positions = state.getPositions(asNumpy=True)
        PDBFile.writeModel(topology, positions, file=file_handle)
        
        if model_num is not None:
            file_handle.write("ENDMDL\n")
    
    # Open PDB file for writing trajectory
    pdb_handle = open(pdb_file, 'w')
    # Write initial frame (with periodic box vectors)
    state_init = context.getState(getPositions=True, getVelocities=False, enforcePeriodicBox=True)
    write_pdb_frame(topology, state_init, pdb_handle, model_num=1)
    frame_count = 2
    
    dipole_times = []
    dipole_values = []
    
    # Open dipole data file for incremental writing
    dipole_data_file = f"{output_prefix}_dipoles_realtime.txt"
    dipole_file_handle = open(dipole_data_file, 'w')
    dipole_file_handle.write("# Real-time dipole moment data\n")
    dipole_file_handle.write("# time_ps dipole_x dipole_y dipole_z\n")
    dipole_file_handle.flush()
    print(f"  Real-time dipole output: {dipole_data_file}")
    
    # Save checkpoint interval (write to disk every N steps)
    checkpoint_interval = 50  # Save every 50 steps for real-time analysis
    
    # Get chemical symbols once
    symbols = [atom.element.symbol for atom in topology.atoms()]
    
    start_time = time.time()
    report_interval = max(1, prod_steps // 20)
    pdb_interval = max(1, int(0.1 / dt_ps))  # Write PDB every ~0.1 ps
    
    for step in range(prod_steps):
        integrator.step(1)
        
        # Record dipole every step
        # For IR spectroscopy, timestep should be 0.5-1.0 fs (0.0005-0.001 ps)
        # to properly sample high-frequency vibrations (e.g., O-H stretch ~3000-3500 cm⁻¹)
        # 
        # Unwrapped positions for dipole (PBC continuity)
        # Do NOT use enforcePeriodicBox=True, otherwise molecules split across
        # boundaries will give artificial large dipoles
        state = context.getState(getPositions=True)
        pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        
        # Convert to ASE for dipole calc
        # Note: CACEDipoleCalculator handles conversion internally
        atoms = Atoms(symbols=symbols, positions=pos_nm * 10.0, 
                     cell=[box_size_nm*10, box_size_nm*10, box_size_nm*10], pbc=True)
        
        dipole, charges = dipole_calc.compute_dipole(atoms)
        
        current_time = step * dt_ps
        dipole_times.append(current_time)
        dipole_values.append(dipole.copy())
        
        # Write dipole to file immediately for real-time access
        dipole_file_handle.write(f"{current_time:.6f} {dipole[0]:.6f} {dipole[1]:.6f} {dipole[2]:.6f}\n")
        
        # Flush to disk every checkpoint_interval steps
        if (step + 1) % checkpoint_interval == 0:
            dipole_file_handle.flush()
            # Also save a checkpoint NPZ file
            if (step + 1) % (checkpoint_interval * 10) == 0:  # Every 500 steps
                checkpoint_file = f"{output_prefix}_checkpoint.npz"
                np.savez(checkpoint_file,
                        time_ps=np.array(dipole_times),
                        dipole_eA=np.array(dipole_values),
                        metadata={'dt_ps': dt_ps, 'num_molecules': num_molecules, 
                                 'steps_completed': step + 1, 'total_steps': prod_steps})
        
        # Write PDB frame periodically
        if (step + 1) % pdb_interval == 0:
            write_pdb_frame(topology, state, pdb_handle, model_num=frame_count)
            frame_count += 1
        
        if (step + 1) % report_interval == 0:
            elapsed = time.time() - start_time
            rate = (step + 1) / max(elapsed, 1e-8)
            ns_per_day = (rate * dt_ps * 86400) / 1000
            # Get current box size for density tracking
            state_box = state.getPeriodicBoxVectors(asNumpy=True)
            box_a = state_box[0].value_in_unit(unit.nanometer)
            current_box_nm = np.linalg.norm(box_a)
            vol_nm3 = current_box_nm ** 3
            density_gcm3 = (num_molecules * 18.01528 / 6.022e23) / (vol_nm3 * 1e-21)
            print(f"  Step {step+1}/{prod_steps} | Dipole: {np.linalg.norm(dipole):.3f} e·A | Box: {current_box_nm:.3f} nm | Density: {density_gcm3:.3f} g/cm³ | Speed: {ns_per_day:.2f} ns/day", flush=True)

    # Close dipole file
    dipole_file_handle.close()
    
    # Final save
    save_data = {
        'time_ps': np.array(dipole_times),
        'dipole_eA': np.array(dipole_values),
        'metadata': {
            'dt_ps': dt_ps,
            'num_molecules': num_molecules,
            'box_size_nm': box_size_nm,
            'temperature_K': temperature_K,
            'model': 'cace-lr-water'
        }
    }
    # Write final PDB frame
    state_final = context.getState(getPositions=True, getVelocities=False, enforcePeriodicBox=True)
    write_pdb_frame(topology, state_final, pdb_handle, model_num=frame_count)
    pdb_handle.close()
    
    np.savez(output_file, **save_data)
    print(f"\nSimulation Complete. Saved {len(dipole_times)} data points.")
    print(f"PDB trajectory saved: {pdb_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='CACE-LR Water MD')
    parser.add_argument('--molecules', type=int, default=32)
    parser.add_argument('--equil', type=float, default=2.0)
    parser.add_argument('--prod', type=float, default=10.0)
    parser.add_argument('--dt', type=float, default=0.001, 
                       help='Timestep in ps (default: 0.001 ps = 1 fs, recommended for IR: 0.0005-0.001 ps)')
    parser.add_argument('--output', type=str, default='water_cace_lr')
    parser.add_argument('--test', action='store_true')
    
    args = parser.parse_args()
    
    if args.test:
        args.molecules = 2
        args.equil = 0.1
        args.prod = 0.5
        
    run_simulation(
        num_molecules=args.molecules,
        equilibration_ps=args.equil,
        production_ps=args.prod,
        dt_ps=args.dt
    )
