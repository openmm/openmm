#!/usr/bin/env python3
"""
Baseline Water MD Simulation - Flexible Bonds for IR Validation
================================================================

System: Flexible TIP4P-Ew water model with intramolecular vibrations
Purpose: Validate IR spectrum shows OH stretch (~3400-3700 cm⁻¹) before cavity experiments
Protocol: NO cavity - pure water baseline for spectroscopic validation

This is a critical validation step. Rigid water models cannot produce OH stretch
vibrations because bonds are constrained. We must confirm flexible bonds work
and produce expected IR peaks before attempting cavity coupling experiments.

Expected IR Peaks (Experimental):
- OH stretch (H-bonded): 3200-3400 cm⁻¹ (broad, strong)
- Free OH stretch: ~3700 cm⁻¹ (shoulder)
- HOH bending: ~1645 cm⁻¹
- Librational modes: 400-900 cm⁻¹
"""

import sys
import numpy as np
import time
from pathlib import Path

try:
    from openmm import openmm
    from openmm import unit
    from openmm import app
    print("OpenMM loaded successfully")
except ImportError as e:
    print(f"Error importing OpenMM: {e}")
    sys.exit(1)


def create_flexible_water_box(num_molecules=1000, box_size_nm=3.0, temperature_K=300.0):
    """
    Create a water box using FLEXIBLE TIP4P-Ew model.
    
    Key difference from rigid models: O-H bonds can vibrate, enabling
    detection of OH stretch frequency in IR spectrum.
    
    Parameters
    ----------
    num_molecules : int
        Number of water molecules
    box_size_nm : float
        Box size in nanometers (cubic box)
    temperature_K : float
        Temperature in Kelvin
        
    Returns
    -------
    system : openmm.System
    topology : openmm.app.Topology
    positions : list of Vec3
    charges : list of float
    num_molecular : int
        Number of molecular atoms (excluding M-sites)
    """
    print(f"\n--- Creating Flexible Water Box ---")
    print(f"  Model: TIP4P-Ew (FLEXIBLE bonds for IR spectroscopy)")
    print(f"  Molecules: {num_molecules}")
    print(f"  Box size: {box_size_nm:.2f} nm")
    print(f"  Temperature: {temperature_K} K")
    
    # Generate PDB for water box
    pdb_lines = ["CRYST1{:9.3f}{:9.3f}{:9.3f}  90.00  90.00  90.00 P 1           1".format(
        box_size_nm * 10, box_size_nm * 10, box_size_nm * 10)]
    
    # TIP4P-Ew geometry
    oh_bond = 0.09572  # nm
    hoh_angle = 104.52 * np.pi / 180.0  # radians
    
    side = int(np.ceil(num_molecules**(1/3)))
    spacing = box_size_nm / side
    
    atom_idx = 1
    mol_count = 0
    
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if mol_count >= num_molecules:
                    break
                
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
                
                # Rotate and translate (convert to Angstrom for PDB)
                o_pos = (R @ o_local + np.array([cx, cy, cz])) * 10
                h1_pos = (R @ h1_local + np.array([cx, cy, cz])) * 10
                h2_pos = (R @ h2_local + np.array([cx, cy, cz])) * 10
                
                # Write to PDB
                resnum = mol_count + 1
                pdb_lines.append(f"HETATM{atom_idx:5d}  O   HOH  {resnum:4d}    {o_pos[0]:8.3f}{o_pos[1]:8.3f}{o_pos[2]:8.3f}  1.00  0.00           O")
                atom_idx += 1
                pdb_lines.append(f"HETATM{atom_idx:5d}  H1  HOH  {resnum:4d}    {h1_pos[0]:8.3f}{h1_pos[1]:8.3f}{h1_pos[2]:8.3f}  1.00  0.00           H")
                atom_idx += 1
                pdb_lines.append(f"HETATM{atom_idx:5d}  H2  HOH  {resnum:4d}    {h2_pos[0]:8.3f}{h2_pos[1]:8.3f}{h2_pos[2]:8.3f}  1.00  0.00           H")
                atom_idx += 1
                
                mol_count += 1
                
            if mol_count >= num_molecules:
                break
        if mol_count >= num_molecules:
            break
    
    pdb_lines.append("END")
    pdb_string = "\n".join(pdb_lines)
    
    # Write and load PDB
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(pdb_string)
        pdb_file = f.name
    
    pdb = app.PDBFile(pdb_file)
    Path(pdb_file).unlink()
    
    print(f"  Created {mol_count} water molecules ({len(pdb.positions)} atoms)")
    
    # Load flexible TIP4P-Ew force field
    print(f"\n  Loading FLEXIBLE TIP4P-Ew forcefield...")
    script_dir = Path(__file__).parent
    flexible_xml = script_dir / 'tip4pew_flexible.xml'
    
    if not flexible_xml.exists():
        print(f"  ERROR: {flexible_xml} not found!")
        print(f"  This file is required for flexible water bonds.")
        sys.exit(1)
    
    forcefield = app.ForceField(str(flexible_xml))
    print(f"  Loaded flexible TIP4P-Ew from {flexible_xml}")
    
    # Add M-sites (virtual particles)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addExtraParticles(forcefield)
    print(f"  Added virtual M-sites for TIP4P-Ew")
    
    # Create system with FLEXIBLE water
    print(f"\n  Creating system with flexible bonds...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=None,              # No constraints on O-H (allow vibration)
        rigidWater=False,
        ewaldErrorTolerance=0.0005
    )
    
    print(f"  System created: FLEXIBLE water (no bond constraints)")
    print(f"  O-H bonds can vibrate → OH stretch will appear in IR")
    
    # Extract charges for dipole calculation
    num_molecular = len(pdb.positions)
    charges = []
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for i in range(num_molecular):
                charge, sigma, epsilon = force.getParticleParameters(i)
                charges.append(charge.value_in_unit(unit.elementary_charge))
            break
    
    # Verify flexible bonds are present
    has_bonds = False
    has_angles = False
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            has_bonds = True
            print(f"  HarmonicBondForce: {force.getNumBonds()} bonds")
        if isinstance(force, openmm.HarmonicAngleForce):
            has_angles = True
            print(f"  HarmonicAngleForce: {force.getNumAngles()} angles")
    
    if not has_bonds or not has_angles:
        print(f"  WARNING: Missing harmonic terms - water may still be rigid!")
    
    print(f"  Total particles: {system.getNumParticles()} (includes {system.getNumParticles() - num_molecular} M-sites)")
    print(f"  Molecular atoms for dipole: {num_molecular}")
    
    return system, modeller.topology, modeller.getPositions(), charges, num_molecular


def setup_thermostat(system, temperature_K, num_molecular_particles, tau_bussi_ps=1.0):
    """
    Set up BussiThermostat for water molecules.
    
    For baseline simulation, all particles are water (no cavity to exclude).
    """
    print(f"\n--- Setting Up Thermostat ---")
    
    try:
        bussi = openmm.BussiThermostat(temperature_K, tau_bussi_ps)
        bussi.setApplyToAllParticles(False)
        
        # Add all molecular particles (O, H atoms - not M-sites)
        for i in range(num_molecular_particles):
            bussi.addParticle(i)
        
        system.addForce(bussi)
        print(f"  BussiThermostat added for {bussi.getNumParticles()} particles")
        print(f"  Temperature: {temperature_K} K, Tau: {tau_bussi_ps} ps")
        return True
    except AttributeError:
        print(f"  ⚠ BussiThermostat not available, will use Langevin only")
        return False


def compute_dipole(state, charges, num_molecular):
    """Compute total molecular dipole moment."""
    positions = state.getPositions(asNumpy=True)
    dipole = np.zeros(3)
    for i in range(num_molecular):
        pos = positions[i].value_in_unit(unit.nanometer)
        dipole += charges[i] * np.array(pos)
    return dipole


def run_baseline_simulation(num_molecules=1000, box_size_nm=3.0, temperature_K=300.0,
                           dt_ps=0.0005, equilibration_ps=100.0, production_ps=500.0,
                           output_prefix="water_baseline"):
    """
    Run baseline water MD with flexible bonds - NO CAVITY.
    
    This validates that:
    1. Flexible water bonds work (O-H can vibrate)
    2. IR spectrum shows OH stretch peak (~3400-3700 cm⁻¹)
    3. System is numerically stable with small timestep
    
    Parameters
    ----------
    num_molecules : int
        Number of water molecules
    box_size_nm : float
        Box size in nm
    temperature_K : float
        Temperature in K
    dt_ps : float
        Timestep in ps (MUST be small for flexible bonds, recommend 0.5 fs)
    equilibration_ps : float
        Equilibration time in ps
    production_ps : float
        Production time in ps (recommend ≥500 ps for good spectral resolution)
    output_prefix : str
        Prefix for output files
    """
    print("=" * 80)
    print("Baseline Water MD Simulation - Flexible Bonds")
    print("=" * 80)
    print("\nPurpose: Validate IR spectrum before cavity experiments")
    print("Expected: OH stretch peak at 3200-3700 cm⁻¹")
    
    print(f"\nSimulation Parameters:")
    print(f"  Water molecules: {num_molecules}")
    print(f"  Box size: {box_size_nm:.2f} nm")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt_ps} ps ({dt_ps*1000:.2f} fs)")
    print(f"  Equilibration: {equilibration_ps} ps")
    print(f"  Production: {production_ps} ps")
    
    if dt_ps > 0.0006:
        print(f"\n  ⚠ WARNING: Timestep {dt_ps} ps may be too large for flexible water!")
        print(f"  Recommend dt ≤ 0.5 fs (0.0005 ps) for numerical stability")
    
    # Calculate steps
    equilibration_steps = int(equilibration_ps / dt_ps)
    production_steps = int(production_ps / dt_ps)
    total_steps = equilibration_steps + production_steps
    
    # Create system
    np.random.seed(42)
    system, topology, positions, charges, num_molecular = create_flexible_water_box(
        num_molecules, box_size_nm, temperature_K
    )
    
    # Setup thermostat (NO CAVITY - thermostat all water)
    bussi_available = setup_thermostat(system, temperature_K, num_molecular, tau_bussi_ps=1.0)
    
    # Create integrator
    print(f"\n--- Creating Integrator ---")
    friction = 1.0  # ps⁻¹
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt_ps * unit.picosecond
    )
    print(f"  LangevinMiddleIntegrator: friction = {friction} ps⁻¹, dt = {dt_ps} ps")
    
    # Create simulation context
    print(f"\n--- Creating Simulation Context ---")
    platform = openmm.Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed'}
    print(f"  Using CUDA platform (GPU acceleration)")
    
    context = openmm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    
    # Minimize energy
    print(f"\n--- Minimizing Energy ---")
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=1000)
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
    print(f"  Streaming updates every 1 ps")
    
    start_time = time.time()
    
    for step in range(total_steps):
        integrator.step(1)
        
        # Record data every step during production
        if step >= equilibration_steps:
            state = context.getState(getPositions=True)
            dipole = compute_dipole(state, charges, num_molecular)
            
            dipole_times.append((step - equilibration_steps) * dt_ps)
            dipole_values.append(dipole.copy())
        
        # Progress report and save every 1 ps
        steps_per_ps = int(1.0 / dt_ps)
        if (step + 1) % steps_per_ps == 0 or (step + 1) == total_steps:
            phase = "EQUIL" if step < equilibration_steps else "PROD"
            pct = 100 * (step + 1) / total_steps
            elapsed = time.time() - start_time
            rate = (step + 1) / max(elapsed, 1e-8)
            eta = (total_steps - step - 1) / max(rate, 1e-8)
            
            # Calculate simulation speed (ns/day)
            ns_per_day = (rate * dt_ps * 86400) / 1000  # steps/s * ps/step * s/day / 1000
            
            state = context.getState(getEnergy=True, getPositions=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            dipole_now = compute_dipole(state, charges, num_molecular)
            
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
                            'flexible_water': True
                        })
            
            print(f"[{pct:5.1f}%] {phase} | PE: {pe:10.1f} kJ/mol | KE: {ke:10.1f} kJ/mol | "
                  f"d_xy: ({dipole_now[0]:6.2f}, {dipole_now[1]:6.2f}) e·nm | "
                  f"Speed: {ns_per_day:5.2f} ns/day | ETA: {eta/60:5.1f}m", flush=True)
    
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
                 'flexible_water': True
             })
    
    print(f"  Saved: {output_file}")
    print(f"  Data points: {len(dipole_times):,}")
    print(f"  Production time: {dipole_times[-1]:.1f} ps")
    print(f"\nNext step: Analyze IR spectrum to verify OH stretch peak")
    print(f"  python analyze_water_ir.py {output_file} --xlim 0 4500")


def main():
    """Main entry point with command-line argument parsing."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Baseline Water MD with Flexible Bonds for IR Validation'
    )
    parser.add_argument('--molecules', type=int, default=1000,
                       help='Number of water molecules (default: 1000)')
    parser.add_argument('--box', type=float, default=3.0,
                       help='Box size in nm (default: 3.0)')
    parser.add_argument('--temp', type=float, default=300.0,
                       help='Temperature in K (default: 300)')
    parser.add_argument('--dt', type=float, default=0.0005,
                       help='Timestep in ps (default: 0.0005 = 0.5 fs)')
    parser.add_argument('--equil', type=float, default=100.0,
                       help='Equilibration time in ps (default: 100)')
    parser.add_argument('--prod', type=float, default=500.0,
                       help='Production time in ps (default: 500)')
    parser.add_argument('--output', type=str, default='water_baseline',
                       help='Output prefix (default: water_baseline)')
    parser.add_argument('--test', action='store_true',
                       help='Quick test mode (216 molecules, 100 ps)')
    
    args = parser.parse_args()
    
    # Override for test mode
    if args.test:
        print("\nQuick test mode\n")
        args.molecules = 216
        args.box = 1.8
        args.equil = 50.0
        args.prod = 100.0
    
    run_baseline_simulation(
        num_molecules=args.molecules,
        box_size_nm=args.box,
        temperature_K=args.temp,
        dt_ps=args.dt,
        equilibration_ps=args.equil,
        production_ps=args.prod,
        output_prefix=args.output
    )


if __name__ == "__main__":
    main()
