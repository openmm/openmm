#!/usr/bin/env python3
"""
Cavity Molecular Dynamics Simulation of Water
==============================================

System: TIP4P-Ew water model (~1000 molecules)
Cavity: Resonant with OH stretch (3656 cm⁻¹)
Protocol: Equilibrate → Activate coupling → Observe nonthermal aging

This script demonstrates cavity-induced nonthermal aging in liquid water,
where strong light-matter coupling selectively pumps OH stretch vibrations
and drives the hydrogen bond network into deeper potential minima.
"""

import sys
import numpy as np
import time
from pathlib import Path

try:
    from openmm import openmm
    from openmm import unit
    from openmm import app
    from openmmforcefields.generators import SystemGenerator
    print("✓ OpenMM and openmmforcefields loaded successfully")
except ImportError as e:
    print(f"Error importing required packages: {e}")
    print("Make sure OpenMM and openmmforcefields are installed.")
    sys.exit(1)

# Unit conversions (matching paper conventions)
BOHR_TO_NM = 0.0529177           # 1 Bohr = 0.0529177 nm
HARTREE_TO_KJMOL = 2625.5        # 1 Hartree = 2625.5 kJ/mol
HARTREE_TO_CM = 219474.63        # 1 Hartree = 219474.63 cm⁻¹
AMU_TO_AU = 1822.888             # 1 amu = 1822.888 electron masses

def wavenumber_to_hartree(omega_cm):
    """Convert wavenumber (cm⁻¹) to Hartree (atomic units)."""
    return omega_cm / HARTREE_TO_CM


def create_water_box(num_molecules=1000, box_size_nm=3.0, temperature_K=300.0):
    """
    Create a water box using TIP4P-Ew model.
    
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
    print(f"\n--- Creating Water Box ---")
    print(f"  Model: TIP4P-Ew")
    print(f"  Molecules: {num_molecules}")
    print(f"  Box size: {box_size_nm:.2f} nm")
    print(f"  Temperature: {temperature_K} K")
    
    # Use the simpler approach: generate a PDB string and load it
    # This avoids the template matching issues
    pdb_lines = ["CRYST1{:9.3f}{:9.3f}{:9.3f}  90.00  90.00  90.00 P 1           1".format(
        box_size_nm * 10, box_size_nm * 10, box_size_nm * 10)]  # Convert nm to Angstrom
    
    # TIP4P-Ew geometry (standard water)
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
                
                # Center of mass position
                cx = (i + 0.5) * spacing
                cy = (j + 0.5) * spacing
                cz = (k + 0.5) * spacing
                
                # Random rotation
                np.random.seed(42 + mol_count)
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                psi = np.random.rand() * 2 * np.pi
                
                # Create water molecule in standard orientation then rotate
                h1_local = np.array([oh_bond * np.sin(hoh_angle/2), 0, oh_bond * np.cos(hoh_angle/2)])
                h2_local = np.array([-oh_bond * np.sin(hoh_angle/2), 0, oh_bond * np.cos(hoh_angle/2)])
                o_local = np.array([0, 0, 0])
                
                # Rotation matrix (Euler angles)
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
                
                # Add to PDB
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
    
    # Write temporary PDB file
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(pdb_string)
        pdb_file = f.name
    
    # Load PDB
    pdb = app.PDBFile(pdb_file)
    
    # Clean up temp file
    Path(pdb_file).unlink()
    
    print(f"  Created {mol_count} water molecules ({len(pdb.positions)} atoms)")
    
    # Create forcefield and modeller
    print(f"\n  Generating system with FLEXIBLE TIP4P-Ew forcefield...")
    # Use flexible TIP4P-Ew XML (same as baseline) to allow O-H vibrations for IR
    script_dir = Path(__file__).parent
    flexible_xml = script_dir / 'tip4pew_flexible.xml'
    
    if not flexible_xml.exists():
        print(f"  ERROR: {flexible_xml} not found!")
        print(f"  This file is required for flexible water bonds.")
        sys.exit(1)
    
    forcefield = app.ForceField(str(flexible_xml))
    print(f"  ✓ Loaded flexible TIP4P-Ew from {flexible_xml}")
    
    # Use Modeller to add M-sites
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addExtraParticles(forcefield)
    
    # Create system with FLEXIBLE water
    print(f"\n  Creating system with flexible bonds...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=None,              # CRITICAL: No constraints on O-H bonds
        rigidWater=False,              # CRITICAL: Allow water to vibrate
        ewaldErrorTolerance=0.0005
    )
    
    print(f"  ✓ System created: FLEXIBLE water (no bond constraints)")
    print(f"  ✓ O-H bonds can vibrate → OH stretch will appear in IR")
    
    # Extract charges for dipole calculation (only O and H atoms, not M-sites)
    num_molecular = len(pdb.positions)  # Original atoms before adding M-sites
    charges = []
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for i in range(num_molecular):
                charge, sigma, epsilon = force.getParticleParameters(i)
                charges.append(charge.value_in_unit(unit.elementary_charge))
            break
    
    print(f"  System created with {system.getNumParticles()} particles (includes M-sites)")
    print(f"  Molecular atoms for dipole: {num_molecular}")
    print(f"  Forces: {[type(f).__name__ for f in system.getForces()]}")
    
    return system, modeller.topology, modeller.getPositions(), charges, num_molecular


def add_cavity_particle(system, positions, omegac_au, photon_mass_amu):
    """
    Add cavity photon mode as a dummy particle.
    
    Parameters
    ----------
    system : openmm.System
    positions : list
        List of positions to append to
    omegac_au : float
        Cavity frequency in atomic units (Hartree)
    photon_mass_amu : float
        Photon mass in amu
        
    Returns
    -------
    cavity_index : int
        Index of the cavity particle
    """
    print(f"\n--- Adding Cavity Particle ---")
    
    # Add cavity particle at box center
    cavity_index = system.addParticle(photon_mass_amu * unit.amu)
    positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)
    
    # Add to nonbonded force with no interactions
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.addParticle(0.0, 0.1 * unit.nanometer, 0.0)
            break
    
    omegac_cm = omegac_au * HARTREE_TO_CM
    print(f"  Cavity particle index: {cavity_index}")
    print(f"  Omega_c: {omegac_au:.6f} a.u. = {omegac_cm:.1f} cm⁻¹")
    print(f"  Photon mass: {photon_mass_amu:.6f} amu")
    
    return cavity_index


def setup_cavity_coupling(system, cavity_index, omegac_au, lambda_coupling, 
                         photon_mass_amu):
    """
    Set up cavity coupling.
    
    Parameters
    ----------
    system : openmm.System
    cavity_index : int
    omegac_au : float
        Cavity frequency in Hartree
    lambda_coupling : float
        Coupling strength (dimensionless)
    photon_mass_amu : float
        Photon mass in amu
        
    Returns
    -------
    cavity_force : openmm.CavityForce
    displacer : openmm.CavityParticleDisplacer
    """
    print(f"\n--- Setting Up Cavity Coupling ---")
    print(f"  Lambda: {lambda_coupling}")
    print(f"  Coupling will be ON from t=0 (no gradual activation)")
    
    # CavityForce - set coupling directly to target value
    cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass_amu)
    system.addForce(cavity_force)
    print(f"  CavityForce added with lambda={lambda_coupling}")
    
    # CavityParticleDisplacer for finite-q correction - also set to target lambda
    displacer = openmm.CavityParticleDisplacer(cavity_index, omegac_au, photon_mass_amu)
    displacer.setSwitchOnStep(0)
    displacer.setSwitchOnLambda(lambda_coupling)  # Set to actual coupling strength!
    system.addForce(displacer)
    print(f"  CavityParticleDisplacer added with lambda={lambda_coupling}")
    
    return cavity_force, displacer


def setup_thermostats(system, temperature_K, num_molecular_particles, tau_bussi_ps=1.0):
    """
    Set up BussiThermostat for water molecules.
    
    Parameters
    ----------
    system : openmm.System
    temperature_K : float
    num_molecular_particles : int
        Number of molecular (non-cavity) particles
    tau_bussi_ps : float
        Bussi thermostat time constant in ps
    """
    print(f"\n--- Setting Up Thermostats ---")
    
    try:
        bussi = openmm.BussiThermostat(temperature_K, tau_bussi_ps)
        bussi.setApplyToAllParticles(False)
        
        # Add all molecular particles (exclude cavity)
        for i in range(num_molecular_particles):
            bussi.addParticle(i)
        
        system.addForce(bussi)
        print(f"  BussiThermostat added for {bussi.getNumParticles()} particles")
        print(f"  Temperature: {temperature_K} K, Tau: {tau_bussi_ps} ps")
        return True
    except AttributeError:
        print(f"  BussiThermostat not available, will use Langevin only")
        return False


def compute_dipole(state, charges, num_molecular):
    """Compute total molecular dipole moment (excluding cavity particle)."""
    positions = state.getPositions(asNumpy=True)
    dipole = np.zeros(3)
    for i in range(num_molecular):
        pos = positions[i].value_in_unit(unit.nanometer)
        dipole += charges[i] * np.array(pos)
    return dipole


def run_simulation(lambda_coupling=0.01, num_molecules=1000, box_size_nm=3.0,
                  temperature_K=300.0, dt_ps=0.001, equilibration_ps=100.0,
                  production_ps=900.0, output_prefix="water_cavity"):
    """
    Run the water cavity MD simulation with streaming output.
    
    Parameters
    ----------
    lambda_coupling : float
        Coupling strength
    num_molecules : int
        Number of water molecules
    box_size_nm : float
        Box size in nm
    temperature_K : float
        Temperature in K
    dt_ps : float
        Timestep in ps
    equilibration_ps : float
        Equilibration time in ps
    production_ps : float
        Production time in ps
    output_prefix : str
        Prefix for output files
    """
    print("=" * 80)
    print("Water Cavity Molecular Dynamics Simulation")
    print("=" * 80)
    
    # Cavity parameters (OH stretch resonance - MESA-determined)
    cavity_freq_cm = 1600  # cm⁻¹ (OH stretch peak from MESA analysis of baseline)
    omegac_au = wavenumber_to_hartree(cavity_freq_cm)
    photon_mass_amu = 1.0 / AMU_TO_AU  # 1 a.u. electron mass -> amu
    
    print(f"\nSimulation Parameters:")
    print(f"  Water molecules: {num_molecules}")
    print(f"  Box size: {box_size_nm:.2f} nm")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Cavity frequency: {cavity_freq_cm} cm⁻¹ (OH stretch - MESA analysis)")
    print(f"  Lambda coupling: {lambda_coupling}")
    print(f"  Timestep: {dt_ps} ps")
    print(f"  Equilibration: {equilibration_ps} ps")
    print(f"  Production: {production_ps} ps")
    
    # Calculate steps
    equilibration_steps = int(equilibration_ps / dt_ps)
    production_steps = int(production_ps / dt_ps)
    total_steps = equilibration_steps + production_steps
    
    # Create system
    np.random.seed(42)
    system, topology, positions, charges, num_molecular = create_water_box(
        num_molecules, box_size_nm, temperature_K
    )
    
    # Add cavity particle
    cavity_index = add_cavity_particle(system, positions, omegac_au, photon_mass_amu)
    
    # Setup cavity coupling
    cavity_force, displacer = setup_cavity_coupling(
        system, cavity_index, omegac_au, lambda_coupling, photon_mass_amu
    )
    
    # Setup thermostats
    bussi_available = setup_thermostats(system, temperature_K, num_molecular, tau_bussi_ps=1.0)
    
    # Create integrator
    print(f"\n--- Creating Integrator ---")
    friction = 0.5  # ps⁻¹ for cavity mode (low to preserve OH stretch vibrations)
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt_ps * unit.picosecond
    )
    print(f"  LangevinMiddleIntegrator: friction = {friction} ps⁻¹, dt = {dt_ps} ps")
    print(f"  Note: Low friction preserves vibrational dynamics for IR spectroscopy")
    
    # Create simulation context
    print(f"\n--- Creating Simulation Context ---")
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed'}
        print(f"  Using CUDA platform (GPU acceleration)")
    except Exception:
        try:
            platform = openmm.Platform.getPlatformByName('OpenCL')
            properties = {}
            print(f"  Using OpenCL platform")
        except Exception:
            platform = openmm.Platform.getPlatformByName('Reference')
            properties = {}
            print(f"  Using Reference platform (CPU, will be slow)")
    
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
    cavity_values = []
    
    output_file = f"{output_prefix}_lambda{lambda_coupling:.4f}.npz"
    
    print(f"\n--- Running Simulation ---")
    print(f"  Total steps: {total_steps}")
    print(f"  Equilibration: {equilibration_steps} steps ({equilibration_ps} ps)")
    print(f"  Production: {production_steps} steps ({production_ps} ps)")
    print(f"  Output file: {output_file}")
    print(f"  Streaming updates every 1 ps (1000 steps)")
    print(f"  Note: Simplified version - coupling constant throughout")
    
    start_time = time.time()
    last_save_time = start_time
    
    for step in range(total_steps):
        integrator.step(1)
        
        # Record data every step during production
        if step >= equilibration_steps:
            state = context.getState(getPositions=True)
            pos = state.getPositions(asNumpy=True)
            dipole = compute_dipole(state, charges, num_molecular)
            q_cav = pos[cavity_index].value_in_unit(unit.nanometer)
            
            dipole_times.append((step - equilibration_steps) * dt_ps)
            dipole_values.append(dipole.copy())
            cavity_values.append(q_cav.copy())
        
        # Progress report and save every 1000 steps (1 ps)
        if (step + 1) % 1000 == 0 or (step + 1) == total_steps:
            phase = "EQUIL" if step < equilibration_steps else "PROD"
            pct = 100 * (step + 1) / total_steps
            elapsed = time.time() - start_time
            rate = (step + 1) / max(elapsed, 1e-8)
            eta = (total_steps - step - 1) / max(rate, 1e-8)
            
            state = context.getState(getEnergy=True, getPositions=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            pos = state.getPositions(asNumpy=True)
            q_cav = pos[cavity_index].value_in_unit(unit.nanometer)
            dipole_now = compute_dipole(state, charges, num_molecular)
            
            # Save NPZ file during production
            if phase == "PROD" and len(dipole_times) > 0:
                np.savez(output_file,
                        time_ps=np.array(dipole_times),
                        dipole_nm=np.array(dipole_values),
                        cavity_nm=np.array(cavity_values),
                        metadata={
                            'lambda_coupling': lambda_coupling,
                            'cavity_freq_cm': cavity_freq_cm,
                            'omega_c_au': omegac_au,
                            'dt_ps': dt_ps,
                            'num_molecules': num_molecules,
                            'box_size_nm': box_size_nm,
                            'temperature_K': temperature_K,
                            'photon_mass_amu': photon_mass_amu,
                            'flexible_water': True,  # Track water model type
                            'status': 'running',
                            'step': step + 1,
                            'total_steps': total_steps
                        })
            
            print(f"[{pct:5.1f}%] {phase} | PE: {pe:10.1f} kJ/mol | "
                  f"d_xy: ({dipole_now[0]:6.2f}, {dipole_now[1]:6.2f}) e·nm | "
                  f"q_cav: ({q_cav[0]:6.3f}, {q_cav[1]:6.3f}) nm | "
                  f"ETA: {eta/60:5.1f}m", flush=True)
    
    # Final save
    elapsed = time.time() - start_time
    print(f"\n--- Simulation Complete ---")
    print(f"  Total time: {elapsed/60:.1f} min ({elapsed/3600:.2f} hr)")
    print(f"  Performance: {(step+1)/elapsed:.1f} steps/s")
    
    np.savez(output_file,
             time_ps=np.array(dipole_times),
             dipole_nm=np.array(dipole_values),
             cavity_nm=np.array(cavity_values),
             metadata={
                 'lambda_coupling': lambda_coupling,
                 'cavity_freq_cm': cavity_freq_cm,
                 'omega_c_au': omegac_au,
                 'dt_ps': dt_ps,
                 'num_molecules': num_molecules,
                 'box_size_nm': box_size_nm,
                 'temperature_K': temperature_K,
                 'photon_mass_amu': photon_mass_amu,
                 'flexible_water': True,  # Track water model type
                 'status': 'complete',
                 'step': total_steps,
                 'total_steps': total_steps,
                 'elapsed_time_s': elapsed
             })
    
    print(f"  Saved: {output_file}")
    print(f"  Data points: {len(dipole_times)}")
    print(f"  Production time: {dipole_times[-1]:.1f} ps")


def main():
    """Main entry point with command-line argument parsing."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Water Cavity MD Simulation')
    parser.add_argument('--lambda', type=float, default=0.01, dest='lambda_coupling',
                       help='Coupling strength (default: 0.01)')
    parser.add_argument('--molecules', type=int, default=1000,
                       help='Number of water molecules (default: 1000)')
    parser.add_argument('--box', type=float, default=3.0,
                       help='Box size in nm (default: 3.0)')
    parser.add_argument('--temp', type=float, default=300.0,
                       help='Temperature in K (default: 300)')
    parser.add_argument('--dt', type=float, default=0.0005,
                       help='Timestep in ps (default: 0.0005 = 0.5 fs, required for flexible water)')
    parser.add_argument('--equil', type=float, default=100.0,
                       help='Equilibration time in ps (default: 100)')
    parser.add_argument('--prod', type=float, default=900.0,
                       help='Production time in ps (default: 900)')
    parser.add_argument('--output', type=str, default='water_cavity',
                       help='Output prefix (default: water_cavity)')
    parser.add_argument('--test', action='store_true',
                       help='Quick test run (216 molecules, 10 ps)')
    
    args = parser.parse_args()
    
    # Override for test mode
    if args.test:
        print("\n*** QUICK TEST MODE ***\n")
        args.molecules = 216
        args.box = 1.8
        args.equil = 5.0
        args.prod = 5.0
    
    run_simulation(
        lambda_coupling=args.lambda_coupling,
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
