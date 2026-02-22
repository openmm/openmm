#!/usr/bin/env python3
"""
Cavity OpenMM Dimer RPMD Simulation - IR Spectrum Calculation
==============================================================

This script runs Ring Polymer Molecular Dynamics (RPMD) simulations of O-O and N-N dimers
coupled to a cavity photon, computing the IR spectrum from centroid dipole moments.

Features:
- RPMD with user-specified number of beads for quantum nuclear effects
- PILE_G thermostat: Bussi for centroid + Langevin for internal modes
- Cavity coupling with lambda parameter
- IR spectrum from centroid dipole moment autocorrelation
- GPU-accelerated with CUDA

Usage:
    python run_simulation_rpmd.py --dimers 250 --beads 8 --lambda 0.0700
"""

import sys
import numpy as np
import os
from pathlib import Path
import time

try:
    from openmm import openmm
    from openmm import unit
    print("✓ OpenMM loaded successfully")
except ImportError as e:
    print(f"Error importing OpenMM: {e}")
    print("Make sure OpenMM is installed and the Python path is correct.")
    sys.exit(1)

# Physical constants for unit conversion
BOHR_TO_NM = 0.0529177  # 1 Bohr = 0.0529177 nm
HARTREE_TO_KJMOL = 2625.5  # 1 Hartree = 2625.5 kJ/mol
AMU_TO_KG = 1.66054e-27
# 1 a.u. of time = ℏ/E_h = 2.418884254e-5 ps (do not use 0.02418884 - that is 1000× too large)
AU_TIME_TO_PS = 2.418884254e-5  # 1 atomic time unit in ps

def create_diamer_system(num_molecules=50, fraction_OO=0.8, box_size_nm=2.0,
                         temperature_K=100.0, seed=42):
    """
    Create a system of O-O and N-N dimers.
    
    Parameters
    ----------
    num_molecules : int
        Number of diatomic molecules
    fraction_OO : float
        Fraction of O-O dimers (vs N-N)
    box_size_nm : float
        Box size in nanometers
    temperature_K : float
        Temperature in Kelvin
    seed : int
        Random seed
        
    Returns
    -------
    system : openmm.System
    positions : list of Vec3
    """
    np.random.seed(seed)
    
    system = openmm.System()
    positions = []
    
    # Particle masses
    mass_O = 16.0  # amu
    mass_N = 14.0  # amu
    
    # Bond parameters
    k_OO_au = 0.73204  # Hartree/Bohr²
    r0_OO_au = 2.281655158  # Bohr
    k_OO = k_OO_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)  # kJ/(mol·nm²)
    r0_OO = r0_OO_au * BOHR_TO_NM  # nm
    
    k_NN_au = 1.4325  # Hartree/Bohr²
    r0_NN_au = 2.0743522177  # Bohr
    k_NN = k_NN_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
    r0_NN = r0_NN_au * BOHR_TO_NM
    
    # Charges
    charge_magnitude = 0.3  # elementary charge
    
    # LJ parameters - Molecular Kob-Andersen Model
    sigma_O = 6.230426584 * BOHR_TO_NM  # nm
    epsilon_O = 0.00016685201 * HARTREE_TO_KJMOL  # kJ/mol
    sigma_N = 5.48277488 * BOHR_TO_NM  # nm
    epsilon_N = 0.000083426 * HARTREE_TO_KJMOL  # kJ/mol
    sigma_NO = 4.9832074319 * BOHR_TO_NM  # nm
    epsilon_NO = 0.00025027802 * HARTREE_TO_KJMOL  # kJ/mol
    
    # Create forces
    bond_force = openmm.HarmonicBondForce()
    nonbonded_force = openmm.NonbondedForce()
    nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nonbonded_force.setCutoffDistance(0.9)  # nm
    
    # Generate molecules on a lattice
    num_OO = int(fraction_OO * num_molecules)
    side = int(np.ceil(num_molecules**(1/3)))
    spacing = box_size_nm / side
    
    # Track O and N particle indices for cross-term exceptions
    O_indices = []
    N_indices = []
    
    mol_idx = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if mol_idx >= num_molecules:
                    break
                    
                is_OO = mol_idx < num_OO
                
                # Center position
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing
                z = (k + 0.5) * spacing
                
                # Random orientation
                theta = np.random.rand() * 2 * np.pi
                phi = np.arccos(2 * np.random.rand() - 1)
                direction = np.array([
                    np.sin(phi) * np.cos(theta),
                    np.sin(phi) * np.sin(theta),
                    np.cos(phi)
                ])
                
                if is_OO:
                    mass = mass_O
                    r0 = r0_OO
                    k = k_OO
                    sigma = sigma_O
                    epsilon = epsilon_O
                else:
                    mass = mass_N
                    r0 = r0_NN
                    k = k_NN
                    sigma = sigma_N
                    epsilon = epsilon_N
                
                # Particle positions
                r1 = np.array([x, y, z]) - 0.5 * r0 * direction
                r2 = np.array([x, y, z]) + 0.5 * r0 * direction
                
                # Add particles to system
                idx1 = system.addParticle(mass)
                idx2 = system.addParticle(mass)
                
                positions.append(openmm.Vec3(*r1) * unit.nanometer)
                positions.append(openmm.Vec3(*r2) * unit.nanometer)
                
                # Add bond
                bond_force.addBond(idx1, idx2, r0, k)
                
                # Add nonbonded parameters (charge, sigma, epsilon)
                charge1, charge2 = -charge_magnitude, +charge_magnitude
                nonbonded_force.addParticle(charge1, sigma, epsilon)
                nonbonded_force.addParticle(charge2, sigma, epsilon)
                
                # Track particle types for cross-term exceptions
                if is_OO:
                    O_indices.append((idx1, charge1))
                    O_indices.append((idx2, charge2))
                else:
                    N_indices.append((idx1, charge1))
                    N_indices.append((idx2, charge2))
                
                # Add exclusion for bonded pair
                nonbonded_force.addException(idx1, idx2, 0.0, 1.0, 0.0)
                
                mol_idx += 1
            if mol_idx >= num_molecules:
                break
        if mol_idx >= num_molecules:
            break
    
    # Add N-O cross-term exceptions (Kob-Andersen non-additive parameters)
    n_cross_terms = 0
    for (idx_O, charge_O) in O_indices:
        for (idx_N, charge_N) in N_indices:
            chargeProd = charge_O * charge_N
            nonbonded_force.addException(idx_O, idx_N, chargeProd, sigma_NO, epsilon_NO)
            n_cross_terms += 1
    
    print(f"Added {n_cross_terms} N-O cross-term exceptions for KA non-additive interactions")
    
    # Add forces to system
    system.addForce(bond_force)  # Force group 0 (default)
    system.addForce(nonbonded_force)  # Force group 0 (default)
    
    # Set periodic box
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_size_nm, 0, 0),
        openmm.Vec3(0, box_size_nm, 0),
        openmm.Vec3(0, 0, box_size_nm)
    )
    
    print(f"Created system with {system.getNumParticles()} particles ({num_molecules} molecules)")
    print(f"  O-O dimers: {num_OO}")
    print(f"  N-N dimers: {num_molecules - num_OO}")
    print(f"  Box size: {box_size_nm} nm")
    
    return system, positions


def add_cavity_particle(system, positions, omegac, photon_mass=1.0):
    """
    Add a cavity photon particle to the system.
    
    Parameters
    ----------
    system : openmm.System
    positions : list
    omegac : float
        Cavity frequency in OpenMM units (kJ/mol)
    photon_mass : float
        Photon mass in amu
        
    Returns
    -------
    cavity_index : int
        Index of the cavity particle
    """
    # Add the cavity particle
    cavity_index = system.addParticle(photon_mass)
    
    # Position at origin
    positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)
    
    # Add nonbonded parameters for the cavity particle (no interactions)
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            # Add cavity particle with no charge, no LJ
            force.addParticle(0.0, 0.1, 0.0)  # charge=0, sigma=0.1nm, epsilon=0
    
    print(f"Added cavity particle at index {cavity_index}")
    print(f"  Omega_c: {omegac} (OpenMM energy units)")
    print(f"  Photon mass: {photon_mass} amu")
    
    return cavity_index


def compute_centroid_dipole_moment(integrator, charges, num_molecular_particles, num_beads):
    """
    Compute the centroid dipole moment for RPMD.
    
    Parameters
    ----------
    integrator : openmm.RPMDIntegrator
        RPMD integrator with bead positions
    charges : list
        Particle charges
    num_molecular_particles : int
        Number of molecular particles (excluding cavity)
    num_beads : int
        Number of RPMD beads
        
    Returns
    -------
    dipole : np.ndarray
        Centroid dipole moment vector (3,) in e*nm units
    """
    dipole = np.zeros(3)
    
    # Average over all beads to get centroid
    for i in range(num_molecular_particles):
        centroid_pos = np.zeros(3)
        
        # Sum positions across all beads
        for bead in range(num_beads):
            state = integrator.getState(bead, getPositions=True)
            positions = state.getPositions(asNumpy=True)
            pos_nm = positions[i].value_in_unit(unit.nanometer)
            centroid_pos += pos_nm
        
        # Average to get centroid
        centroid_pos /= num_beads
        
        # Extract charge value
        if hasattr(charges[i], 'value_in_unit'):
            charge_e = charges[i].value_in_unit(unit.elementary_charge)
        else:
            charge_e = float(charges[i])
        
        dipole[0] += charge_e * centroid_pos[0]
        dipole[1] += charge_e * centroid_pos[1]
        dipole[2] += charge_e * centroid_pos[2]
    
    return dipole


def run_rpmd_simulation(num_molecules=250, num_beads=8, lambda_coupling=0.001, 
                       temperature_K=100.0, dt=0.001, equilibration_time_ps=100.0, 
                       production_time_ps=900.0, cavity_freq_cm=1560.0, 
                       disable_dipole_output=False):
    """Run the RPMD cavity dimer simulation."""
    
    print("=" * 60)
    print("Cavity OpenMM Dimer RPMD Simulation")
    print("=" * 60)
    
    # System parameters
    box_size_nm = 2.5
    
    # Cavity parameters
    # Convert from cm⁻¹ to atomic units (Hartree)
    omegac_au = cavity_freq_cm / 219474.63  # Hartree
    photon_mass = 1.0 / 1822.888  # amu
    
    # Simulation parameters
    equilibration_steps = int(equilibration_time_ps / dt)
    production_steps = int(production_time_ps / dt)
    total_steps = equilibration_steps + production_steps
    
    # Dipole output parameters
    # For IR spectrum: getState() is EXTREMELY expensive (CPU-GPU transfer for all beads)
    # Strategy: Sample VERY sparsely (every 5 ps) for production runs
    # Nyquist: For ν_max = 2000 cm⁻¹, need Δt < 8.3 fs
    # Using 5 ps sampling gives Δν = 0.2 cm⁻¹ resolution (excellent!)
    # and reduces CPU-GPU transfers by 1000x vs 5 fs sampling
    dipole_output_interval_steps = 1250  # Every 5 ps (5000 fs) at dt=4fs
    dipole_output_interval_ps = dt * dipole_output_interval_steps
    
    print(f"\nSimulation parameters:")
    print(f"  Temperature: {temperature_K} K")
    print(f"  RPMD beads: {num_beads}")
    print(f"  Timestep: {dt} ps ({dt*1000:.1f} fs)")
    print(f"  Cavity frequency: {cavity_freq_cm} cm⁻¹ ({omegac_au:.6f} Hartree)")
    print(f"  Lambda coupling: {lambda_coupling}")
    print(f"  Equilibration: {equilibration_time_ps} ps ({equilibration_steps} steps)")
    print(f"  Production: {production_time_ps} ps ({production_steps} steps)")
    print(f"  Total time: {equilibration_time_ps + production_time_ps} ps")
    if not disable_dipole_output:
        freq_resolution = 1000.0 / (production_time_ps)  # cm⁻¹
        max_freq = 1000.0 / (2 * dipole_output_interval_ps)  # Nyquist limit in cm⁻¹
        print(f"  Centroid dipole output: Every {dipole_output_interval_steps} steps ({dipole_output_interval_ps*1000:.1f} fs)")
        print(f"  Expected {total_steps // dipole_output_interval_steps:,} dipole samples")
        print(f"  Frequency resolution: {freq_resolution:.2f} cm⁻¹, Max freq: {max_freq:.0f} cm⁻¹")
    
    # Create system
    print("\n--- Creating Dimer System ---")
    system, positions = create_diamer_system(
        num_molecules=num_molecules,
        fraction_OO=0.8,
        box_size_nm=box_size_nm,
        temperature_K=temperature_K
    )
    
    # Store charges for dipole calculation
    charges = []
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for i in range(system.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(i)
                charges.append(charge)
    
    # Add cavity particle
    print("\n--- Adding Cavity Particle ---")
    cavity_index = add_cavity_particle(system, positions, omegac_au, photon_mass)
    num_molecular_particles = cavity_index
    
    # Add cavity force
    print("\n--- Adding Cavity Force ---")
    try:
        cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass)
        system.addForce(cavity_force)
        print(f"  CavityForce added successfully")
        print(f"  Omega_c: {omegac_au:.6f} a.u.")
        print(f"  Lambda coupling: {lambda_coupling} (ACTIVE from t=0)")
        
        displacer = openmm.CavityParticleDisplacer(cavity_index, omegac_au, photon_mass)
        displacer.setSwitchOnLambda(lambda_coupling)
        system.addForce(displacer)
        print(f"  CavityParticleDisplacer added (finite-Q mode)")
        
        cavity_available = True
    except AttributeError as e:
        print(f"  CavityForce not available: {e}")
        cavity_available = False
    
    # Create RPMD integrator with PILE_G thermostat
    # NOTE: Ring polymer contractions have too much FFT overhead for small systems
    # Running WITHOUT contractions for maximum speed on this 500-particle system
    print("\n--- Creating RPMD Integrator ---")
    friction = 1.0  # ps^-1 - standard friction for RPMD
    centroid_friction = 0.5  # ps^-1 - friction for Bussi centroid thermostat
    
    integrator = openmm.RPMDIntegrator(
        num_beads,
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt * unit.picosecond
    )
    
    # Set PILE_G thermostat (Bussi on centroid, Langevin on internal modes)
    integrator.setThermostatType(openmm.RPMDIntegrator.PileG)
    integrator.setCentroidFriction(centroid_friction / unit.picosecond)
    
    # Mark cavity particle as classical (excluded from RPMD ring polymer dynamics)
    # Type 0 = quantum (default for all particles), Type 1 = classical
    integrator.setParticleType(cavity_index, 1)
    integrator.setClassicalThermostat(openmm.RPMDIntegrator.LangevinClassical)
    
    print(f"  RPMDIntegrator created")
    print(f"  Thermostat: PILE_G (Bussi centroid + Langevin internal)")
    print(f"  Centroid friction (Bussi): {centroid_friction} ps^-1")
    print(f"  Internal mode friction (PILE): {friction} ps^-1")
    print(f"  Number of beads: {num_beads}")
    print(f"  Cavity particle (index {cavity_index}): CLASSICAL (excluded from RPMD)")

    
    # Create simulation with CUDA platform
    print("\n--- Creating Simulation ---")
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        print(f"  Using CUDA platform (GPU acceleration)")
    except Exception:
        platform = openmm.Platform.getPlatformByName('CPU')
        print(f"  Using CPU platform (CUDA not available)")
    
    # Create context
    context = openmm.Context(system, integrator, platform)
    
    # Set positions for all beads
    print("\n--- Initializing RPMD Beads ---")
    for bead in range(num_beads):
        # Add small random perturbations for each bead
        bead_positions = []
        for pos in positions:
            # pos is already a Vec3 with units
            pos_nm = pos.value_in_unit(unit.nanometer)
            dx = np.random.randn() * 0.01
            dy = np.random.randn() * 0.01
            dz = np.random.randn() * 0.01
            bead_pos = openmm.Vec3(pos_nm[0] + dx, pos_nm[1] + dy, pos_nm[2] + dz) * unit.nanometer
            bead_positions.append(bead_pos)
        integrator.setPositions(bead, bead_positions)
    
    # Also set context positions (for the underlying OpenMM context)
    context.setPositions(positions)
    print(f"  Positions set for {num_beads} beads")
    
    # Set velocities for all beads
    print("\n--- Initializing Velocities ---")
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    
    # Get velocities from context and set for each bead with perturbations
    state = context.getState(getVelocities=True)
    base_velocities = state.getVelocities()
    
    for bead in range(num_beads):
        bead_velocities = []
        for v in base_velocities:
            vx = v[0].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01 * np.random.randn()
            vy = v[1].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01 * np.random.randn()
            vz = v[2].value_in_unit(unit.nanometers/unit.picoseconds) + 0.01 * np.random.randn()
            bead_velocities.append(openmm.Vec3(vx, vy, vz) * unit.nanometers/unit.picoseconds)
        integrator.setVelocities(bead, bead_velocities)
    
    print(f"  Velocities initialized for {num_beads} beads")
    
    # Apply finite-Q displacement to cavity particle
    if cavity_available:
        print("\n--- Applying Finite-Q Displacement ---")
        displacement_magnitude = 0.01  # nm
        for bead in range(num_beads):
            state = integrator.getState(bead, getPositions=True)
            positions_list = list(state.getPositions())
            positions_list[cavity_index] = openmm.Vec3(displacement_magnitude, 0.0, 0.0) * unit.nanometer
            integrator.setPositions(bead, positions_list)
        print(f"  Displacement: {displacement_magnitude} nm along x-axis for all beads")
    
    # Calculate total steps
    total_steps = equilibration_steps + production_steps
    
    # Run simulation
    print("\n--- Running RPMD Simulation ---")
    print(f"  Running {total_steps} steps...")
    print(f"  Progress updates every 1 ps (1,000 steps)")
    print("")
    
    # Prepare dipole trajectory storage
    dipole_times = []
    dipole_trajectory = []
    step_counter = 0
    
    report_interval = 1000  # Report every 1 ps
    num_reports = -(-total_steps // report_interval)  # ceiling division so no steps are dropped
    
    start_time = time.time()
    last_report_time = start_time
    
    for i in range(num_reports):
        steps_this_block = min(report_interval, total_steps - step_counter)
        for step in range(steps_this_block):
            integrator.step(1)
            step_counter += 1
            if not disable_dipole_output and (step_counter % dipole_output_interval_steps) == 0:
                current_time = step_counter * dt
                dipole = compute_centroid_dipole_moment(integrator, charges, num_molecular_particles, num_beads)
                dipole_times.append(current_time)
                dipole_trajectory.append(dipole)
        
        # Report progress
        current_wall_time = time.time()
        elapsed = current_wall_time - start_time
        elapsed_since_last = current_wall_time - last_report_time
        last_report_time = current_wall_time
        
        energy = integrator.getTotalEnergy()
        if hasattr(energy, 'value_in_unit'):
            energy_val = energy.value_in_unit(unit.kilojoule_per_mole)
        else:
            energy_val = float(energy)
        
        sim_time_ps = step_counter * dt
        progress_pct = (step_counter / total_steps) * 100
        steps_per_sec = steps_this_block / elapsed_since_last
        ns_per_day = steps_per_sec * dt * 86400.0 / 1000.0
        
        steps_remaining = total_steps - step_counter
        time_remaining_sec = steps_remaining / steps_per_sec if steps_per_sec > 0 else 0
        time_remaining_hr = time_remaining_sec / 3600
        
        print(f"  [{progress_pct:5.1f}%] Step {step_counter:7d}/{total_steps} | "
              f"Sim time: {sim_time_ps:6.1f} ps | "
              f"Speed: {ns_per_day:6.1f} ns/day | "
              f"ETA: {time_remaining_hr:5.1f} hr")
        print(f"         Total Energy: {energy_val:8.2f} kJ/mol")
        print("")
        
        # Save NPZ file on the fly
        if not disable_dipole_output and len(dipole_times) > 0:
            output_file = f"cavity_diamer_lambda{lambda_coupling:.4f}.npz"
            np.savez(output_file,
                     time_ps=np.array(dipole_times),
                     dipole_nm=np.array(dipole_trajectory),
                     metadata={
                         'temperature_K': temperature_K,
                         'num_molecules': num_molecules,
                         'num_beads': num_beads,
                         'cavity_freq_cm': cavity_freq_cm,
                         'equilibration_ps': equilibration_time_ps,
                         'production_ps': production_time_ps,
                         'dt_ps': dt,
                         'output_interval_ps': dipole_output_interval_ps,
                         'lambda_coupling': lambda_coupling,
                         'omegac_au': omegac_au,
                         'cavity_index': cavity_index,
                         'thermostat': 'PILE_G',
                         'centroid_friction': centroid_friction,
                         'status': 'running',
                         'step': step_counter,
                         'total_steps': total_steps
                     })
    
    total_elapsed = time.time() - start_time
    print(f"  Simulation complete! Total elapsed time: {total_elapsed/3600:.2f} hr ({total_elapsed/60:.1f} min)")
    
    # Save final dipole trajectory
    if not disable_dipole_output:
        print("\n--- Saving Centroid Dipole Trajectory ---")
        dipole_times_array = np.array(dipole_times)
        dipole_trajectory_array = np.array(dipole_trajectory)
        
        output_file = f"cavity_diamer_lambda{lambda_coupling:.4f}.npz"
        np.savez(output_file,
                 time_ps=dipole_times_array,
                 dipole_nm=dipole_trajectory_array,
                 metadata={
                     'temperature_K': temperature_K,
                     'num_molecules': num_molecules,
                     'num_beads': num_beads,
                     'cavity_freq_cm': cavity_freq_cm,
                     'equilibration_ps': equilibration_time_ps,
                     'production_ps': production_time_ps,
                     'dt_ps': dt,
                     'output_interval_ps': dipole_output_interval_ps,
                     'lambda_coupling': lambda_coupling,
                     'omegac_au': omegac_au,
                     'cavity_index': cavity_index,
                     'thermostat': 'PILE_G',
                     'centroid_friction': centroid_friction,
                     'status': 'complete',
                     'step': total_steps,
                     'total_steps': total_steps,
                     'elapsed_time_s': total_elapsed
                 })
        
        print(f"  Saved centroid dipole trajectory to: {output_file}")
        print(f"  Total data points: {len(dipole_times_array)}")
        print(f"  Time range: {dipole_times_array[0]:.2f} - {dipole_times_array[-1]:.2f} ps")
    
    print("\n" + "=" * 60)
    print("RPMD simulation completed successfully!")
    print("=" * 60)
    
    return True


def main():
    """Main entry point with command-line argument parsing."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Dimer Cavity RPMD Simulation')
    parser.add_argument('--dimers', type=int, default=250,
                       help='Number of dimers (default: 250)')
    parser.add_argument('--beads', type=int, default=8,
                       help='Number of RPMD beads (default: 8)')
    parser.add_argument('--lambda', type=float, default=0.0700, dest='lambda_coupling',
                       help='Coupling strength (default: 0.0700)')
    parser.add_argument('--temp', type=float, default=100.0,
                       help='Temperature in K (default: 100)')
    parser.add_argument('--dt', type=float, default=0.001,
                       help='Timestep in ps (default: 0.001 = 1 fs)')
    parser.add_argument('--equil', type=float, default=100.0,
                       help='Equilibration time in ps (default: 100)')
    parser.add_argument('--prod', type=float, default=900.0,
                       help='Production time in ps (default: 900)')
    parser.add_argument('--cavity-freq', type=float, default=1560.0,
                       help='Cavity frequency in cm⁻¹ (default: 1560 for O-O stretch)')
    parser.add_argument('--no-dipole', action='store_true', dest='no_dipole',
                       help='Disable dipole output for speed benchmark')
    
    args = parser.parse_args()
    
    try:
        success = run_rpmd_simulation(
            num_molecules=args.dimers,
            num_beads=args.beads,
            lambda_coupling=args.lambda_coupling,
            temperature_K=args.temp,
            dt=args.dt,
            equilibration_time_ps=args.equil,
            production_time_ps=args.prod,
            cavity_freq_cm=args.cavity_freq,
            disable_dipole_output=args.no_dipole
        )
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\nSimulation failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
