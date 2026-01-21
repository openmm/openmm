#!/usr/bin/env python3
"""
Cavity-coupled diamer simulation with strong coupling for Rabi splitting.
Version 2: With proper equilibration phase before turning on coupling.

Parameters:
- Cavity frequency: 1560 cm⁻¹ (resonant with O-O stretch)
- Lambda coupling: 0.07 (strong coupling for visible Rabi splitting)
"""

import numpy as np
import time
import sys

# OpenMM imports
try:
    from openmm import openmm
    from openmm import unit
    print("✓ OpenMM loaded successfully")
except ImportError as e:
    print(f"Error importing OpenMM: {e}")
    sys.exit(1)

# Unit conversions
BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5
HARTREE_TO_CM = 219474.63  # 1 Hartree = 219474.63 cm⁻¹

def create_diamer_system(num_molecules=50, fraction_OO=0.8, box_size_nm=2.5, 
                         temperature_K=100.0, seed=42):
    """Create a system of O-O and N-N dimers."""
    
    np.random.seed(seed)
    system = openmm.System()
    
    # Set periodic box
    box_vec = openmm.Vec3(box_size_nm, 0, 0) * unit.nanometer
    system.setDefaultPeriodicBoxVectors(box_vec, 
                                        openmm.Vec3(0, box_size_nm, 0) * unit.nanometer,
                                        openmm.Vec3(0, 0, box_size_nm) * unit.nanometer)
    
    # Masses
    mass_O = 16.0  # amu
    mass_N = 14.0  # amu
    
    # Bond parameters - CORRECT values with 2x factor for OpenMM
    k_OO_au = 2 * 0.73204  # Hartree/Bohr² (doubled for OpenMM convention)
    r0_OO_au = 2.281655158  # Bohr
    k_OO = k_OO_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
    r0_OO = r0_OO_au * BOHR_TO_NM
    
    k_NN_au = 2 * 1.4325  # Hartree/Bohr² (doubled for OpenMM)
    r0_NN_au = 2.0743522177  # Bohr
    k_NN = k_NN_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
    r0_NN = r0_NN_au * BOHR_TO_NM
    
    # Charges
    charge_magnitude = 0.3  # elementary charge
    
    # Create forces
    bond_force = openmm.HarmonicBondForce()
    nonbonded_force = openmm.NonbondedForce()
    nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nonbonded_force.setCutoffDistance(1.0 * unit.nanometer)
    
    positions = []
    charges = []
    
    n_OO = int(num_molecules * fraction_OO)
    n_NN = num_molecules - n_OO
    
    # Place dimers
    for i in range(num_molecules):
        is_OO = i < n_OO
        mass = mass_O if is_OO else mass_N
        r0 = r0_OO if is_OO else r0_NN
        k = k_OO if is_OO else k_NN
        sigma = 0.3 if is_OO else 0.28
        epsilon = 0.5 if is_OO else 0.4
        
        # Random position in box
        center = np.random.uniform(0.2, box_size_nm - 0.2, 3)
        
        # Random orientation
        theta = np.random.uniform(0, 2*np.pi)
        phi = np.random.uniform(0, np.pi)
        direction = np.array([np.sin(phi)*np.cos(theta), 
                             np.sin(phi)*np.sin(theta), 
                             np.cos(phi)])
        
        # Two atoms
        pos1 = center - 0.5 * r0 * direction
        pos2 = center + 0.5 * r0 * direction
        
        idx1 = system.addParticle(mass)
        idx2 = system.addParticle(mass)
        
        positions.append(openmm.Vec3(*pos1) * unit.nanometer)
        positions.append(openmm.Vec3(*pos2) * unit.nanometer)
        
        # Add bond
        bond_force.addBond(idx1, idx2, r0, k)
        
        # Add nonbonded - opposite charges for dipole
        nonbonded_force.addParticle(-charge_magnitude, sigma, epsilon)
        nonbonded_force.addParticle(+charge_magnitude, sigma, epsilon)
        charges.append(-charge_magnitude)
        charges.append(+charge_magnitude)
        
        # Exclude bonded pairs from nonbonded
        nonbonded_force.addException(idx1, idx2, 0.0, 1.0, 0.0)
    
    system.addForce(bond_force)
    system.addForce(nonbonded_force)
    
    print(f"Created system with {system.getNumParticles()} particles ({num_molecules} molecules)")
    print(f"  O-O dimers: {n_OO}")
    print(f"  N-N dimers: {n_NN}")
    print(f"  Box size: {box_size_nm} nm")
    
    return system, positions, charges


def compute_dipole_moment(state, charges, num_molecular_particles):
    """Compute total dipole moment of molecular particles."""
    positions = state.getPositions(asNumpy=True)
    dipole = np.zeros(3)
    
    for i in range(num_molecular_particles):
        pos_nm = positions[i].value_in_unit(unit.nanometer)
        charge_e = float(charges[i])
        dipole[0] += charge_e * pos_nm[0]
        dipole[1] += charge_e * pos_nm[1]
        dipole[2] += charge_e * pos_nm[2]
    
    return dipole


def run_simulation():
    """Run cavity-coupled simulation with strong coupling for Rabi splitting."""
    
    print("=" * 70)
    print("Cavity Simulation for Rabi Splitting (v2 - with equilibration)")
    print("=" * 70)
    
    # System parameters
    num_molecules = 50
    box_size_nm = 2.5
    temperature_K = 100.0
    
    # CAVITY PARAMETERS FOR RABI SPLITTING
    target_freq_cm = 1560  # cm⁻¹ - resonant with O-O stretch!
    omegac_au = target_freq_cm / HARTREE_TO_CM  # ~0.00711 a.u.
    lambda_coupling = 0.07  # STRONG coupling for visible Rabi splitting
    photon_mass = 1.0  # amu
    
    print(f"\nCavity parameters:")
    print(f"  Cavity frequency: {target_freq_cm} cm⁻¹ (resonant with O-O stretch)")
    print(f"  ω_c (a.u.): {omegac_au:.6f}")
    print(f"  λ coupling: {lambda_coupling} (STRONG)")
    print(f"  Photon mass: {photon_mass} amu")
    
    # Simulation parameters
    dt = 0.001  # ps (1 fs)
    equilibration_time_ps = 100.0  # 100 ps equilibration with coupling OFF
    production_time_ps = 400.0     # 400 ps production with coupling ON
    
    equilibration_steps = int(equilibration_time_ps / dt)
    production_steps = int(production_time_ps / dt)
    total_steps = equilibration_steps + production_steps
    
    print(f"\nSimulation parameters:")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt*1000:.0f} fs")
    print(f"  Equilibration: {equilibration_time_ps} ps (coupling OFF)")
    print(f"  Production: {production_time_ps} ps (coupling ON)")
    print(f"  Total: {equilibration_time_ps + production_time_ps} ps")
    
    # Create system
    print("\n--- Creating Diamer System ---")
    system, positions, charges = create_diamer_system(
        num_molecules=num_molecules,
        box_size_nm=box_size_nm,
        temperature_K=temperature_K
    )
    num_molecular_particles = len(positions)
    
    # Add cavity particle
    print("\n--- Adding Cavity Particle ---")
    cavity_index = system.addParticle(photon_mass)
    positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)
    print(f"Cavity particle added at index {cavity_index}")
    
    # Add cavity force - START WITH COUPLING OFF
    print("\n--- Adding Cavity Force ---")
    try:
        # Start with lambda=0, turn on after equilibration
        cavity_force = openmm.CavityForce(cavity_index, omegac_au, 0.0, photon_mass)
        # Schedule coupling to turn on after equilibration
        cavity_force.setCouplingOnStep(equilibration_steps, lambda_coupling)
        system.addForce(cavity_force)
        print(f"  CavityForce added with initial λ = 0")
        print(f"  Coupling will turn ON at step {equilibration_steps} with λ = {lambda_coupling}")
    except AttributeError as e:
        print(f"  ERROR: CavityForce not available: {e}")
        return False
    
    # Add nonbonded exception for cavity
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.addParticle(0.0, 0.1, 0.0)  # No charge, no LJ for cavity
    
    # Add Bussi thermostat for molecules only
    print("\n--- Adding Bussi Thermostat ---")
    try:
        tau = 0.1  # ps
        bussi = openmm.BussiThermostat(temperature_K, tau)
        bussi.setApplyToAllParticles(False)
        for i in range(num_molecular_particles):
            bussi.addParticle(i)
        system.addForce(bussi)
        print(f"  Bussi thermostat added for {num_molecular_particles} particles")
    except AttributeError:
        print("  Bussi thermostat not available")
    
    # Create integrator with low friction
    print("\n--- Creating Integrator ---")
    friction = 0.01  # ps^-1 - very low to preserve vibrational dynamics
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt * unit.picoseconds
    )
    print(f"  LangevinMiddleIntegrator with friction = {friction} ps⁻¹")
    
    # Create simulation
    print("\n--- Creating Simulation ---")
    platform = openmm.Platform.getPlatformByName('Reference')
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    print("  Using Reference platform")
    
    # Energy minimization
    print("\n--- Energy Minimization ---")
    state = context.getState(getEnergy=True)
    print(f"  Initial PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=1000)
    state = context.getState(getEnergy=True)
    print(f"  Final PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")
    
    # Storage for dipole trajectory (only during production)
    dipole_times = []
    dipole_values = []
    
    # Run simulation
    print("\n--- Running Simulation ---")
    print(f"  Total steps: {total_steps}")
    
    start_time = time.time()
    report_interval = 10000  # Report every 10 ps
    save_interval = 1  # Save dipole every step (1 fs) during production
    
    for step in range(total_steps):
        integrator.step(1)
        
        # Save dipole only during production phase
        in_production = (step >= equilibration_steps)
        if in_production and (step % save_interval == 0):
            state = context.getState(getPositions=True)
            dipole = compute_dipole_moment(state, charges, num_molecular_particles)
            dipole_times.append((step - equilibration_steps) * dt)  # Time from start of production
            dipole_values.append(dipole.copy())
        
        # Progress report
        if (step + 1) % report_interval == 0:
            elapsed = time.time() - start_time
            speed = (step + 1) / elapsed
            remaining = (total_steps - step - 1) / speed
            sim_time = (step + 1) * dt
            state = context.getState(getEnergy=True, getPositions=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            
            # Get cavity position
            pos = state.getPositions(asNumpy=True)
            cavity_pos = pos[cavity_index].value_in_unit(unit.nanometer)
            
            phase = "EQUILIBRATION" if step < equilibration_steps else "PRODUCTION"
            print(f"  [{100*(step+1)/total_steps:5.1f}%] t={sim_time:6.1f} ps | {phase} | "
                  f"PE={pe:7.1f} kJ/mol | Cavity z: {cavity_pos[2]:.4f} nm | "
                  f"ETA: {remaining/60:.1f} min")
    
    total_time = time.time() - start_time
    print(f"\n  Simulation complete in {total_time/60:.1f} min")
    
    # Save data
    print("\n--- Saving Data ---")
    dipole_times = np.array(dipole_times)
    dipole_values = np.array(dipole_values)
    
    output_file = 'rabi_splitting_dipole.npz'
    metadata = {
        'temperature_K': temperature_K,
        'num_molecules': num_molecules,
        'equilibration_ps': equilibration_time_ps,
        'production_ps': production_time_ps,
        'dt_ps': dt,
        'lambda_coupling': lambda_coupling,
        'cavity_freq_cm': target_freq_cm,
        'omegac_au': omegac_au
    }
    
    np.savez(output_file,
             time_ps=dipole_times,
             dipole_nm=dipole_values,
             metadata=metadata)
    
    print(f"  Saved: {output_file}")
    print(f"  Data points: {len(dipole_times)}")
    print(f"  Time range: {dipole_times[0]:.2f} - {dipole_times[-1]:.2f} ps")
    
    print("\n" + "=" * 70)
    print("Simulation complete!")
    print("=" * 70)
    
    return True


if __name__ == "__main__":
    success = run_simulation()
    sys.exit(0 if success else 1)
