#!/usr/bin/env python3
"""
Cavity-coupled diamer simulation with strong coupling for Rabi splitting.
Final version with lattice placement and proper finite-q equilibration.

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
HARTREE_TO_CM = 219474.63


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
    print("Cavity Simulation for Rabi Splitting")
    print("(with lattice placement and equilibration)")
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
    equilibration_time_ps = 50.0   # 50 ps equilibration with coupling OFF
    production_time_ps = 200.0     # 200 ps production with coupling ON
    
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
    system = openmm.System()
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_size_nm, 0, 0) * unit.nanometer,
        openmm.Vec3(0, box_size_nm, 0) * unit.nanometer,
        openmm.Vec3(0, 0, box_size_nm) * unit.nanometer
    )
    
    # Bond parameters
    k_OO = 1374026  # kJ/(mol*nm²) - correct for O-O stretch at ~1560 cm⁻¹
    r0_OO = 0.1207  # nm
    charge_magnitude = 0.3  # e
    mass_O = 16.0  # amu
    
    positions = []
    charges = []
    bond_force = openmm.HarmonicBondForce()
    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(1.0 * unit.nanometer)
    
    # LATTICE placement
    np.random.seed(42)
    side = int(np.ceil(num_molecules**(1/3)))
    spacing = box_size_nm / side
    print(f"Lattice: {side}x{side}x{side}, spacing = {spacing:.3f} nm")
    
    mol_idx = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if mol_idx >= num_molecules:
                    break
                
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing
                z = (k + 0.5) * spacing
                
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                direction = np.array([np.sin(phi)*np.cos(theta), 
                                      np.sin(phi)*np.sin(theta), 
                                      np.cos(phi)])
                
                pos1 = np.array([x, y, z]) - 0.5 * r0_OO * direction
                pos2 = np.array([x, y, z]) + 0.5 * r0_OO * direction
                
                idx1 = system.addParticle(mass_O)
                idx2 = system.addParticle(mass_O)
                
                positions.append(openmm.Vec3(*pos1) * unit.nanometer)
                positions.append(openmm.Vec3(*pos2) * unit.nanometer)
                
                bond_force.addBond(idx1, idx2, r0_OO, k_OO)
                nonbonded.addParticle(-charge_magnitude, 0.3, 0.5)
                nonbonded.addParticle(+charge_magnitude, 0.3, 0.5)
                nonbonded.addException(idx1, idx2, 0.0, 1.0, 0.0)
                charges.extend([-charge_magnitude, +charge_magnitude])
                
                mol_idx += 1
            if mol_idx >= num_molecules:
                break
        if mol_idx >= num_molecules:
            break
    
    system.addForce(bond_force)
    system.addForce(nonbonded)
    num_molecular_particles = len(positions)
    print(f"Created {mol_idx} molecules ({num_molecular_particles} particles)")
    
    # Add cavity particle at origin (will be displaced by displacer)
    print("\n--- Adding Cavity ---")
    cavity_index = system.addParticle(photon_mass)
    positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)
    nonbonded.addParticle(0.0, 0.1, 0.0)
    print(f"Cavity particle added at index {cavity_index}")
    
    # CavityForce: start with coupling OFF, turn on after equilibration
    cavity_force = openmm.CavityForce(cavity_index, omegac_au, 0.0, photon_mass)
    cavity_force.setCouplingOnStep(equilibration_steps, lambda_coupling)
    system.addForce(cavity_force)
    print(f"CavityForce: λ=0 initially, will switch to λ={lambda_coupling} at step {equilibration_steps}")
    
    # CavityParticleDisplacer for finite-q correction
    displacer = openmm.CavityParticleDisplacer(cavity_index, omegac_au, photon_mass)
    displacer.setSwitchOnStep(equilibration_steps)
    displacer.setSwitchOnLambda(lambda_coupling)
    system.addForce(displacer)
    print("CavityParticleDisplacer added for finite-q correction")
    
    # Create integrator
    friction = 0.01  # ps^-1 - very low to preserve vibrational dynamics
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt * unit.picoseconds
    )
    print(f"LangevinMiddleIntegrator with friction = {friction} ps⁻¹")
    
    # Create context
    platform = openmm.Platform.getPlatformByName('CPU')
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    print(f"Using {platform.getName()} platform")
    
    # Energy minimization
    print("\n--- Energy Minimization ---")
    state = context.getState(getEnergy=True)
    print(f"Initial PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=1000)
    state = context.getState(getEnergy=True)
    print(f"Final PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")
    
    # Storage for dipole trajectory (only during production)
    dipole_times = []
    dipole_values = []
    
    # Run simulation
    print("\n--- Running Simulation ---")
    print(f"Total steps: {total_steps}")
    
    start_time = time.time()
    report_interval = 10000  # Report every 10 ps
    save_interval = 1  # Save dipole every step during production
    
    for step in range(total_steps):
        integrator.step(1)
        
        # Save dipole only during production phase
        in_production = (step >= equilibration_steps)
        if in_production and (step % save_interval == 0):
            state = context.getState(getPositions=True)
            dipole = compute_dipole_moment(state, charges, num_molecular_particles)
            dipole_times.append((step - equilibration_steps) * dt)
            dipole_values.append(dipole.copy())
        
        # Progress report
        if (step + 1) % report_interval == 0:
            elapsed = time.time() - start_time
            speed = (step + 1) / elapsed
            remaining = (total_steps - step - 1) / speed
            sim_time = (step + 1) * dt
            state = context.getState(getEnergy=True, getPositions=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            
            pos = state.getPositions(asNumpy=True)
            cavity_pos = pos[cavity_index].value_in_unit(unit.nanometer)
            
            phase = "EQUIL" if step < equilibration_steps else "PROD"
            print(f"[{100*(step+1)/total_steps:5.1f}%] t={sim_time:6.1f} ps | {phase} | "
                  f"PE={pe:7.1f} kJ/mol | Cavity: ({cavity_pos[0]:+.2f}, {cavity_pos[1]:+.2f}, {cavity_pos[2]:+.2f}) | "
                  f"ETA: {remaining/60:.1f} min")
    
    total_time = time.time() - start_time
    print(f"\nSimulation complete in {total_time/60:.1f} min")
    
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
        'omegac_au': omegac_au,
        'finite_q_correction': True
    }
    
    np.savez(output_file,
             time_ps=dipole_times,
             dipole_nm=dipole_values,
             metadata=metadata)
    
    print(f"Saved: {output_file}")
    print(f"Data points: {len(dipole_times)}")
    print(f"Time range: {dipole_times[0]:.2f} - {dipole_times[-1]:.2f} ps")
    
    print("\n" + "=" * 70)
    print("Run IR spectrum analysis with:")
    print("  python3 calculate_ir_spectrum.py rabi_splitting_dipole.npz")
    print("=" * 70)
    
    return True


if __name__ == "__main__":
    success = run_simulation()
    sys.exit(0 if success else 1)
