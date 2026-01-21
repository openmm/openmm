#!/usr/bin/env python3
"""
Cavity Performance Test - Compare WITH and WITHOUT Cavity Particle
===================================================================

This test measures performance WITH a cavity particle to compare against
the baseline (no cavity) test.

System:
- 250 O-O and N-N dimers (500 atoms + 1 cavity particle = 501 total)
- Harmonic bonds
- Lennard-Jones interactions  
- Coulomb interactions
- Bussi thermostat for molecules
- Langevin for cavity particle
- CavityForce coupling

Goal: Measure performance and compare to baseline (~0.7 μs/day expected).
"""

import sys
import numpy as np
import time

try:
    from openmm import openmm
    from openmm import unit
    print("✓ OpenMM loaded successfully")
except ImportError as e:
    print(f"Error importing OpenMM: {e}")
    sys.exit(1)

# Physical constants
BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5
# Time conversion: 1 a.u. of time = 0.02418884254 ps
AU_TIME_TO_PS = 0.02418884254

def create_diamer_system(num_molecules=250, fraction_OO=0.8, box_size_nm=4.0, 
                         temperature_K=100.0, seed=42):
    """Create a system of O-O and N-N dimers."""
    np.random.seed(seed)
    
    system = openmm.System()
    positions = []
    
    # Particle masses
    mass_O = 16.0  # amu
    mass_N = 14.0  # amu
    
    # Bond parameters
    k_OO_au = 0.732  # Hartree/Bohr²
    r0_OO_au = 2.28  # Bohr
    k_OO = k_OO_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
    r0_OO = r0_OO_au * BOHR_TO_NM
    
    k_NN_au = 1.43
    r0_NN_au = 2.07
    k_NN = k_NN_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
    r0_NN = r0_NN_au * BOHR_TO_NM
    
    # Charges
    charge_magnitude = 0.3  # elementary charge
    
    # LJ parameters
    sigma_O = 0.3  # nm
    epsilon_O = 0.5  # kJ/mol
    sigma_N = 0.28
    epsilon_N = 0.4
    
    # Create forces
    bond_force = openmm.HarmonicBondForce()
    nonbonded_force = openmm.NonbondedForce()
    nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nonbonded_force.setCutoffDistance(0.9)  # nm
    
    # Generate molecules on a lattice
    num_OO = int(fraction_OO * num_molecules)
    side = int(np.ceil(num_molecules**(1/3)))
    spacing = box_size_nm / side
    
    charges = []
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
                
                # Add particles
                idx1 = system.addParticle(mass)
                idx2 = system.addParticle(mass)
                
                positions.append(openmm.Vec3(*r1) * unit.nanometer)
                positions.append(openmm.Vec3(*r2) * unit.nanometer)
                
                # Add bond
                bond_force.addBond(idx1, idx2, r0, k)
                
                # Add nonbonded parameters
                nonbonded_force.addParticle(-charge_magnitude, sigma, epsilon)
                nonbonded_force.addParticle(+charge_magnitude, sigma, epsilon)
                charges.append(-charge_magnitude)
                charges.append(+charge_magnitude)
                
                # Add exclusion for bonded pair
                nonbonded_force.addException(idx1, idx2, 0.0, 1.0, 0.0)
                
                mol_idx += 1
            if mol_idx >= num_molecules:
                break
        if mol_idx >= num_molecules:
            break
    
    # Add forces to system
    system.addForce(bond_force)
    system.addForce(nonbonded_force)
    
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
    
    return system, positions, charges


def add_cavity_particle(system, positions, omegac, photon_mass=1.0):
    """Add a cavity photon particle to the system."""
    # Add the cavity particle
    cavity_index = system.addParticle(photon_mass)
    
    # Position at origin
    positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)
    
    # Add nonbonded parameters for the cavity particle (no interactions)
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            # Add cavity particle with no charge, no LJ
            force.addParticle(0.0, 0.1, 0.0)
    
    print(f"Added cavity particle at index {cavity_index}")
    print(f"  Omega_c: {omegac} kJ/mol")
    print(f"  Photon mass: {photon_mass} amu")
    
    return cavity_index


def run_cavity_test():
    """Run performance test WITH cavity particle."""
    
    print("=" * 70)
    print("Cavity Performance Test (WITH Cavity Particle)")
    print("=" * 70)
    
    # System parameters - match baseline
    num_molecules = 250  # 500 atoms + 1 cavity = 501 total
    box_size_nm = 4.0
    temperature_K = 100.0
    
    # Cavity parameters
    omegac_au = 0.00913  # Hartree (passed in atomic units, converted internally)
    lambda_coupling = 0.001  # Small coupling
    photon_mass = 1.0  # amu
    
    # Simulation parameters - match baseline
    dt = 0.002  # ps (2 fs)
    equilibration_time_ps = 10.0
    test_time_ps = 100.0
    
    equilibration_steps = int(equilibration_time_ps / dt)
    test_steps = int(test_time_ps / dt)
    
    print(f"\nSimulation parameters:")
    print(f"  Number of molecules: {num_molecules} ({num_molecules*2} atoms + 1 cavity)")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt} ps (2 fs)")
    print(f"  Cavity frequency: {omegac_au:.6f} a.u. (Hartree)")
    AMU_TO_AU = 1822.888
    K_au = photon_mass * AMU_TO_AU * omegac_au**2  # Correct: convert mass to a.u.
    K_openmm = K_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
    print(f"  Spring constant K: {K_openmm:.0f} kJ/(mol·nm^2)")
    print(f"  Coupling: {lambda_coupling}")
    print(f"  Equilibration: {equilibration_time_ps} ps ({equilibration_steps} steps)")
    print(f"  Test run: {test_time_ps} ps ({test_steps} steps)")
    
    # Create system
    print("\n--- Creating Diamer System ---")
    system, positions, charges = create_diamer_system(
        num_molecules=num_molecules,
        fraction_OO=0.8,
        box_size_nm=box_size_nm,
        temperature_K=temperature_K
    )
    
    # Add cavity particle
    print("\n--- Adding Cavity Particle ---")
    cavity_index = add_cavity_particle(system, positions, omegac_au, photon_mass)
    num_molecular_particles = cavity_index
    
    # Add CavityForce
    print("\n--- Adding Cavity Force ---")
    try:
        cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass)
        system.addForce(cavity_force)
        print(f"  CavityForce added successfully")
        print(f"  Lambda: {lambda_coupling}")
        epsilon_au = lambda_coupling * omegac_au
        epsilon_openmm = epsilon_au * HARTREE_TO_KJMOL / BOHR_TO_NM
        print(f"  Epsilon (lambda*omega): {epsilon_openmm:.4f} kJ/(mol·nm·e)")
    except Exception as e:
        print(f"  ERROR: Could not add CavityForce: {e}")
        return False
    
    # Add Bussi thermostat for molecules only
    print("\n--- Adding Bussi Thermostat ---")
    try:
        tau = 0.1  # ps
        bussi = openmm.BussiThermostat(temperature_K, tau)
        bussi.setApplyToAllParticles(False)
        
        # Add all molecular particles (not the cavity)
        for i in range(num_molecular_particles):
            bussi.addParticle(i)
        
        system.addForce(bussi)
        print(f"  BussiThermostat added for {num_molecular_particles} particles")
        print(f"  Temperature: {temperature_K} K, Tau: {tau} ps")
    except Exception as e:
        print(f"  ERROR: Could not add BussiThermostat: {e}")
        return False
    
    # Create integrator - Langevin handles the cavity particle
    print("\n--- Creating Integrator ---")
    friction = 0.01  # ps^-1 - low friction
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt * unit.picosecond
    )
    print(f"  LangevinMiddleIntegrator created (friction: {friction} ps^-1)")
    print(f"  Bussi thermostat handles molecular particles, Langevin handles cavity")
    
    # Try to use CUDA
    print("\n--- Creating Simulation ---")
    platform_name = None
    for platform_try in ['CUDA', 'OpenCL', 'CPU']:
        try:
            platform = openmm.Platform.getPlatformByName(platform_try)
            platform_name = platform_try
            print(f"  Using {platform_name} platform")
            break
        except Exception:
            continue
    
    if platform_name is None:
        print("  ERROR: No suitable platform available!")
        return False
    
    # Set platform properties for GPU
    properties = {}
    if platform_name in ['CUDA', 'OpenCL']:
        properties['Precision'] = 'mixed'
        print(f"  Precision: mixed")
    
    # Create context
    context = openmm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    
    # Minimize energy
    print("\n--- Energy Minimization ---")
    state = context.getState(getEnergy=True)
    print(f"  Initial energy: {state.getPotentialEnergy()}")
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=100)
    state = context.getState(getEnergy=True)
    print(f"  Final energy: {state.getPotentialEnergy()}")
    
    # Set velocities
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    
    # Run equilibration
    print("\n--- Equilibration Phase ---")
    print(f"  Running {equilibration_steps} steps ({equilibration_time_ps} ps)...")
    
    start_time = time.time()
    integrator.step(equilibration_steps)
    equil_time = time.time() - start_time
    
    state = context.getState(getEnergy=True, getPositions=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    positions_now = state.getPositions()
    cavity_pos = positions_now[cavity_index]
    
    print(f"  Completed in {equil_time:.2f} seconds")
    print(f"  Final PE: {pe:.2f} kJ/mol")
    print(f"  Cavity position: ({cavity_pos.x:.4f}, {cavity_pos.y:.4f}, {cavity_pos.z:.4f}) nm")
    
    # Run test for performance measurement
    print("\n--- Performance Test Run ---")
    print(f"  Running {test_steps} steps ({test_time_ps} ps)...")
    
    report_interval = test_steps // 10  # 10 progress reports
    
    test_start_time = time.time()
    
    for i in range(10):
        integrator.step(report_interval)
        
        elapsed = time.time() - test_start_time
        steps_done = (i + 1) * report_interval
        progress_pct = (steps_done / test_steps) * 100
        
        state = context.getState(getEnergy=True, getPositions=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        positions_now = state.getPositions()
        cavity_pos = positions_now[cavity_index]
        
        # Calculate performance metrics
        steps_per_sec = steps_done / elapsed if elapsed > 0 else 0
        ns_per_day = (steps_per_sec * dt * 86400) / 1000
        us_per_day = ns_per_day / 1000
        
        print(f"  [{progress_pct:5.1f}%] Step {steps_done:6d}/{test_steps} | "
              f"PE: {pe:8.2f} kJ/mol | "
              f"Speed: {steps_per_sec:7.1f} steps/s")
        print(f"         Cavity: ({cavity_pos.x:7.4f}, {cavity_pos.y:7.4f}, {cavity_pos.z:7.4f}) nm | "
              f"Perf: {ns_per_day:6.1f} ns/day ({us_per_day:5.2f} μs/day)")
    
    total_test_time = time.time() - test_start_time
    
    # Final statistics
    print("\n" + "=" * 70)
    print("PERFORMANCE SUMMARY")
    print("=" * 70)
    
    final_steps_per_sec = test_steps / total_test_time
    final_ns_per_day = (final_steps_per_sec * dt * 86400) / 1000
    final_us_per_day = final_ns_per_day / 1000
    
    print(f"Platform: {platform_name}")
    print(f"System size: {num_molecules} molecules ({num_molecular_particles} atoms + 1 cavity)")
    print(f"Timestep: {dt} ps")
    print(f"Test duration: {test_time_ps} ps ({test_steps} steps)")
    print(f"Wall time: {total_test_time:.2f} seconds")
    print(f"")
    print(f"Speed: {final_steps_per_sec:.1f} steps/second")
    print(f"Performance: {final_ns_per_day:.1f} ns/day")
    print(f"             {final_us_per_day:.2f} μs/day")
    print(f"")
    
    # Compare to baseline expectation
    baseline_us_per_day = 0.71
    slowdown_factor = baseline_us_per_day / final_us_per_day if final_us_per_day > 0 else float('inf')
    
    print(f"Baseline (no cavity): ~{baseline_us_per_day} μs/day")
    print(f"With cavity: {final_us_per_day:.2f} μs/day")
    print(f"Slowdown factor: {slowdown_factor:.2f}x")
    print(f"")
    
    if slowdown_factor < 2.0:
        print(f"✓ PERFORMANCE GOOD (< 2x slowdown)")
    elif slowdown_factor < 5.0:
        print(f"⚠ PERFORMANCE ACCEPTABLE (2-5x slowdown)")
    else:
        print(f"✗ PERFORMANCE POOR (> 5x slowdown)")
        print(f"  This suggests the cavity force may not be running on GPU!")
    
    print("=" * 70)
    
    return True


if __name__ == "__main__":
    try:
        success = run_cavity_test()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
