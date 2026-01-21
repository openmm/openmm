#!/usr/bin/env python3
"""
Baseline performance test: 50 molecules, 1 fs timestep, NO cavity
This matches the parameters of the running IR simulation to compare speed
"""
import sys
import numpy as np
from openmm import openmm, unit
import time

# Physical constants
BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5

def create_diamer_system(num_molecules=50, fraction_OO=0.8, box_size_nm=2.5, 
                         temperature_K=100.0, seed=42):
    """Create a system of O-O and N-N dimers."""
    np.random.seed(seed)
    
    system = openmm.System()
    positions = []
    
    # Particle masses
    mass_O = 16.0  # amu
    mass_N = 14.0  # amu
    
    # Bond parameters (doubled for OpenMM convention)
    k_OO_au = 2 * 0.73204
    r0_OO_au = 2.281655158
    k_OO = k_OO_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
    r0_OO = r0_OO_au * BOHR_TO_NM
    
    k_NN_au = 2 * 1.4325
    r0_NN_au = 2.067661567
    k_NN = k_NN_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
    r0_NN = r0_NN_au * BOHR_TO_NM
    
    # Charges (from cav-hoomd)
    charge_O = -0.3  # elementary charge
    charge_N = -0.2
    
    # Create molecules
    num_OO = int(num_molecules * fraction_OO)
    num_NN = num_molecules - num_OO
    
    harmonic_force = openmm.HarmonicBondForce()
    nonbonded_force = openmm.NonbondedForce()
    nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nonbonded_force.setCutoffDistance(1.0 * unit.nanometer)
    
    particle_index = 0
    
    # Add O-O dimers
    for i in range(num_OO):
        x, y, z = np.random.uniform(-box_size_nm/2, box_size_nm/2, 3)
        
        system.addParticle(mass_O)
        system.addParticle(mass_O)
        positions.append(openmm.Vec3(x, y, z) * unit.nanometer)
        positions.append(openmm.Vec3(x + r0_OO, y, z) * unit.nanometer)
        
        harmonic_force.addBond(particle_index, particle_index + 1, r0_OO, k_OO)
        nonbonded_force.addParticle(charge_O, 0.1, 0.0)
        nonbonded_force.addParticle(charge_O, 0.1, 0.0)
        nonbonded_force.addException(particle_index, particle_index + 1, 0.0, 1.0, 0.0)
        
        particle_index += 2
    
    # Add N-N dimers
    for i in range(num_NN):
        x, y, z = np.random.uniform(-box_size_nm/2, box_size_nm/2, 3)
        
        system.addParticle(mass_N)
        system.addParticle(mass_N)
        positions.append(openmm.Vec3(x, y, z) * unit.nanometer)
        positions.append(openmm.Vec3(x + r0_NN, y, z) * unit.nanometer)
        
        harmonic_force.addBond(particle_index, particle_index + 1, r0_NN, k_NN)
        nonbonded_force.addParticle(charge_N, 0.1, 0.0)
        nonbonded_force.addParticle(charge_N, 0.1, 0.0)
        nonbonded_force.addException(particle_index, particle_index + 1, 0.0, 1.0, 0.0)
        
        particle_index += 2
    
    system.addForce(harmonic_force)
    system.addForce(nonbonded_force)
    
    # Periodic box
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_size_nm, 0, 0),
        openmm.Vec3(0, box_size_nm, 0),
        openmm.Vec3(0, 0, box_size_nm)
    )
    
    print(f"Created system with {system.getNumParticles()} particles ({num_molecules} molecules)")
    print(f"  O-O dimers: {num_OO}")
    print(f"  N-N dimers: {num_NN}")
    print(f"  Box size: {box_size_nm} nm")
    
    return system, positions

def run_baseline_test():
    """Run baseline test matching IR simulation parameters."""
    print("=" * 70)
    print("Baseline Performance Test (NO Cavity)")
    print("Matching IR simulation: 50 mol, 1 fs timestep")
    print("=" * 70)
    
    # Parameters matching the IR simulation
    num_molecules = 50
    box_size_nm = 2.5
    temperature_K = 100.0
    dt = 0.001  # 1 fs timestep
    
    # Short test
    test_time_ps = 50.0
    test_steps = int(test_time_ps / dt)
    
    print(f"\nSimulation parameters:")
    print(f"  Number of molecules: {num_molecules} ({num_molecules * 2} atoms)")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt} ps (1 fs)")
    print(f"  Test run: {test_time_ps} ps ({test_steps} steps)")
    
    # Create system
    print("\n--- Creating Diamer System ---")
    system, positions = create_diamer_system(
        num_molecules=num_molecules,
        fraction_OO=0.8,
        box_size_nm=box_size_nm,
        temperature_K=temperature_K
    )
    
    # Add Bussi thermostat (matching IR simulation)
    print("\n--- Adding Bussi Thermostat ---")
    try:
        bussi = openmm.BussiThermostat(
            temperature_K * unit.kelvin,
            0.1 * unit.picoseconds,
            system.getNumParticles()
        )
        system.addForce(bussi)
        print(f"  BussiThermostat added for {system.getNumParticles()} particles")
        print(f"  Temperature: {temperature_K} K, Tau: 0.1 ps")
    except:
        print("  BussiThermostat not available, using Langevin only")
    
    # Create integrator (matching IR simulation)
    print("\n--- Creating Integrator ---")
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        0.01 / unit.picoseconds,  # Very low friction (like IR sim)
        dt * unit.picoseconds
    )
    print(f"  LangevinMiddleIntegrator created (friction: 0.01 ps^-1)")
    
    # Create simulation on CUDA
    print("\n--- Creating Simulation ---")
    platform = openmm.Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed'}
    context = openmm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    print(f"  Using CUDA platform")
    print(f"  Precision: mixed")
    
    # Energy minimization
    print("\n--- Energy Minimization ---")
    initial_state = context.getState(getEnergy=True)
    print(f"  Initial energy: {initial_state.getPotentialEnergy()}")
    openmm.LocalEnergyMinimizer.minimize(context, 1.0, 100)
    final_state = context.getState(getEnergy=True)
    print(f"  Final energy: {final_state.getPotentialEnergy()}")
    
    # Equilibration (10 ps)
    print("\n--- Equilibration Phase ---")
    equil_steps = int(10.0 / dt)
    print(f"  Running {equil_steps} steps (10.0 ps)...")
    start = time.time()
    integrator.step(equil_steps)
    equil_time = time.time() - start
    print(f"  Completed in {equil_time:.2f} seconds")
    state = context.getState(getEnergy=True, getPositions=True)
    print(f"  Final PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.2f} kJ/mol")
    
    # Performance test
    print("\n--- Performance Test Run ---")
    print(f"  Running {test_steps} steps ({test_time_ps} ps)...")
    
    report_interval = 10000  # Every 10 ps
    start_time = time.time()
    
    for i in range(0, test_steps, report_interval):
        integrator.step(report_interval)
        state = context.getState(getEnergy=True)
        elapsed = time.time() - start_time
        steps_done = i + report_interval
        progress = 100.0 * steps_done / test_steps
        
        # Calculate performance
        steps_per_sec = steps_done / elapsed
        ns_per_day = (steps_per_sec * dt * 86400)  # dt in ps, convert to ns
        us_per_day = ns_per_day / 1000.0
        
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        print(f"  [{progress:5.1f}%] Step {steps_done:7d}/{test_steps} | PE: {pe:8.2f} kJ/mol | Speed: {steps_per_sec:6.1f} steps/s")
        print(f"         Perf: {ns_per_day:6.1f} ns/day ({us_per_day:5.2f} μs/day)")
    
    # Final summary
    total_time = time.time() - start_time
    avg_steps_per_sec = test_steps / total_time
    avg_ns_per_day = avg_steps_per_sec * dt * 86400
    avg_us_per_day = avg_ns_per_day / 1000.0
    
    print("\n" + "=" * 70)
    print("PERFORMANCE SUMMARY")
    print("=" * 70)
    print(f"Platform: CUDA")
    print(f"System size: {num_molecules} molecules ({num_molecules * 2} atoms)")
    print(f"Timestep: {dt} ps (1 fs)")
    print(f"Test duration: {test_time_ps} ps ({test_steps} steps)")
    print(f"Wall time: {total_time:.2f} seconds")
    print()
    print(f"Speed: {avg_steps_per_sec:.1f} steps/second")
    print(f"Performance: {avg_ns_per_day:.1f} ns/day")
    print(f"             {avg_us_per_day:.2f} μs/day")
    print("=" * 70)
    
    return True

if __name__ == '__main__':
    try:
        import openmm
        print("✓ OpenMM loaded successfully")
    except ImportError as e:
        print(f"✗ Failed to import OpenMM: {e}")
        sys.exit(1)
    
    run_baseline_test()
