#!/usr/bin/env python3
"""
Example: Cavity-coupled dimer with direct laser-molecule coupling

Demonstrates:
1. Laser directly drives molecular dipoles (bypasses cavity filtering)
2. Square pulse envelope
3. Direct energy injection into molecular vibrations

This script shows Case (1) from the theory: the external laser field couples
directly to the molecular dipole, bypassing the cavity. The molecule sees the
raw laser spectrum at ω_L without frequency filtering.
"""

import sys
import numpy as np
import argparse
from pathlib import Path

try:
    # Find OpenMM module - check build directory first (where it was actually built)
    import os
    from pathlib import Path
    
    # Get the project root (3 levels up from this script)
    project_root = Path(__file__).parent.parent.parent.parent
    build_lib = project_root / "build" / "python" / "build" / "lib.linux-x86_64-cpython-312"
    
    # Remove problematic build/python paths that have incomplete modules
    sys.path = [p for p in sys.path if not (('build/python' in p or 'build/lib' in p) and str(project_root) not in p)]
    
    # Add the correct build directory
    if str(build_lib) not in sys.path:
        sys.path.insert(0, str(build_lib))
    
    from openmm import openmm
    from openmm import unit
    
    # Load CUDA platform plugin if available
    cuda_lib = project_root / "build" / "libOpenMMCUDA.so"
    if cuda_lib.exists():
        try:
            openmm.Platform.loadPluginLibrary(str(cuda_lib))
            print("CUDA platform plugin loaded")
        except Exception as e:
            print(f"⚠ Could not load CUDA plugin: {e}")
    
    print("OpenMM loaded successfully")
except ImportError as e:
    print(f"Error importing OpenMM: {e}")
    sys.exit(1)

# Import helper functions from the main simulation script
sys.path.insert(0, str(Path(__file__).parent))
from run_simulation import create_diamer_system_in_code, add_cavity_particle, compute_dipole_moment

# Physical constants
BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5

def run_direct_laser_simulation(
    num_molecules=50,
    box_size_nm=2.5,
    temperature_K=100.0,
    dt=0.001,
    total_time_ps=100.0,
    equilibration_time_ps=None,
    production_time_ps=None,
    fraction_OO=0.8,
    cavity_freq_cm=1560.0,
    lambda_coupling=0.01,
    laser_amplitude=0.005,
    laser_freq_cm=2000.0,
    laser_phase=0.0,
    envelope_type="square",
    envelope_param1=0.0,
    envelope_param2=50.0,
    friction=0.01,
    minimize=True,
    output_file="direct_laser_dipole_trajectory.npz"
):
    """Run simulation with direct laser-molecule coupling (Case 1).
    
    Parameters
    ----------
    num_molecules : int
        Number of dimers
    box_size_nm : float
        Box size in nanometers
    temperature_K : float
        Temperature in Kelvin
    dt : float
        Timestep in picoseconds
    total_time_ps : float
        Total simulation time in picoseconds (used if equilibration/production not specified)
    equilibration_time_ps : float, optional
        Equilibration time in picoseconds (laser forces disabled during equilibration)
    production_time_ps : float, optional
        Production time in picoseconds (laser forces enabled during production)
    fraction_OO : float
        Fraction of O-O dimers (vs N-N)
    cavity_freq_cm : float
        Cavity frequency in cm⁻¹
    lambda_coupling : float
        Coupling strength λ in atomic units
    laser_amplitude : float
        Electric field amplitude E₀ in atomic units
    laser_freq_cm : float
        Laser frequency ω_L in cm⁻¹
    laser_phase : float
        Phase φ in radians
    envelope_type : str
        Envelope type: "constant", "gaussian", "square", or "exponential"
    envelope_param1 : float
        First envelope parameter (peak_time for gaussian, start_time for square, tau for exponential)
    envelope_param2 : float
        Second envelope parameter (width for gaussian, stop_time for square, unused for exponential)
    friction : float
        Friction coefficient in ps⁻¹
    minimize : bool
        Whether to perform energy minimization
    output_file : str
        Output filename for dipole trajectory
    """
    
    print("=" * 60)
    print("Direct Laser-Molecule Coupling Simulation")
    print("=" * 60)
    
    # Validate envelope type
    valid_envelopes = ["constant", "gaussian", "square", "exponential"]
    if envelope_type not in valid_envelopes:
        raise ValueError(f"Invalid envelope type: {envelope_type}. Must be one of {valid_envelopes}")
    
    # Handle equilibration/production split
    if equilibration_time_ps is not None and production_time_ps is not None:
        # Use equilibration + production mode
        equil_time = equilibration_time_ps
        prod_time = production_time_ps
        total_time_ps = equil_time + prod_time
        use_equilibration = True
    else:
        # Use single simulation mode (backward compatible)
        equil_time = 0.0
        prod_time = total_time_ps
        use_equilibration = False
    
    # Convert cavity frequency to atomic units
    omegac_au = cavity_freq_cm / 219474.63  # Hartree
    photon_mass = 1.0 / 1822.888  # amu
    
    # Convert laser frequency to atomic units
    omega_L_au = laser_freq_cm / 219474.63  # Hartree
    
    # Laser parameters (Case 1: direct molecule-laser coupling)
    E0 = laser_amplitude
    omega_L = omega_L_au
    phase_L = laser_phase
    
    # Set envelope parameters based on type
    if envelope_type == "square":
        start_time_ps = envelope_param1
        stop_time_ps = envelope_param2
    elif envelope_type == "gaussian":
        start_time_ps = envelope_param1  # peak_time
        stop_time_ps = envelope_param2   # width
    elif envelope_type == "exponential":
        start_time_ps = envelope_param1  # tau
        stop_time_ps = envelope_param2   # unused
    else:  # constant
        start_time_ps = 0.0
        stop_time_ps = 0.0  # Will be ignored for constant
    
    print(f"\nSimulation parameters:")
    print(f"  Number of dimers: {num_molecules}")
    print(f"  Box size: {box_size_nm} nm")
    print(f"  Fraction O-O: {fraction_OO}")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt} ps")
    if use_equilibration:
        print(f"  Equilibration time: {equil_time} ps (laser forces disabled)")
        print(f"  Production time: {prod_time} ps (laser forces enabled)")
        print(f"  Total time: {total_time_ps} ps")
    else:
        print(f"  Total time: {total_time_ps} ps")
    print(f"  Cavity frequency: {cavity_freq_cm} cm⁻¹ ({omegac_au:.6f} Hartree)")
    print(f"  Lambda coupling: {lambda_coupling}")
    
    print(f"\nLaser parameters (Case 1: Direct molecule-laser coupling):")
    print(f"  Field amplitude E_0: {E0} a.u.")
    print(f"  Laser frequency ω_L: {laser_freq_cm} cm⁻¹ ({omega_L:.6f} a.u.)")
    print(f"  Phase: {phase_L}")
    if envelope_type == "gaussian":
        print(f"  Envelope: Gaussian (peak={start_time_ps} ps, width={stop_time_ps} ps)")
    elif envelope_type == "square":
        print(f"  Envelope: Square (start={start_time_ps} ps, stop={stop_time_ps} ps)")
    elif envelope_type == "exponential":
        print(f"  Envelope: Exponential (tau={start_time_ps} ps)")
    else:
        print(f"  Envelope: Constant (always on)")
    
    # Create system
    print("\n--- Creating System ---")
    system, positions = create_diamer_system_in_code(
        num_molecules=num_molecules,
        fraction_OO=fraction_OO,
        box_size_nm=box_size_nm,
        temperature_K=temperature_K
    )
    
    # Store charges
    charges = []
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for i in range(system.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(i)
                charges.append(charge)
    
    # Add cavity particle
    cavity_index = add_cavity_particle(system, positions, omegac_au, photon_mass)
    num_molecular_particles = cavity_index
    
    # Add CavityForce
    print("\n--- Adding Cavity Force with Direct Laser Coupling ---")
    cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass)
    
    # Configure direct laser coupling (Case 1)
    cavity_force.setDirectLaserAmplitude(E0)
    cavity_force.setDirectLaserFrequency(omega_L)
    cavity_force.setDirectLaserPhase(phase_L)
    cavity_force.setDirectLaserEnvelope(envelope_type, start_time_ps, stop_time_ps)
    # Disable laser forces initially if using equilibration mode
    if use_equilibration:
        cavity_force.setDirectLaserCouplingEnabled(False)
    else:
        cavity_force.setDirectLaserCouplingEnabled(True)
    cavity_force.setCavityDriveEnabled(False)  # No cavity driving
    
    system.addForce(cavity_force)
    print("  Direct molecule-laser coupling enabled")
    print("  Cavity-mode driving disabled")
    
    # Create integrator
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt * unit.picosecond
    )
    
    # Create simulation
    # Check available platforms
    num_platforms = openmm.Platform.getNumPlatforms()
    platform_names = [openmm.Platform.getPlatform(i).getName() for i in range(num_platforms)]
    print(f"\n  Available platforms: {platform_names}")
    
    # Try CUDA first, fall back to Reference if not available
    platform = None
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        print("  Using CUDA platform")
    except Exception as e:
        print(f"  ⚠ CUDA not available: {e}")
        platform = openmm.Platform.getPlatformByName('Reference')
        print("  Using Reference platform")
    
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)
    
    # Set initial time to 0 (important for time-dependent forces)
    context.setTime(0.0 * unit.picosecond)
    
    # Minimize with laser forces temporarily disabled
    if minimize:
        print("\n--- Energy Minimization ---")
        
        # Store current laser enabled state
        cavity_drive_was_enabled = cavity_force.getCavityDriveEnabled()
        direct_laser_was_enabled = cavity_force.getDirectLaserCouplingEnabled()
        
        # Temporarily disable laser forces for minimization
        cavity_force.setCavityDriveEnabled(False)
        cavity_force.setDirectLaserCouplingEnabled(False)
        # Update context to apply the changes
        cavity_force.updateParametersInContext(context)
        print("  Temporarily disabled laser forces for minimization...")
        
        try:
            openmm.LocalEnergyMinimizer.minimize(context, maxIterations=100)
            print("  Minimization completed")
        except Exception as e:
            print(f"  ⚠ Minimization failed: {e}")
            print("  Continuing with initial positions...")
        
        # Re-enable laser forces for simulation
        cavity_force.setCavityDriveEnabled(cavity_drive_was_enabled)
        cavity_force.setDirectLaserCouplingEnabled(direct_laser_was_enabled)
        cavity_force.updateParametersInContext(context)
        if cavity_drive_was_enabled or direct_laser_was_enabled:
            print("  Laser forces re-enabled for simulation")
    else:
        print("\n--- Skipping Energy Minimization ---")
    
    # Set velocities
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    
    # Run simulation
    if use_equilibration:
        print(f"\n--- Running Equilibration ({equil_time} ps) ---")
        equil_steps = int(equil_time / dt)
        for step in range(equil_steps):
            integrator.step(1)
            if step % 1000 == 0 and step > 0:
                time_ps = context.getState().getTime().value_in_unit(unit.picosecond)
                print(f"  Equilibration: t={time_ps:.1f} ps")
        
        # Enable laser forces for production
        print(f"\n--- Starting Production ({prod_time} ps) ---")
        cavity_force.setDirectLaserCouplingEnabled(True)
        cavity_force.updateParametersInContext(context)
        print("  Laser forces enabled for production")
        
        total_steps = int(prod_time / dt)
        start_step = 0
    else:
        print(f"\n--- Running Simulation ({total_time_ps} ps) ---")
        total_steps = int(total_time_ps / dt)
        start_step = 0
    
    report_interval = 1000  # Every 1 ps for console output
    dipole_save_interval = 4  # Every 4 fs (0.004 ps) for IR spectrum analysis
    
    # Track energies
    times = []
    cavity_energies = []
    laser_energies = []
    
    # Track dipole moments for IR spectrum (every 4 fs)
    dipole_times = []
    dipole_moments = []
    
    # Track simulation speed (instantaneous)
    import time
    simulation_start_time = time.time()
    last_report_wall_time = simulation_start_time
    last_report_sim_time_ps = 0.0
    
    # Run simulation with fine-grained dipole sampling
    # Only save dipole during production phase
    for step in range(start_step, start_step + total_steps):
        integrator.step(1)  # Step one timestep at a time
        
        # Save dipole moment every 4 fs (0.004 ps = 4 steps)
        if step % dipole_save_interval == 0:
            try:
                state = context.getState(getPositions=True)
                positions = state.getPositions()
                time_ps = state.getTime().value_in_unit(unit.picosecond)
                
                # Compute dipole moment (excluding cavity particle)
                # compute_dipole_moment expects (state, charges, num_molecular_particles)
                dipole = compute_dipole_moment(state, charges, num_molecular_particles)
                dipole_times.append(time_ps)
                dipole_moments.append(dipole)
            except Exception as e:
                print(f"  ⚠ Error computing dipole at step {step}: {e}")
        
        # Report energies every report_interval (1 ps)
        if step % report_interval == 0 and step > 0:
            try:
                state = context.getState(getEnergy=True)
                time_ps = state.getTime().value_in_unit(unit.picosecond)
                
                # Get energy components
                harmonic = cavity_force.getHarmonicEnergy(context)
                coupling = cavity_force.getCouplingEnergy(context)
                laser = cavity_force.getDirectLaserEnergy(context)
                total_cavity = cavity_force.getTotalCavityEnergy(context)
                
                times.append(time_ps)
                cavity_energies.append(total_cavity)
                laser_energies.append(laser)
                
                # Calculate instantaneous simulation speed (since last report)
                current_wall_time = time.time()
                time_val = time_ps.value_in_unit(unit.picosecond) if hasattr(time_ps, 'value_in_unit') else time_ps
                
                # Instantaneous speed: (sim_time_interval) / (wall_time_interval)
                wall_time_interval_seconds = current_wall_time - last_report_wall_time
                wall_time_interval_days = wall_time_interval_seconds / (24.0 * 3600.0)
                sim_time_interval_ps = time_val - last_report_sim_time_ps
                sim_time_interval_ns = sim_time_interval_ps / 1000.0  # Convert ps to ns
                
                if wall_time_interval_days > 0:
                    instantaneous_speed_ns_per_day = sim_time_interval_ns / wall_time_interval_days
                else:
                    instantaneous_speed_ns_per_day = 0.0
                
                # Update for next interval
                last_report_wall_time = current_wall_time
                last_report_sim_time_ps = time_val
                
                total_cavity_val = total_cavity.value_in_unit(unit.kilojoule_per_mole) if hasattr(total_cavity, 'value_in_unit') else total_cavity
                laser_val = laser.value_in_unit(unit.kilojoule_per_mole) if hasattr(laser, 'value_in_unit') else laser
                print(f"  t={time_val:.1f} ps: E_cavity={total_cavity_val:.4f} kJ/mol, "
                      f"E_laser={laser_val:.4f} kJ/mol, Speed={instantaneous_speed_ns_per_day:.1f} ns/day")
            except Exception as e:
                print(f"  ⚠ Error retrieving energy at step {step}: {e}")
    
    # Save dipole moment trajectory for IR spectrum analysis
    dipole_times_array = np.array(dipole_times)
    dipole_moments_array = np.array(dipole_moments)
    np.savez(output_file, 
             time_ps=dipole_times_array,  # Match analyze_spectrum.py expected key
             dipole_nm=dipole_moments_array)  # Match analyze_spectrum.py expected key
    print(f"\n  Saved dipole trajectory: {len(dipole_times)} points ({output_file})")
    
    # Calculate simulation speed
    simulation_end_time = time.time()
    wall_clock_time_seconds = simulation_end_time - simulation_start_time
    wall_clock_time_days = wall_clock_time_seconds / (24.0 * 3600.0)
    # Use production time for speed calculation
    sim_time_for_speed = prod_time if use_equilibration else total_time_ps
    simulation_time_ns = sim_time_for_speed / 1000.0  # Convert ps to ns
    speed_ns_per_day = simulation_time_ns / wall_clock_time_days if wall_clock_time_days > 0 else 0.0
    
    print("\n--- Simulation Complete ---")
    if len(cavity_energies) > 0:
        final_cavity = cavity_energies[-1].value_in_unit(unit.kilojoule_per_mole) if hasattr(cavity_energies[-1], 'value_in_unit') else cavity_energies[-1]
        final_laser = laser_energies[-1].value_in_unit(unit.kilojoule_per_mole) if hasattr(laser_energies[-1], 'value_in_unit') else laser_energies[-1]
        print(f"  Final cavity energy: {final_cavity:.4f} kJ/mol")
        print(f"  Final laser energy: {final_laser:.4f} kJ/mol")
    else:
        # Get final energy if no reports were made (very short simulation)
        try:
            state = context.getState(getEnergy=True)
            harmonic = cavity_force.getHarmonicEnergy(context)
            coupling = cavity_force.getCouplingEnergy(context)
            laser = cavity_force.getDirectLaserEnergy(context)
            total_cavity = cavity_force.getTotalCavityEnergy(context)
            final_cavity = total_cavity.value_in_unit(unit.kilojoule_per_mole) if hasattr(total_cavity, 'value_in_unit') else total_cavity
            final_laser = laser.value_in_unit(unit.kilojoule_per_mole) if hasattr(laser, 'value_in_unit') else laser
            print(f"  Final cavity energy: {final_cavity:.4f} kJ/mol")
            print(f"  Final laser energy: {final_laser:.4f} kJ/mol")
        except Exception as e:
            print(f"  Could not retrieve final energies: {e}")
    print(f"\n  Simulation speed: {speed_ns_per_day:.2f} ns/day")
    print(f"  Wall-clock time: {wall_clock_time_seconds:.2f} seconds")
    print("\nNote: The laser directly forces molecular dipoles at ω_L.")
    print("The cavity does NOT filter the laser spectrum - any frequency")
    print("component in E_ext(t) acts directly on nuclei.")

def main():
    """Main entry point with command-line argument parsing."""
    parser = argparse.ArgumentParser(
        description='Direct Laser-Molecule Coupling Simulation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default parameters
  python run_simulation_direct_laser.py

  # Custom parameters with constant envelope (always on)
  python run_simulation_direct_laser.py --laser-freq 1555 --laser-amplitude 1.0 --envelope-type constant

  # With equilibration and production phases
  python run_simulation_direct_laser.py --equil 20 --prod 80 --envelope-type constant
        """
    )
    
    # System parameters
    parser.add_argument('--dimers', type=int, default=50,
                       help='Number of dimers (default: 50)')
    parser.add_argument('--box-size', type=float, default=2.5,
                       help='Box size in nm (default: 2.5)')
    parser.add_argument('--temp', type=float, default=100.0,
                       help='Temperature in K (default: 100.0)')
    parser.add_argument('--dt', type=float, default=0.001,
                       help='Timestep in ps (default: 0.001 = 1 fs)')
    parser.add_argument('--time', type=float, default=100.0,
                       help='Total simulation time in ps (default: 100.0, ignored if --equil and --prod are specified)')
    parser.add_argument('--equil', type=float, default=None,
                       help='Equilibration time in ps (laser forces disabled, default: None, uses --time if not specified)')
    parser.add_argument('--prod', type=float, default=None,
                       help='Production time in ps (laser forces enabled, default: None, uses --time if not specified)')
    parser.add_argument('--fraction-OO', type=float, default=0.8,
                       help='Fraction of O-O dimers (default: 0.8)')
    
    # Cavity parameters
    parser.add_argument('--cavity-freq', type=float, default=1560.0,
                       help='Cavity frequency in cm⁻¹ (default: 1560.0)')
    parser.add_argument('--lambda', type=float, default=0.01, dest='lambda_coupling',
                       help='Coupling strength λ in atomic units (default: 0.01)')
    
    # Laser parameters (direct molecule-laser coupling)
    parser.add_argument('--laser-amplitude', '--drive-amplitude', type=float, default=0.005,
                       help='Electric field amplitude E₀ in atomic units (default: 0.005)')
    parser.add_argument('--laser-freq', '--drive-freq', type=float, default=2000.0,
                       help='Laser frequency ω_L in cm⁻¹ (default: 2000.0)')
    parser.add_argument('--laser-phase', type=float, default=0.0,
                       help='Phase φ in radians (default: 0.0)')
    parser.add_argument('--envelope-type', type=str, default='square',
                       choices=['constant', 'gaussian', 'square', 'exponential'],
                       help='Envelope type: constant, gaussian, square, or exponential (default: square)')
    parser.add_argument('--envelope-param1', type=float, default=0.0,
                       help='First envelope parameter: peak_time (gaussian), start_time (square), or tau (exponential) (default: 0.0)')
    parser.add_argument('--envelope-param2', type=float, default=50.0,
                       help='Second envelope parameter: width (gaussian), stop_time (square), or unused (exponential) (default: 50.0)')
    
    # Simulation control
    parser.add_argument('--friction', type=float, default=0.01,
                       help='Friction coefficient in ps⁻¹ (default: 0.01)')
    parser.add_argument('--no-minimize', action='store_true',
                       help='Skip energy minimization')
    parser.add_argument('--output', type=str, default='direct_laser_dipole_trajectory.npz',
                       help='Output filename for dipole trajectory (default: direct_laser_dipole_trajectory.npz)')
    
    args = parser.parse_args()
    
    try:
        run_direct_laser_simulation(
            num_molecules=args.dimers,
            box_size_nm=args.box_size,
            temperature_K=args.temp,
            dt=args.dt,
            total_time_ps=args.time,
            equilibration_time_ps=args.equil,
            production_time_ps=args.prod,
            fraction_OO=args.fraction_OO,
            cavity_freq_cm=args.cavity_freq,
            lambda_coupling=args.lambda_coupling,
            laser_amplitude=args.laser_amplitude,
            laser_freq_cm=args.laser_freq,
            laser_phase=args.laser_phase,
            envelope_type=args.envelope_type,
            envelope_param1=args.envelope_param1,
            envelope_param2=args.envelope_param2,
            friction=args.friction,
            minimize=not args.no_minimize,
            output_file=args.output
        )
        sys.exit(0)
    except Exception as e:
        print(f"\nSimulation failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
