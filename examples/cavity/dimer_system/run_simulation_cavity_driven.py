#!/usr/bin/env python3
"""
Example: Cavity-coupled dimer with cavity-mode driving (input-output theory)

Demonstrates:
1. Laser drives cavity mode (cavity acts as frequency filter)
2. Time-dependent Gaussian pulse
3. Cavity-mediated energy transfer to molecules

This script shows Case (2) from the theory: the external laser injects energy
into the cavity through boundaries (mirrors), and the cavity acts as a frequency
filter. The molecule experiences a field shaped by the cavity response rather
than the raw laser spectrum.
"""

import sys
import numpy as np
from pathlib import Path
import json
import os
import argparse
from datetime import datetime

try:
    # Find OpenMM module - check build directory first (where it was actually built)
    import sys
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
    import os
    project_root = Path(__file__).parent.parent.parent.parent
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

def run_cavity_driven_simulation(
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
    drive_amplitude=0.01,
    drive_freq_cm=None,
    drive_phase=0.0,
    envelope_type="gaussian",
    envelope_param1=50.0,
    envelope_param2=10.0,
    friction=0.01,
    minimize=True,
    output_file="cavity_driven_dipole_trajectory.npz"
):
    """Run simulation with cavity-mode driving (Case 2).
    
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
    drive_amplitude : float
        Driving force amplitude f₀ in atomic units
    drive_freq_cm : float, optional
        Drive frequency ω_d in cm⁻¹ (if None, uses cavity frequency)
    drive_phase : float
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
    print("Cavity-Driven Simulation (Input-Output Theory)")
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
    
    # Handle drive frequency default (use cavity frequency if not specified)
    if drive_freq_cm is None:
        drive_freq_cm = cavity_freq_cm
    omega_d_au = drive_freq_cm / 219474.63  # Convert to atomic units
    
    # Laser parameters (Case 2: cavity-mode driving)
    f0 = drive_amplitude
    omega_d = omega_d_au
    phase_d = drive_phase
    
    # Set envelope parameters based on type
    if envelope_type == "gaussian":
        peak_time_ps = envelope_param1
        width_ps = envelope_param2
    elif envelope_type == "square":
        peak_time_ps = envelope_param1  # start_time
        width_ps = envelope_param2    # stop_time
    elif envelope_type == "exponential":
        peak_time_ps = envelope_param1  # tau
        width_ps = envelope_param2      # unused
    else:  # constant
        peak_time_ps = 0.0
        width_ps = 0.0
    
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
    
    print(f"\nLaser parameters (Case 2: Cavity-mode driving):")
    print(f"  Drive amplitude f_0: {f0} a.u.")
    print(f"  Drive frequency ω_d: {drive_freq_cm} cm⁻¹ ({omega_d:.6f} a.u.)")
    print(f"  Phase: {phase_d}")
    if envelope_type == "gaussian":
        print(f"  Envelope: Gaussian (peak={peak_time_ps} ps, width={width_ps} ps)")
    elif envelope_type == "square":
        print(f"  Envelope: Square (start={peak_time_ps} ps, stop={width_ps} ps)")
    elif envelope_type == "exponential":
        print(f"  Envelope: Exponential (tau={peak_time_ps} ps)")
    else:
        print(f"  Envelope: Constant")
    
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
    print("\n--- Adding Cavity Force with Laser Driving ---")
    try:
        cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass)
    except Exception as e:
        raise
    
    # Configure cavity driving (Case 2)
    cavity_force.setCavityDriveAmplitude(f0)
    cavity_force.setCavityDriveFrequency(omega_d)
    cavity_force.setCavityDrivePhase(phase_d)
    cavity_force.setCavityDriveEnvelope(envelope_type, peak_time_ps, width_ps)
    # Disable laser forces initially if using equilibration mode
    if use_equilibration:
        cavity_force.setCavityDriveEnabled(False)
    else:
        cavity_force.setCavityDriveEnabled(True)
    cavity_force.setDirectLaserCouplingEnabled(False)  # Only cavity driving

    system.addForce(cavity_force)
    print("  Cavity-mode driving enabled")
    print("  Direct molecule-laser coupling disabled")
    
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

    try:
        context = openmm.Context(system, integrator, platform)
    except Exception as e:
        raise

    try:
        context.setPositions(positions)
    except Exception as e:
        raise

    # Set initial time to 0 (important for time-dependent forces)
    try:
        context.setTime(0.0 * unit.picosecond)
    except Exception as e:
        raise
    
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
    try:
        context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    except Exception as e:
        raise
    
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
        cavity_force.setCavityDriveEnabled(True)
        cavity_force.updateParametersInContext(context)
        print("  Laser forces enabled for production")
        
        # Reset time for production (optional - comment out if you want continuous time)
        # context.setTime(0.0 * unit.picosecond)
        
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
    drive_energies = []
    
    # Track dipole moments for IR spectrum (every 4 fs)
    dipole_times = []
    dipole_moments = []
    
    # Track cavity positions (q) for spectrum analysis (every 4 fs)
    cavity_positions = []
    
    # Track simulation speed (instantaneous)
    import time
    simulation_start_time = time.time()
    last_report_wall_time = simulation_start_time
    last_report_sim_time_ps = 0.0

    # Run simulation with fine-grained dipole sampling
    # Only save dipole during production phase
    for step in range(start_step, start_step + total_steps):
        try:
            integrator.step(1)
        except Exception as e:
            raise
        
        # Save dipole moment and cavity position every 4 fs (0.004 ps = 4 steps)
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
                
                # Get cavity position (q) - typically we use x,y components for cavity mode
                cavity_pos = positions[cavity_index]
                # Convert to numpy array in nm
                cavity_q = np.array([
                    cavity_pos[0].value_in_unit(unit.nanometer),
                    cavity_pos[1].value_in_unit(unit.nanometer),
                    cavity_pos[2].value_in_unit(unit.nanometer)
                ])
                cavity_positions.append(cavity_q)
            except Exception as e:
                print(f"  ⚠ Error computing dipole/cavity position at step {step}: {e}")
        
        # Report energies every report_interval (1 ps)
        if step % report_interval == 0 and step > 0:
            try:
                state = context.getState(getEnergy=True)
                time_ps = state.getTime().value_in_unit(unit.picosecond)
                
                # Get energy components
                harmonic = cavity_force.getHarmonicEnergy(context)
                coupling = cavity_force.getCouplingEnergy(context)
                drive = cavity_force.getCavityDriveEnergy(context)
                total_cavity = cavity_force.getTotalCavityEnergy(context)

                times.append(time_ps)
                cavity_energies.append(total_cavity)
                drive_energies.append(drive)
                
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
                drive_val = drive.value_in_unit(unit.kilojoule_per_mole) if hasattr(drive, 'value_in_unit') else drive
                print(f"  t={time_val:.1f} ps: E_cavity={total_cavity_val:.4f} kJ/mol, "
                      f"E_drive={drive_val:.4f} kJ/mol, Speed={instantaneous_speed_ns_per_day:.1f} ns/day")
            except Exception as e:
                print(f"  ⚠ Error retrieving energy at step {step}: {e}")
    
    # Calculate simulation speed
    simulation_end_time = time.time()
    wall_clock_time_seconds = simulation_end_time - simulation_start_time
    wall_clock_time_days = wall_clock_time_seconds / (24.0 * 3600.0)
    # Use production time for speed calculation
    sim_time_for_speed = prod_time if use_equilibration else total_time_ps
    simulation_time_ns = sim_time_for_speed / 1000.0  # Convert ps to ns
    speed_ns_per_day = simulation_time_ns / wall_clock_time_days if wall_clock_time_days > 0 else 0.0
    
    # Save dipole moment and cavity position trajectories for spectrum analysis
    dipole_times_array = np.array(dipole_times)
    dipole_moments_array = np.array(dipole_moments)
    cavity_positions_array = np.array(cavity_positions)
    np.savez(output_file, 
             time_ps=dipole_times_array,  # Match analyze_spectrum.py expected key
             dipole_nm=dipole_moments_array,  # Match analyze_spectrum.py expected key
             cavity_q_nm=cavity_positions_array)  # Cavity position q in nm
    print(f"\n  Saved dipole trajectory: {len(dipole_times)} points ({output_file})")
    print(f"  Saved cavity position trajectory: {len(cavity_positions)} points ({output_file})")
    
    print("\n--- Simulation Complete ---")
    if len(cavity_energies) > 0:
        final_cavity = cavity_energies[-1].value_in_unit(unit.kilojoule_per_mole) if hasattr(cavity_energies[-1], 'value_in_unit') else cavity_energies[-1]
        final_drive = drive_energies[-1].value_in_unit(unit.kilojoule_per_mole) if hasattr(drive_energies[-1], 'value_in_unit') else drive_energies[-1]
        print(f"  Final cavity energy: {final_cavity:.4f} kJ/mol")
        print(f"  Final drive energy: {final_drive:.4f} kJ/mol")
    else:
        # Get final energy if no reports were made (very short simulation)
        try:
            state = context.getState(getEnergy=True)
            harmonic = cavity_force.getHarmonicEnergy(context)
            coupling = cavity_force.getCouplingEnergy(context)
            drive = cavity_force.getCavityDriveEnergy(context)
            total_cavity = cavity_force.getTotalCavityEnergy(context)
            final_cavity = total_cavity.value_in_unit(unit.kilojoule_per_mole) if hasattr(total_cavity, 'value_in_unit') else total_cavity
            final_drive = drive.value_in_unit(unit.kilojoule_per_mole) if hasattr(drive, 'value_in_unit') else drive
            print(f"  Final cavity energy: {final_cavity:.4f} kJ/mol")
            print(f"  Final drive energy: {final_drive:.4f} kJ/mol")
        except Exception as e:
            print(f"  Could not retrieve final energies: {e}")
    print(f"\n  Simulation speed: {speed_ns_per_day:.2f} ns/day")
    print(f"  Wall-clock time: {wall_clock_time_seconds:.2f} seconds")
    print("\nNote: The cavity acts as a frequency filter - the molecule")
    print("experiences a field shaped by the cavity response at ω_c.")

def main():
    """Main entry point with command-line argument parsing."""
    parser = argparse.ArgumentParser(
        description='Cavity-Driven Dimer Simulation (Input-Output Theory)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default parameters
  python run_simulation_cavity_driven.py

  # Custom parameters
  python run_simulation_cavity_driven.py --dimers 100 --time 200 --cavity-freq 1560 --drive-freq 1600

  # With equilibration and production phases
  python run_simulation_cavity_driven.py --equil 20 --prod 80

  # Square pulse
  python run_simulation_cavity_driven.py --envelope-type square --envelope-param1 0 --envelope-param2 50
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
    
    # Laser parameters (cavity-mode driving)
    parser.add_argument('--drive-amplitude', type=float, default=0.01,
                       help='Driving force amplitude f₀ in atomic units (default: 0.01)')
    parser.add_argument('--drive-freq', type=float, default=None,
                       help='Drive frequency ω_d in cm⁻¹ (default: uses cavity frequency)')
    parser.add_argument('--drive-phase', type=float, default=0.0,
                       help='Phase φ in radians (default: 0.0)')
    parser.add_argument('--envelope-type', type=str, default='gaussian',
                       choices=['constant', 'gaussian', 'square', 'exponential'],
                       help='Envelope type: constant, gaussian, square, or exponential (default: gaussian)')
    parser.add_argument('--envelope-param1', type=float, default=50.0,
                       help='First envelope parameter: peak_time (gaussian), start_time (square), or tau (exponential) (default: 50.0)')
    parser.add_argument('--envelope-param2', type=float, default=10.0,
                       help='Second envelope parameter: width (gaussian), stop_time (square), or unused (exponential) (default: 10.0)')
    
    # Simulation control
    parser.add_argument('--friction', type=float, default=0.01,
                       help='Friction coefficient in ps⁻¹ (default: 0.01)')
    parser.add_argument('--no-minimize', action='store_true',
                       help='Skip energy minimization')
    parser.add_argument('--output', type=str, default='cavity_driven_dipole_trajectory.npz',
                       help='Output filename for dipole trajectory (default: cavity_driven_dipole_trajectory.npz)')
    
    args = parser.parse_args()
    
    try:
        run_cavity_driven_simulation(
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
            drive_amplitude=args.drive_amplitude,
            drive_freq_cm=args.drive_freq,
            drive_phase=args.drive_phase,
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
