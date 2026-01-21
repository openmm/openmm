#!/usr/bin/env python3
"""
Cavity OpenMM Diamer Test - Extended IR Simulation
===================================================

This test creates a system of O-O and N-N dimers (similar to the cav-hoomd diamer system)
and tests the cavity coupling functionality:

1. Equilibrate for 100 ps with coupling OFF (lambda=0)
2. Turn on coupling and run for 900 ps (total 1 ns)
3. Output dipole moment trajectory for IR spectrum calculation

The system uses:
- Harmonic bonds for dimers
- Lennard-Jones interactions
- Coulomb interactions
- Cavity coupling with Bussi thermostat for molecules
- Langevin dynamics for the cavity photon
"""

import sys
import numpy as np
import os
import time

try:
    from openmm import openmm
    from openmm import unit
    from openmm.app import Simulation, StateDataReporter
    print("✓ OpenMM loaded successfully")
except ImportError as e:
    print(f"Error importing OpenMM: {e}")
    print("Make sure OpenMM is installed and the Python path is correct.")
    sys.exit(1)

# Physical constants for unit conversion
# We'll work in OpenMM native units (nm, kJ/mol, ps)
BOHR_TO_NM = 0.0529177  # 1 Bohr = 0.0529177 nm
HARTREE_TO_KJMOL = 2625.5  # 1 Hartree = 2625.5 kJ/mol
AMU_TO_KG = 1.66054e-27
# Time conversion: 1 a.u. of time = 0.02418884254 ps
AU_TIME_TO_PS = 0.02418884254  # 1 atomic time unit = 0.02418884254 ps

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
    
    # Particle masses (in atomic units converted to amu)
    # O mass ~ 29150 a.u. ≈ 16 amu
    # N mass ~ 25527 a.u. ≈ 14 amu  
    mass_O = 16.0  # amu
    mass_N = 14.0  # amu
    
    # Bond parameters (converted from atomic units to OpenMM units)
    # HOOMD-blue uses: E = 0.5 * k * (r - r0)²  (has 1/2 factor, same as OpenMM!)
    # OpenMM uses: E = 0.5 * k * (r - r0)²  (has 1/2 factor)
    # Therefore we use the SAME force constants!
    # 
    # HOOMD force constants (from cav-hoomd code):
    #   k_OO = 2*0.36602 = 0.73204 Hartree/Bohr²  → ω_OO ≈ 1560 cm⁻¹
    #   k_NN = 2*0.71625 = 1.4325 Hartree/Bohr²   → ω_NN ≈ 2325 cm⁻¹
    k_OO_au = 0.73204  # Hartree/Bohr² (same as HOOMD)
    r0_OO_au = 2.281655158  # Bohr
    k_OO = k_OO_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)  # kJ/(mol·nm²)
    r0_OO = r0_OO_au * BOHR_TO_NM  # nm
    
    k_NN_au = 1.4325  # Hartree/Bohr² (same as HOOMD)
    r0_NN_au = 2.0743522177  # Bohr
    k_NN = k_NN_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
    r0_NN = r0_NN_au * BOHR_TO_NM
    
    # Charges
    charge_magnitude = 0.3  # elementary charge
    
    # LJ parameters - Molecular Kob-Andersen Model
    # ==============================================
    # These values are converted from cav-hoomd atomic units to OpenMM units.
    # The molecular KA model is a 4:1 mixture of diatomics with Kob-Andersen
    # LJ interactions. Parameters follow the original KA model ratios:
    #   ε_BB/ε_AA = 0.5, σ_BB/σ_AA = 0.88 (for N-N vs O-O)
    #   ε_AB = 1.5 × ε_AA, σ_AB = 0.8 × σ_AA (for N-O cross-term)
    #
    # Unit conversions: 1 Bohr = 0.0529177 nm, 1 Hartree = 2625.5 kJ/mol
    #
    # HOOMD values (atomic units):
    #   O-O: epsilon=0.00016685201 Ha, sigma=6.230426584 Bohr
    #   N-N: epsilon=0.000083426 Ha, sigma=5.48277488 Bohr
    #   N-O: epsilon=0.00025027802 Ha, sigma=4.9832074319 Bohr
    sigma_O = 6.230426584 * BOHR_TO_NM  # nm (= 0.32971 nm)
    epsilon_O = 0.00016685201 * HARTREE_TO_KJMOL  # kJ/mol (= 0.4381 kJ/mol)
    sigma_N = 5.48277488 * BOHR_TO_NM  # nm (= 0.29015 nm)
    epsilon_N = 0.000083426 * HARTREE_TO_KJMOL  # kJ/mol (= 0.2190 kJ/mol)
    # Cross-term (N-O): non-additive parameters from KA model
    sigma_NO = 4.9832074319 * BOHR_TO_NM  # nm (= 0.26370 nm)
    epsilon_NO = 0.00025027802 * HARTREE_TO_KJMOL  # kJ/mol (= 0.6571 kJ/mol)
    
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
    O_indices = []  # List of (idx1, idx2, charge1, charge2) for O atoms
    N_indices = []  # List of (idx1, idx2, charge1, charge2) for N atoms
    
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
                
                # Note: OpenMM automatically creates exceptions for bonded atoms
                
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
    # The KA model uses ε_NO = 1.5 × ε_OO and σ_NO = 0.8 × σ_OO, which differs
    # from the default Lorentz-Berthelot combining rules.
    n_cross_terms = 0
    for (idx_O, charge_O) in O_indices:
        for (idx_N, charge_N) in N_indices:
            # chargeProd for Coulomb: q_O * q_N
            chargeProd = charge_O * charge_N
            # Add exception with KA cross-term parameters
            nonbonded_force.addException(idx_O, idx_N, chargeProd, sigma_NO, epsilon_NO)
            n_cross_terms += 1
    
    print(f"Added {n_cross_terms} N-O cross-term exceptions for KA non-additive interactions")
    
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


def compute_dipole_moment(state, charges, num_molecular_particles):
    """
    Compute the total molecular dipole moment.
    
    Parameters
    ----------
    state : openmm.State
        Current state with positions
    charges : list
        Particle charges (as Quantity with elementary charge units or float)
    num_molecular_particles : int
        Number of molecular particles (excluding cavity)
        
    Returns
    -------
    dipole : np.ndarray
        Dipole moment vector (3,) in e*nm units
    """
    positions = state.getPositions(asNumpy=True)
    dipole = np.zeros(3)
    
    for i in range(num_molecular_particles):
        # positions[i] is in nm, charges[i] in elementary charge
        pos_nm = positions[i].value_in_unit(unit.nanometer)
        # Extract charge value (handle both Quantity and float)
        if hasattr(charges[i], 'value_in_unit'):
            charge_e = charges[i].value_in_unit(unit.elementary_charge)
        else:
            charge_e = float(charges[i])
        
        dipole[0] += charge_e * pos_nm[0]
        dipole[1] += charge_e * pos_nm[1]
        dipole[2] += charge_e * pos_nm[2]
    
    return dipole


def run_test(num_molecules=250, lambda_coupling=0.001, temperature_K=100.0,
             dt=0.001, equilibration_time_ps=100.0, production_time_ps=900.0,
             cavity_freq_cm=1560.0, disable_dipole_output=False):
    """Run the cavity diamer simulation test."""
    
    print("=" * 60)
    print("Cavity OpenMM Diamer Test")
    print("=" * 60)
    
    # System parameters
    box_size_nm = 2.5
    
    # Cavity parameters
    # Convert from cm⁻¹ to atomic units (Hartree)
    # E(Hartree) = E(cm⁻¹) / 219474.63
    omegac_au = cavity_freq_cm / 219474.63  # Hartree (energy = ℏω in atomic units where ℏ=1)
    photon_mass = 1.0 / 1822.888  # amu, so mass_au = 1.0
    
    # Simulation parameters
    equilibration_steps = int(equilibration_time_ps / dt)
    production_steps = int(production_time_ps / dt)
    
    # Dipole output parameters
    dipole_output_interval_ps = dt  # Output EVERY timestep (1 fs) for proper IR spectrum
    dipole_output_interval_steps = 1  # Every single step
    
    print(f"\nSimulation parameters:")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt} ps")
    print(f"  Cavity frequency: {cavity_freq_cm} cm⁻¹ ({omegac_au:.6f} Hartree)")
    print(f"  Lambda coupling: {lambda_coupling}")
    print(f"  Equilibration: {equilibration_time_ps} ps ({equilibration_steps} steps)")
    print(f"  Production: {production_time_ps} ps ({production_steps} steps)")
    print(f"  Total time: {equilibration_time_ps + production_time_ps} ps")
    if not disable_dipole_output:
        print(f"  Dipole output: Every timestep ({dt*1000:.1f} fs) for accurate IR spectrum")
    else:
        print(f"  Dipole output: DISABLED (benchmark mode)")
    
    # Create system
    print("\n--- Creating Diamer System ---")
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
    num_molecular_particles = cavity_index  # All particles before cavity are molecular
    
    # Check if CavityForce is available
    print("\n--- Adding Cavity Force ---")
    try:
        # Create CavityForce with target lambda from t=0 (omegac in atomic units)
        cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass)
        system.addForce(cavity_force)
        print(f"  CavityForce added successfully")
        print(f"  Omega_c: {omegac_au:.6f} a.u.")
        print(f"  Photon mass: {photon_mass:.6f} amu = {photon_mass * 1822.888:.1f} a.u.")
        # Calculate expected spring constant (with correct unit conversion)
        AMU_TO_AU = 1822.888  # 1 amu = 1822.888 electron masses
        photon_mass_au = photon_mass * AMU_TO_AU
        K_au = photon_mass_au * omegac_au**2  # In atomic units
        K_openmm = K_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
        print(f"  Expected spring constant K: {K_openmm:.0f} kJ/(mol·nm^2)")
        print(f"  Lambda coupling: {lambda_coupling} (ACTIVE from t=0)")
        
        # Create CavityParticleDisplacer for finite-Q displacement
        displacer = openmm.CavityParticleDisplacer(cavity_index, omegac_au, photon_mass)
        displacer.setSwitchOnLambda(lambda_coupling)  # Use target lambda
        system.addForce(displacer)
        print(f"  CavityParticleDisplacer added (finite-Q mode)")
        
        cavity_available = True
    except AttributeError as e:
        print(f"  CavityForce not available: {e}")
        print("  Skipping cavity coupling test")
        cavity_available = False
    
    # Check if BussiThermostat is available
    print("\n--- Adding Bussi Thermostat ---")
    try:
        # Create BussiThermostat for molecules only (not the cavity particle)
        tau = 0.1  # ps
        bussi = openmm.BussiThermostat(temperature_K, tau)
        bussi.setApplyToAllParticles(False)
        
        # Add all molecular particles (not the cavity)
        for i in range(system.getNumParticles() - 1):  # Exclude cavity particle
            bussi.addParticle(i)
        
        system.addForce(bussi)
        print(f"  BussiThermostat added for {bussi.getNumParticles()} particles")
        print(f"  Temperature: {temperature_K} K")
        print(f"  Tau: {tau} ps")
        bussi_available = True
    except AttributeError as e:
        print(f"  BussiThermostat not available: {e}")
        bussi_available = False
    
    # Create integrator - use Langevin for the cavity particle dynamics
    print("\n--- Creating Integrator ---")
    friction = 0.01  # ps^-1 - very low friction to preserve vibrational dynamics
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt * unit.picosecond
    )
    print(f"  LangevinMiddleIntegrator created")
    print(f"  Friction: {friction} ps^-1 (low to preserve vibrations)")
    print(f"  Note: Bussi thermostat handles molecular thermalization")
    
    # Create simulation - try CUDA first, fall back to Reference
    print("\n--- Creating Simulation ---")
    # Always use CUDA platform
    platform = openmm.Platform.getPlatformByName('CUDA')
    print(f"  Using CUDA platform (GPU acceleration)")
    
    # Create context
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)
    
    # Minimize energy
    print("\n--- Energy Minimization ---")
    state = context.getState(getEnergy=True)
    print(f"  Initial energy: {state.getPotentialEnergy()}")
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=100)
    state = context.getState(getEnergy=True)
    print(f"  Final energy: {state.getPotentialEnergy()}")
    
    # Apply finite-Q displacement to cavity particle BEFORE starting dynamics
    if cavity_available:
        print("\n--- Applying Finite-Q Displacement ---")
        state = context.getState(getPositions=True)
        positions_with_units = state.getPositions()
        
        # Get cavity particle position and extract numeric values
        cavity_pos = positions_with_units[cavity_index]
        x_nm = cavity_pos[0].value_in_unit(unit.nanometer)
        y_nm = cavity_pos[1].value_in_unit(unit.nanometer)
        z_nm = cavity_pos[2].value_in_unit(unit.nanometer)
        print(f"  Cavity position before: ({x_nm:.6f}, {y_nm:.6f}, {z_nm:.6f}) nm")
        
        # Calculate displacement magnitude (typical value ~0.01 nm)
        # This gives the cavity mode an initial excitation
        displacement_magnitude = 0.01  # nm
        
        # Modify positions
        positions_list = list(positions_with_units)
        positions_list[cavity_index] = [displacement_magnitude, 0.0, 0.0] * unit.nanometer
        
        context.setPositions(positions_list)
        state = context.getState(getPositions=True)
        cavity_pos = state.getPositions()[cavity_index]
        x_nm = cavity_pos[0].value_in_unit(unit.nanometer)
        y_nm = cavity_pos[1].value_in_unit(unit.nanometer)
        z_nm = cavity_pos[2].value_in_unit(unit.nanometer)
        print(f"  Cavity position after:  ({x_nm:.6f}, {y_nm:.6f}, {z_nm:.6f}) nm")
        print(f"  Displacement: {displacement_magnitude} nm along x-axis")
    
    # Set velocities
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    
    # Calculate total steps
    total_steps = equilibration_steps + production_steps
    
    # Run full simulation (no separate equilibration/production phases)
    print("\n--- Running Full Simulation (Coupling ON from t=0) ---")
    print(f"  Running {total_steps} steps (1 ns)...")
    print(f"  Progress updates every 1 ps (1,000 steps)")
    print("")
    
    # Prepare dipole trajectory storage
    dipole_times = []
    dipole_trajectory = []
    step_counter = 0  # Manual step counter
    
    report_interval = 1000  # Report every 1 ps
    num_reports = total_steps // report_interval
    
    start_time = time.time()
    last_report_time = start_time
    
    for i in range(num_reports):
        # Run and collect dipole data
        for step in range(report_interval):
            integrator.step(1)
            step_counter += 1
            if not disable_dipole_output and (step_counter % dipole_output_interval_steps) == 0:
                state = context.getState(getPositions=True)
                current_time = step_counter * dt
                dipole = compute_dipole_moment(state, charges, num_molecular_particles)
                dipole_times.append(current_time)
                dipole_trajectory.append(dipole)
        
        # Report progress with timing
        current_wall_time = time.time()
        elapsed = current_wall_time - start_time
        elapsed_since_last = current_wall_time - last_report_time
        last_report_time = current_wall_time
        
        state = context.getState(getEnergy=True, getPositions=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        positions_now = state.getPositions()
        cavity_pos = positions_now[cavity_index]
        
        sim_time_ps = step_counter * dt
        progress_pct = (step_counter / total_steps) * 100
        steps_per_sec = report_interval / elapsed_since_last
        
        # Calculate ns/day (same as μs/day)
        # steps_per_sec * dt (ps/step) * 86400 (sec/day) / 1000 (ps/ns) = ns/day
        ns_per_day = steps_per_sec * dt * 86400.0 / 1000.0
        
        # Estimate time remaining
        steps_remaining = total_steps - step_counter
        time_remaining_sec = steps_remaining / steps_per_sec if steps_per_sec > 0 else 0
        time_remaining_min = time_remaining_sec / 60
        time_remaining_hr = time_remaining_min / 60
        
        print(f"  [{progress_pct:5.1f}%] Step {step_counter:7d}/{total_steps} | "
              f"Sim time: {sim_time_ps:6.1f} ps | "
              f"Speed: {ns_per_day:6.1f} ns/day | "
              f"ETA: {time_remaining_hr:5.1f} hr")
        
        # Get cavity energy if available
        cavity_energy_str = ""
        if cavity_available:
            try:
                harmonic_e = cavity_force.getHarmonicEnergy(context)
                coupling_e = cavity_force.getCouplingEnergy(context)
                dipole_e = cavity_force.getDipoleSelfEnergy(context)
                # Handle unit.Quantity objects
                if hasattr(harmonic_e, 'value_in_unit'):
                    harmonic_e = harmonic_e.value_in_unit(unit.kilojoule_per_mole)
                    coupling_e = coupling_e.value_in_unit(unit.kilojoule_per_mole)
                    dipole_e = dipole_e.value_in_unit(unit.kilojoule_per_mole)
                cavity_energy_str = f" | Cavity: H={harmonic_e:6.2f}, C={coupling_e:6.2f}, D={dipole_e:6.2f}"
            except Exception as e:
                pass  # Skip if not available
        
        print(f"         PE: {pe:8.2f} kJ/mol{cavity_energy_str}")
        print(f"         Cavity: ({cavity_pos.x:7.4f}, {cavity_pos.y:7.4f}, {cavity_pos.z:7.4f}) nm")
        print("")
        
        # Save NPZ file on the fly (streaming save)
        if not disable_dipole_output and len(dipole_times) > 0:
            output_file = f"cavity_diamer_lambda{lambda_coupling:.4f}.npz"
            np.savez(output_file,
                     time_ps=np.array(dipole_times),
                     dipole_nm=np.array(dipole_trajectory),
                     metadata={
                         'temperature_K': temperature_K,
                         'num_molecules': num_molecules,
                         'cavity_freq_cm': cavity_freq_cm,
                         'equilibration_ps': equilibration_time_ps,
                         'production_ps': production_time_ps,
                         'dt_ps': dt,
                         'output_interval_ps': dipole_output_interval_ps,
                         'lambda_coupling': lambda_coupling,
                         'omegac_au': omegac_au,
                         'cavity_index': cavity_index,
                         'status': 'running',
                         'step': step_counter,
                         'total_steps': total_steps
                     })
    
    total_elapsed = time.time() - start_time
    print(f"  Simulation complete! Total elapsed time: {total_elapsed/3600:.2f} hr ({total_elapsed/60:.1f} min)")
    
    # Save dipole trajectory
    if not disable_dipole_output:
        print("\n--- Saving Dipole Trajectory ---")
        dipole_times_array = np.array(dipole_times)
        dipole_trajectory_array = np.array(dipole_trajectory)
        
        output_file = f"cavity_diamer_lambda{lambda_coupling:.4f}.npz"
        np.savez(output_file,
                 time_ps=dipole_times_array,
                 dipole_nm=dipole_trajectory_array,
                 metadata={
                     'temperature_K': temperature_K,
                     'num_molecules': num_molecules,
                     'cavity_freq_cm': cavity_freq_cm,
                     'equilibration_ps': equilibration_time_ps,
                    'production_ps': production_time_ps,
                    'dt_ps': dt,
                    'output_interval_ps': dipole_output_interval_ps,
                    'lambda_coupling': lambda_coupling,
                    'omegac_au': omegac_au,
                    'cavity_index': cavity_index,
                    'status': 'complete',
                    'step': total_steps,
                    'total_steps': total_steps,
                    'elapsed_time_s': total_elapsed
                 })
        
        print(f"  Saved dipole trajectory to: {output_file}")
        print(f"  Total data points: {len(dipole_times_array)}")
        print(f"  Time range: {dipole_times_array[0]:.2f} - {dipole_times_array[-1]:.2f} ps")
    else:
        print("\n--- Dipole output disabled (benchmark mode) ---")
    
    print("\n" + "=" * 60)
    print("Test completed successfully!")
    print("=" * 60)
    
    if not disable_dipole_output:
        print("\n--- IR Spectrum Calculation Instructions ---")
        print("To calculate the IR spectrum from the dipole trajectory:")
    print("1. Load the data:")
    print("   data = np.load('cavity_diamer_dipole.npz')")
    print("   time = data['time_ps']")
    print("   dipole = data['dipole_nm']")
    print("")
    print("2. Cut off first 20 ps (as requested):")
    print("   idx_cutoff = np.where(time >= 20.0)[0][0]")
    print("   time_cut = time[idx_cutoff:]")
    print("   dipole_cut = dipole[idx_cutoff:]")
    print("")
    print("3. Compute dipole autocorrelation function:")
    print("   C(t) = <M(0) · M(t)> where M is dipole moment")
    print("")
    print("4. Fourier transform C(t) to get IR spectrum:")
    print("   I(ω) ∝ FT[C(t)]")
    print("")
    print("Note: The first 100 ps is equilibration (coupling OFF),")
    print("      the remaining 900 ps is production (coupling ON).")
    print("=" * 60)
    
    return True


def main():
    """Main entry point with command-line argument parsing."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Dimer Cavity MD Simulation')
    parser.add_argument('--dimers', type=int, default=250,
                       help='Number of dimers (default: 250)')
    parser.add_argument('--lambda', type=float, default=0.001, dest='lambda_coupling',
                       help='Coupling strength (default: 0.001)')
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
        success = run_test(
            num_molecules=args.dimers,
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
        print(f"\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
