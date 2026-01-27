#!/usr/bin/env python3
"""
Cavity Particle Utilities
==========================

Shared utilities for cavity-coupled molecular dynamics simulations.
Provides common functions for adding cavity particles, setting up coupling,
and managing cavity-related calculations.
"""

import numpy as np
import openmm
from openmm import unit

# Physical constants
HARTREE_TO_KJ_MOL = 2625.5  # 1 Hartree = 2625.5 kJ/mol
HARTREE_TO_CM = 219474.63   # 1 Hartree = 219474.63 cm⁻¹
AMU_TO_AU = 1822.888        # 1 amu = 1822.888 atomic units of mass


def wavenumber_to_hartree(wavenumber_cm):
    """
    Convert wavenumber (cm⁻¹) to energy in Hartree (atomic units).
    
    Parameters
    ----------
    wavenumber_cm : float
        Wavenumber in cm⁻¹
        
    Returns
    -------
    float
        Energy in Hartree (atomic units)
    """
    return wavenumber_cm / HARTREE_TO_CM


def hartree_to_wavenumber(hartree):
    """
    Convert energy in Hartree to wavenumber (cm⁻¹).
    
    Parameters
    ----------
    hartree : float
        Energy in Hartree (atomic units)
        
    Returns
    -------
    float
        Wavenumber in cm⁻¹
    """
    return hartree * HARTREE_TO_CM


def add_cavity_particle(system, positions, omegac_au, photon_mass_amu=None):
    """
    Add a cavity photon particle to the system.
    
    This function adds a cavity particle at the origin with the specified
    photon mass. The particle is added to the NonbondedForce with no
    interactions (charge=0, epsilon=0).
    
    Parameters
    ----------
    system : openmm.System
        The OpenMM system to add the cavity particle to
    positions : list
        List of positions to append to
    omegac_au : float
        Cavity frequency in atomic units (Hartree)
    photon_mass_amu : float, optional
        Photon mass in amu. If None, uses 1.0/1822.888 (1 a.u.)
        
    Returns
    -------
    cavity_index : int
        Index of the added cavity particle
        
    Examples
    --------
    >>> system = openmm.System()
    >>> positions = []
    >>> omegac_au = wavenumber_to_hartree(3663.0)  # OH stretch
    >>> cavity_idx = add_cavity_particle(system, positions, omegac_au)
    """
    if photon_mass_amu is None:
        photon_mass_amu = 1.0 / AMU_TO_AU  # 1 a.u. of mass
    
    # Add cavity particle at origin
    cavity_index = system.addParticle(photon_mass_amu * unit.amu)
    positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)
    
    # Add to nonbonded force with no interactions
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.addParticle(0.0, 0.1 * unit.nanometer, 0.0)
            break
    
    omegac_cm = hartree_to_wavenumber(omegac_au)
    print(f"  Cavity particle index: {cavity_index}")
    print(f"  Omega_c: {omegac_au:.6f} a.u. = {omegac_cm:.1f} cm⁻¹")
    print(f"  Photon mass: {photon_mass_amu:.6f} amu")
    
    return cavity_index


def setup_cavity_coupling(system, cavity_index, omegac_au, lambda_coupling, 
                         photon_mass_amu=None):
    """
    Set up cavity-matter coupling forces.
    
    This function adds two forces:
    1. CavityForce: Implements the light-matter coupling Hamiltonian
    2. CavityParticleDisplacer: Handles finite-q corrections
    
    Parameters
    ----------
    system : openmm.System
        The OpenMM system
    cavity_index : int
        Index of the cavity particle
    omegac_au : float
        Cavity frequency in Hartree (atomic units)
    lambda_coupling : float
        Coupling strength (dimensionless)
    photon_mass_amu : float, optional
        Photon mass in amu. If None, uses 1.0/1822.888 (1 a.u.)
        
    Returns
    -------
    cavity_force : openmm.CavityForce
        The cavity force object
    displacer : openmm.CavityParticleDisplacer
        The cavity particle displacer object
        
    Examples
    --------
    >>> cavity_force, displacer = setup_cavity_coupling(
    ...     system, cavity_idx, omegac_au, lambda_coupling=0.01
    ... )
    """
    if photon_mass_amu is None:
        photon_mass_amu = 1.0 / AMU_TO_AU
    
    print(f"\n--- Setting Up Cavity Coupling ---")
    print(f"  Lambda: {lambda_coupling}")
    print(f"  Coupling is ON from t=0")
    
    # CavityForce - dipole-cavity coupling
    cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, 
                                     photon_mass_amu)
    system.addForce(cavity_force)
    print(f"  CavityForce added with lambda={lambda_coupling}")
    
    # CavityParticleDisplacer - finite-q correction
    displacer = openmm.CavityParticleDisplacer(cavity_index, omegac_au, 
                                               photon_mass_amu)
    displacer.setSwitchOnStep(0)
    displacer.setSwitchOnLambda(lambda_coupling)
    system.addForce(displacer)
    print(f"  CavityParticleDisplacer added with lambda={lambda_coupling}")
    
    return cavity_force, displacer


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
        
    Examples
    --------
    >>> state = context.getState(getPositions=True)
    >>> dipole = compute_dipole_moment(state, charges, n_atoms)
    >>> dipole_magnitude = np.linalg.norm(dipole)
    """
    positions = state.getPositions(asNumpy=True)
    dipole = np.zeros(3)
    
    for i in range(num_molecular_particles):
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


def setup_bussi_thermostat(system, temperature_K, num_molecular_particles, 
                           tau_ps=1.0):
    """
    Set up Bussi thermostat for molecular particles (excluding cavity).
    
    Parameters
    ----------
    system : openmm.System
        The OpenMM system
    temperature_K : float
        Target temperature in Kelvin
    num_molecular_particles : int
        Number of molecular (non-cavity) particles
    tau_ps : float, optional
        Bussi thermostat time constant in picoseconds (default: 1.0)
        
    Returns
    -------
    bussi : openmm.BussiThermostat
        The Bussi thermostat object
        
    Examples
    --------
    >>> bussi = setup_bussi_thermostat(system, 300.0, n_water_atoms)
    """
    print(f"\n--- Setting Up Bussi Thermostat ---")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Time constant: {tau_ps} ps")
    
    bussi = openmm.BussiThermostat(temperature_K, tau_ps)
    bussi.setApplyToAllParticles(False)
    
    # Add all molecular particles (exclude cavity)
    for i in range(num_molecular_particles):
        bussi.addParticle(i)
    
    system.addForce(bussi)
    print(f"  Applied to {num_molecular_particles} molecular particles")
    print(f"  Cavity particle excluded from thermostat")
    
    return bussi


def get_cavity_displacement(state, cavity_index):
    """
    Get the current displacement of the cavity particle from origin.
    
    Parameters
    ----------
    state : openmm.State
        Current state with positions
    cavity_index : int
        Index of the cavity particle
        
    Returns
    -------
    displacement : np.ndarray
        Displacement vector (3,) in nm
    magnitude : float
        Displacement magnitude in nm
        
    Examples
    --------
    >>> state = context.getState(getPositions=True)
    >>> disp, mag = get_cavity_displacement(state, cavity_idx)
    >>> print(f"Cavity displacement: {mag:.4f} nm")
    """
    positions = state.getPositions(asNumpy=True)
    cavity_pos = positions[cavity_index].value_in_unit(unit.nanometer)
    displacement = np.array(cavity_pos)
    magnitude = np.linalg.norm(displacement)
    
    return displacement, magnitude


def cavity_info_dict(omegac_au=None, omegac_cm=None, lambda_coupling=None, 
                     photon_mass_amu=None):
    """
    Create a dictionary with cavity parameters for logging/metadata.
    
    Parameters
    ----------
    omegac_au : float, optional
        Cavity frequency in Hartree
    omegac_cm : float, optional
        Cavity frequency in cm⁻¹
    lambda_coupling : float, optional
        Coupling strength
    photon_mass_amu : float, optional
        Photon mass in amu
        
    Returns
    -------
    dict
        Dictionary with cavity parameters
    """
    info = {}
    
    if omegac_au is not None:
        info['omega_cavity_au'] = omegac_au
        info['omega_cavity_cm'] = hartree_to_wavenumber(omegac_au)
    elif omegac_cm is not None:
        info['omega_cavity_cm'] = omegac_cm
        info['omega_cavity_au'] = wavenumber_to_hartree(omegac_cm)
    
    if lambda_coupling is not None:
        info['lambda_coupling'] = lambda_coupling
    
    if photon_mass_amu is not None:
        info['photon_mass_amu'] = photon_mass_amu
        info['photon_mass_au'] = photon_mass_amu * AMU_TO_AU
    
    return info
