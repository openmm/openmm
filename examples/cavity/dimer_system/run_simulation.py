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
from pathlib import Path
import time
from datetime import date

try:
    from openmm import openmm
    from openmm import unit
    from openmm.app import Simulation, StateDataReporter, ForceField, Topology, Element, CutoffPeriodic, PDBFile
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

# Reference density for constant-density scaling: N_ref dimers in cubic box of side L_ref Bohr
DIAMER_REF_N = 250
DIAMER_REF_L_BOHR = 40


def box_size_nm_at_constant_density(num_molecules, N_ref=DIAMER_REF_N, L_ref_bohr=DIAMER_REF_L_BOHR):
    """
    Box edge length in nm so that number density is constant.

    Reference: N_ref dimers in a cubic box of side L_ref_bohr Bohr. Then
    L(N) = L_ref_bohr * (N / N_ref)^(1/3) Bohr, converted to nm.

    Parameters
    ----------
    num_molecules : int
        Number of dimers
    N_ref : int
        Reference number of dimers (default 250)
    L_ref_bohr : float
        Reference box length in Bohr (default 40)

    Returns
    -------
    float
        Box edge length in nm
    """
    L_bohr = L_ref_bohr * (num_molecules / N_ref) ** (1 / 3)
    return L_bohr * BOHR_TO_NM


def _system_xml_paths(num_molecules, fraction_OO, box_size_nm, seed):
    tag = f"diamer_{num_molecules}_OO{fraction_OO:.2f}_box{box_size_nm:.3f}_seed{seed}"
    base_dir = Path(__file__).parent
    return base_dir / f"{tag}.xml", base_dir / f"{tag}_positions.npz"


def create_diamer_topology_and_positions(num_molecules=50, fraction_OO=0.8, box_size_nm=2.0, seed=42,
                                         include_cavity=False):
    """
    Create topology and positions for O-O and N-N dimers (same layout as create_diamer_system_in_code).
    For use with ForceField('diamer_forcefield.xml').createSystem(topology).
    If include_cavity=True, append one residue "CAV" with one atom at (0,0,0) for LennardJonesForce compatibility.
    Returns (topology, positions) with positions as list of Vec3 in nm.
    """
    np.random.seed(seed)
    # Bond lengths in nm (same as create_diamer_system_in_code)
    r0_OO_au = 2.281655158  # Bohr
    r0_NN_au = 2.0743522177
    r0_OO = r0_OO_au * BOHR_TO_NM
    r0_NN = r0_NN_au * BOHR_TO_NM

    topology = Topology()
    chain = topology.addChain()
    positions = []
    num_OO = int(fraction_OO * num_molecules)
    side = int(np.ceil(num_molecules**(1/3)))
    spacing = box_size_nm / side

    mol_idx = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if mol_idx >= num_molecules:
                    break
                is_OO = mol_idx < num_OO
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing
                z = (k + 0.5) * spacing
                theta = np.random.rand() * 2 * np.pi
                phi = np.arccos(2 * np.random.rand() - 1)
                direction = np.array([
                    np.sin(phi) * np.cos(theta),
                    np.sin(phi) * np.sin(theta),
                    np.cos(phi)
                ])
                r0 = r0_OO if is_OO else r0_NN
                r1 = np.array([x, y, z]) - 0.5 * r0 * direction
                r2 = np.array([x, y, z]) + 0.5 * r0 * direction
                res_name = "OO" if is_OO else "NN"
                elem = Element.getBySymbol("O") if is_OO else Element.getBySymbol("N")
                residue = topology.addResidue(res_name, chain)
                a1 = topology.addAtom("A", elem, residue)
                a2 = topology.addAtom("B", elem, residue)
                topology.addBond(a1, a2)
                positions.append(openmm.Vec3(*r1) * unit.nanometer)
                positions.append(openmm.Vec3(*r2) * unit.nanometer)
                mol_idx += 1
            if mol_idx >= num_molecules:
                break
        if mol_idx >= num_molecules:
            break
    if include_cavity:
        residue_cav = topology.addResidue("CAV", chain)
        topology.addAtom("Q", Element.getBySymbol("He"), residue_cav)
        positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)
    return topology, positions


def create_diamer_system_in_code(num_molecules=50, fraction_OO=0.8, box_size_nm=2.0,
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


def create_diamer_system(num_molecules=50, fraction_OO=0.8, box_size_nm=2.0,
                         temperature_K=100.0, seed=42):
    """
    Create or load a system of O-O and N-N dimers from XML.
    """
    xml_path, pos_path = _system_xml_paths(num_molecules, fraction_OO, box_size_nm, seed)
    if not (xml_path.exists() and pos_path.exists()):
        raise FileNotFoundError(
            f"Missing system files. Expected {xml_path.name} and {pos_path.name} in {xml_path.parent}"
        )
    print(f"Loading system from XML: {xml_path.name}")
    with open(xml_path, 'r') as f:
        system = openmm.XmlSerializer.deserialize(f.read())
    data = np.load(pos_path, allow_pickle=True)
    positions_nm = data["positions_nm"]
    positions = [openmm.Vec3(*pos) * unit.nanometer for pos in positions_nm]
    if system.getNumParticles() != len(positions):
        raise ValueError("XML system and positions count do not match.")
    return system, positions


def create_diamer_system_from_forcefield(num_molecules=50, fraction_OO=0.8, box_size_nm=2.0,
                                         seed=42, ff_dir=None, include_cavity=True):
    """
    Build dimer system from OpenMM ForceField XML (diamer_forcefield.xml).
    N-O Kob-Andersen cross terms are defined in XML via LennardJonesForce/NBFixPair
    (no in-code exceptions). If include_cavity=True, topology includes one CAV residue
    at (0,0,0) so LennardJonesForce/CustomNonbondedForce get the right particle count.
    Returns (system, positions, topology).
    """
    if ff_dir is None:
        ff_dir = Path(__file__).parent
    ff_path = ff_dir / "diamer_forcefield.xml"
    if not ff_path.exists():
        raise FileNotFoundError(f"ForceField XML not found: {ff_path}")
    forcefield = ForceField(str(ff_path))
    topology, positions = create_diamer_topology_and_positions(
        num_molecules=num_molecules, fraction_OO=fraction_OO,
        box_size_nm=box_size_nm, seed=seed, include_cavity=include_cavity
    )
    vectors = (
        openmm.Vec3(box_size_nm, 0, 0) * unit.nanometer,
        openmm.Vec3(0, box_size_nm, 0) * unit.nanometer,
        openmm.Vec3(0, 0, box_size_nm) * unit.nanometer,
    )
    topology.setPeriodicBoxVectors(vectors)
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=CutoffPeriodic,
        nonbondedCutoff=0.9 * unit.nanometer,
    )
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_size_nm, 0, 0),
        openmm.Vec3(0, box_size_nm, 0),
        openmm.Vec3(0, 0, box_size_nm),
    )
    num_OO = int(fraction_OO * num_molecules)
    print(f"Created system with {system.getNumParticles()} particles ({num_molecules} molecules)")
    print(f"  O-O dimers: {num_OO}")
    print(f"  N-N dimers: {num_molecules - num_OO}")
    print(f"  Box size: {box_size_nm} nm")
    print(f"  N-O cross terms from XML (LennardJonesForce NBFixPair)")
    if include_cavity:
        cavity_index = 2 * num_molecules
        # Cavity feels only CavityForce: no bonded/nonbonded with anybody.
        # Exclude cavity from LennardJones; add same pairs as exceptions in NonbondedForce
        # so OpenMM's "identical exceptions" requirement is satisfied.
        nbf = None
        ljf = None
        for f in system.getForces():
            if isinstance(f, openmm.NonbondedForce):
                nbf = f
            if isinstance(f, openmm.CustomNonbondedForce) and f.getName() == 'LennardJones':
                ljf = f
        for i in range(system.getNumParticles()):
            if i == cavity_index:
                continue
            if ljf is not None:
                ljf.addExclusion(cavity_index, i)
            if nbf is not None:
                nbf.addException(cavity_index, i, 0.0, 0.1, 0.0)  # no Coulomb, no LJ
        if ljf is not None or nbf is not None:
            print(f"  Cavity excluded from nonbonded (only feels CavityForce)")
        print(f"  Cavity particle at index {cavity_index} (from topology)")
        return system, positions, topology, cavity_index
    return system, positions, topology


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
             cavity_freq_cm=1560.0, disable_dipole_output=False,
             box_size_nm=2.5, fraction_OO=0.8, friction=0.01, minimize=True,
             output_file=None, seed=42, report_interval_steps=1000,
             pdb_file=None, pdb_interval_steps=1000,
             enable_fkt=False, fkt_kmag=113.4, fkt_num_wavevectors=50,
             fkt_reference_interval_ps=1.0, fkt_max_refs=10,
             fkt_output_period_ps=1.0, fkt_output_prefix=None):
    """Run the cavity diamer simulation test."""
    
    print("=" * 60)
    print("Cavity OpenMM Diamer Test")
    print("=" * 60)
    
    # System parameters (box_size_nm, fraction_OO come from arguments)
    
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
    
    g_collective = lambda_coupling * np.sqrt(num_molecules)
    # Strong coupling + coarse dt can blow up the cavity mode (light particle, stiff coupling)
    if g_collective > 0.5 and dt > 0.001:
        print(f"\n  *** WARNING: g = {g_collective:.3g} with dt = {dt*1000:.0f} fs may be unstable. Use --dt 0.001 (1 fs) for strong coupling. ***")
    print(f"\nSimulation parameters:")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt} ps")
    print(f"  Cavity frequency: {cavity_freq_cm} cm⁻¹ ({omegac_au:.6f} Hartree)")
    print(f"  Coupling: λ = {lambda_coupling:.6g} (g = λ√N = {g_collective:.6g})")
    print(f"  Equilibration: {equilibration_time_ps} ps ({equilibration_steps} steps)")
    print(f"  Production: {production_time_ps} ps ({production_steps} steps)")
    print(f"  Total time: {equilibration_time_ps + production_time_ps} ps")
    if not disable_dipole_output:
        print(f"  Dipole output: Every timestep ({dt*1000:.1f} fs) from t=0 for full run (equilibration + production)")
    else:
        print(f"  Dipole output: DISABLED (benchmark mode)")
    print(f"  Console report interval: every {report_interval_steps} steps")
    if pdb_file:
        print(f"  PDB trajectory: {pdb_file} (every {pdb_interval_steps} steps)")
    else:
        print(f"  PDB trajectory: disabled")
    fkt_prefix_resolved = None
    if enable_fkt:
        fkt_prefix_resolved = (fkt_output_prefix if fkt_output_prefix is not None else
                               (Path(output_file).with_suffix("").name if output_file else f"fkt_lambda{lambda_coupling:.4f}"))
        print(f"  F(k,t): enabled, prefix={fkt_prefix_resolved}, k={fkt_kmag} nm⁻¹")
        print(f"    reference interval (new file every): {fkt_reference_interval_ps} ps")
        print(f"    max reference files (FIFO drop): {fkt_max_refs}")
        print(f"    F(k,t) output period: {fkt_output_period_ps} ps")
    
    # Create system from ForceField XML (diamer_forcefield.xml); cavity in topology for LJ force compatibility
    print("\n--- Creating Diamer System ---")
    result = create_diamer_system_from_forcefield(
        num_molecules=num_molecules,
        fraction_OO=fraction_OO,
        box_size_nm=box_size_nm,
        seed=seed,
        include_cavity=True,
    )
    if len(result) == 4:
        system, positions, topology, cavity_index = result
        num_molecular_particles = cavity_index
        system.setParticleMass(cavity_index, photon_mass)
        print("\n--- Cavity (from topology) ---")
        print(f"  Cavity particle at index {cavity_index}, mass set to {photon_mass:.6f} amu")
    else:
        system, positions, topology = result
        print("\n--- Adding Cavity Particle ---")
        cavity_index = add_cavity_particle(system, positions, omegac_au, photon_mass)
        num_molecular_particles = cavity_index
    
    # Store charges for dipole calculation
    charges = []
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for i in range(system.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(i)
                charges.append(charge)
    
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
    
    # Minimize energy (optional)
    if minimize:
        print("\n--- Energy Minimization ---")
        state = context.getState(getEnergy=True)
        print(f"  Initial energy: {state.getPotentialEnergy()}")
        openmm.LocalEnergyMinimizer.minimize(context, maxIterations=100)
        state = context.getState(getEnergy=True)
        print(f"  Final energy: {state.getPotentialEnergy()}")
    else:
        print("\n--- Skipping energy minimization (--no-minimize) ---")
    
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
    
    report_interval = report_interval_steps
    num_reports = total_steps // report_interval
    
    # PDB output (molecular atoms only; exclude cavity)
    pdb_handle = None
    pdb_model_index = 0
    if pdb_file:
        pdb_handle = open(pdb_file, 'w')
        # Minimal header (skip PDBFile.writeHeader to avoid Quantity formatting of box vectors)
        try:
            from openmm import Platform
            _v = Platform.getOpenMMVersion()
        except Exception:
            _v = "unknown"
        pdb_handle.write("REMARK   1 CREATED WITH OPENMM %s, %s\n" % (_v, str(date.today())))
        pdb_handle.write("REMARK   2 Molecular atoms only (cavity particle omitted)\n")
    
    start_time = time.time()
    last_report_time = start_time
    
    fkt_tracker = None
    fkt_interval_steps = 1
    if enable_fkt:
        _here = Path(__file__).resolve().parent
        if str(_here) not in sys.path:
            sys.path.insert(0, str(_here))
        from fkt_tracker import FKTTracker
        fkt_tracker = FKTTracker(
            kmag_nm_inv=fkt_kmag,
            num_wavevectors=fkt_num_wavevectors,
            reference_interval_ps=fkt_reference_interval_ps,
            max_references=fkt_max_refs,
            output_period_ps=fkt_output_period_ps,
            output_prefix=fkt_prefix_resolved,
        )
        fkt_interval_steps = max(1, int(round(fkt_output_period_ps / dt)))
    
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
            # PDB frame output (molecular atoms only; PDBFile expects angstroms)
            if pdb_handle and (step_counter % pdb_interval_steps) == 0:
                _state = context.getState(getPositions=True)
                _pos_all = _state.getPositions(asNumpy=True)
                _pos_mol = _pos_all[:num_molecular_particles].value_in_unit(unit.angstroms)
                PDBFile.writeModel(topology, _pos_mol, pdb_handle, modelIndex=pdb_model_index + 1)
                pdb_model_index += 1
            # F(k,t) intermediate scattering function (molecular positions only)
            if enable_fkt and fkt_tracker is not None and (step_counter % fkt_interval_steps) == 0:
                _st = context.getState(getPositions=True)
                _pos_all = _st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                pos_mol = _pos_all[:num_molecular_particles]
                current_time_ps = step_counter * dt
                fkt_tracker.update(current_time_ps, pos_mol)
        
        # Report progress with timing
        current_wall_time = time.time()
        elapsed = current_wall_time - start_time
        elapsed_since_last = current_wall_time - last_report_time
        last_report_time = current_wall_time
        
        state = context.getState(getEnergy=True, getPositions=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        positions_now = state.getPositions()
        cavity_pos = positions_now[cavity_index]
        cx = cavity_pos[0].value_in_unit(unit.nanometer) if hasattr(cavity_pos[0], 'value_in_unit') else float(cavity_pos[0])
        cy = cavity_pos[1].value_in_unit(unit.nanometer) if hasattr(cavity_pos[1], 'value_in_unit') else float(cavity_pos[1])
        cz = cavity_pos[2].value_in_unit(unit.nanometer) if hasattr(cavity_pos[2], 'value_in_unit') else float(cavity_pos[2])
        
        sim_time_ps = step_counter * dt
        progress_pct = (step_counter / total_steps) * 100
        steps_per_sec = report_interval / elapsed_since_last
        
        # Simulation time per day: steps_per_sec * dt (ps/step) * 86400 (s/day) = ps/day
        # 1 μs = 1000 ns = 1e6 ps, so μs/day = (ps/day) / 1e6
        us_per_day = steps_per_sec * dt * 86400.0 / 1e6
        
        # Estimate time remaining
        steps_remaining = total_steps - step_counter
        time_remaining_sec = steps_remaining / steps_per_sec if steps_per_sec > 0 else 0
        time_remaining_min = time_remaining_sec / 60
        time_remaining_hr = time_remaining_min / 60
        
        print(f"  [{progress_pct:5.1f}%] Step {step_counter:7d}/{total_steps} | "
              f"Sim time: {sim_time_ps:6.1f} ps | "
              f"Speed: {us_per_day:6.2f} μs/day | "
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
        print(f"         Cavity: ({cx:7.4f}, {cy:7.4f}, {cz:7.4f}) nm")
        print("")
        
        # Save NPZ file on the fly (streaming save)
        if not disable_dipole_output and len(dipole_times) > 0:
            _out = output_file or f"cavity_diamer_lambda{lambda_coupling:.4f}.npz"
            np.savez(_out,
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
    
    # Close PDB trajectory
    if pdb_handle:
        PDBFile.writeFooter(topology, pdb_handle)
        pdb_handle.close()
        print(f"  PDB trajectory written to {pdb_file} ({pdb_model_index} frames)")
    
    if enable_fkt and fkt_tracker is not None:
        fkt_tracker.finalize()
    
    total_elapsed = time.time() - start_time
    print(f"  Simulation complete! Total elapsed time: {total_elapsed/3600:.2f} hr ({total_elapsed/60:.1f} min)")
    
    # Save dipole trajectory
    if not disable_dipole_output:
        print("\n--- Saving Dipole Trajectory ---")
        dipole_times_array = np.array(dipole_times)
        dipole_trajectory_array = np.array(dipole_trajectory)
        _out = output_file or f"cavity_diamer_lambda{lambda_coupling:.4f}.npz"
        np.savez(_out,
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
        
        print(f"  Saved dipole trajectory to: {_out}")
        print(f"  Total data points: {len(dipole_times_array)}")
        if len(dipole_times_array) > 0:
            print(f"  Time range: {dipole_times_array[0]:.2f} - {dipole_times_array[-1]:.2f} ps")
        else:
            print("  Time range: (no dipole data recorded)")
    else:
        print("\n--- Dipole output disabled (benchmark mode) ---")
    
    print("\n" + "=" * 60)
    print("Test completed successfully!")
    print("=" * 60)
    
    if not disable_dipole_output:
        print("\n--- IR Spectrum Calculation Instructions ---")
        print("The dipole trajectory includes equilibration and production (saved from t=0).")
        print("To calculate the IR spectrum:")
        print("1. Load the data (use your output filename or default cavity_diamer_lambda<λ>.npz):")
        print("   data = np.load('cavity_diamer_lambda0.0010.npz', allow_pickle=True)")
        print("   time = data['time_ps']")
        print("   dipole = data['dipole_nm']")
        print("   meta = data['metadata'].item() if 'metadata' in data else {}")
        print("   equil_ps = meta.get('equilibration_ps', 100.0)")
        print("")
        print("2. Optionally use only production: trim the first equil_ps when computing spectra:")
        print("   idx_cutoff = np.where(time >= equil_ps)[0][0]")
        print("   time_cut = time[idx_cutoff:]")
        print("   dipole_cut = dipole[idx_cutoff:]")
        print("")
        print("3. Compute dipole autocorrelation function:")
        print("   C(t) = <M(0) · M(t)> where M is dipole moment")
        print("")
        print("4. Fourier transform C(t) to get IR spectrum:")
        print("   I(ω) ∝ FT[C(t)]")
        print("")
        print("Note: metadata in the .npz contains 'equilibration_ps' and 'production_ps' for trimming.")
    print("=" * 60)
    
    return True


def main():
    """Main entry point with command-line argument parsing."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Cavity OpenMM Diamer Simulation (CavityForce + Bussi thermostat)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default parameters (250 dimers, 100 ps equil + 900 ps prod)
  python run_simulation.py

  # Short run for testing
  python run_simulation.py --dimers 8 --equil 0.1 --prod 0.1

  # Constant density (250 dimers in 40 Bohr box reference); box size computed from --dimers
  python run_simulation.py --dimers 100 --constant-density

  # Custom system and coupling (--lambda sets λ in a.u.; else use --g for g=λ√N)
  python run_simulation.py --dimers 100 --box-size 2.0 --fraction-OO 0.75 --lambda 0.0005

  # Total time only (10% equil, 90% prod)
  python run_simulation.py --time 200

  # Skip minimization, custom output
  python run_simulation.py --no-minimize --output my_dipole.npz

  # F(k,t) ISF: new reference every 200 ps, keep up to 5, output every 1 ps
  python run_simulation.py --enable-fkt --fkt-ref-interval-ps 200 --fkt-max-refs 5 --fkt-output-period-ps 1
        """,
    )
    # System parameters
    parser.add_argument('--dimers', type=int, default=250,
                        help='Number of dimers (default: 250)')
    parser.add_argument('--constant-density', action='store_true', dest='constant_density',
                        help='Scale box so density is constant: 250 dimers in 40 Bohr box (--box-size ignored)')
    parser.add_argument('--box-size', type=float, default=DIAMER_REF_L_BOHR * BOHR_TO_NM,
                        help='Box edge length in nm (default: 2.117 = 40 Bohr, cav-hoomd parity). Ignored if --constant-density.')
    parser.add_argument('--fraction-OO', type=float, default=0.8, dest='fraction_OO',
                        help='Fraction of O-O dimers vs N-N (default: 0.8)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for initial positions (default: 42)')
    # Time parameters
    parser.add_argument('--dt', type=float, default=0.001,
                        help='Timestep in ps (default: 0.001 = 1 fs). Use 1 fs or smaller for strong coupling (g > 0.5)')
    parser.add_argument('--time', type=float, default=None,
                        help='Total simulation time in ps; if set, uses 10%% equil + 90%% prod (overrides --equil/--prod)')
    parser.add_argument('--equil', type=float, default=100.0,
                        help='Equilibration time in ps (default: 100)')
    parser.add_argument('--prod', type=float, default=900.0,
                        help='Production time in ps (default: 900)')
    # Cavity parameters
    parser.add_argument('--cavity-freq', type=float, default=1560.0,
                        help='Cavity frequency in cm⁻¹ (default: 1560 for O-O stretch)')
    parser.add_argument('--lambda', type=float, default=None, dest='lambda_au',
                        help='Coupling λ in a.u.; overrides --g when set (e.g. 0.001 for cav-hoomd parity)')
    parser.add_argument('--g', type=float, default=0.0158,
                        help='Coupling g = λ√N when --lambda not set; λ = g/√N (default: 0.0158 ⇒ λ≈0.001 for N=250)')
    # Simulation control
    parser.add_argument('--temp', type=float, default=100.0,
                        help='Temperature in K (default: 100)')
    parser.add_argument('--friction', type=float, default=0.01,
                        help='Langevin friction in ps⁻¹ (default: 0.01)')
    parser.add_argument('--no-minimize', action='store_true',
                        help='Skip energy minimization')
    parser.add_argument('--no-dipole', action='store_true', dest='no_dipole',
                        help='Disable dipole output for speed benchmark')
    parser.add_argument('--output', type=str, default=None,
                        help='Output filename for dipole trajectory .npz (default: cavity_diamer_lambda<λ>.npz)')
    # Console and PDB output
    parser.add_argument('--report-interval', type=int, default=1000, dest='report_interval',
                        help='Console progress report interval in steps (default: 1000)')
    parser.add_argument('--pdb', type=str, default=None, dest='pdb_file',
                        help='Write PDB trajectory to this file (molecular atoms only; disable if not set)')
    parser.add_argument('--pdb-interval', type=int, default=1000, dest='pdb_interval',
                        help='PDB frame write interval in steps (default: 1000)')
    # F(k,t) intermediate scattering function
    parser.add_argument('--enable-fkt', action='store_true', dest='enable_fkt',
                        help='Enable F(k,t) density correlation output')
    parser.add_argument('--fkt-kmag', type=float, default=113.4, dest='fkt_kmag',
                        help='F(k,t) k-magnitude in nm⁻¹ (default: 113.4 ~ 6.0 a.u.)')
    parser.add_argument('--fkt-num-wavevectors', type=int, default=50, dest='fkt_num_wavevectors',
                        help='Number of F(k,t) wavevectors (default: 50)')
    parser.add_argument('--fkt-ref-interval-ps', type=float, default=1.0, dest='fkt_reference_interval_ps',
                        help='F(k,t) interval between reference times t_w in ps. A new reference (and new F(k,t) file) is added every this many ps (default: 1.0)')
    parser.add_argument('--fkt-max-refs', type=int, default=10, dest='fkt_max_refs',
                        help='F(k,t) maximum number of reference states kept. One F(k,t) file per reference; when this many exist, the oldest is dropped and we stop writing to it. Total files created = ceil(runtime_ps / fkt-ref-interval-ps) (default: 10)')
    parser.add_argument('--fkt-output-period-ps', type=float, default=1.0, dest='fkt_output_period_ps',
                        help='F(k,t) output period in ps (default: 1.0)')
    parser.add_argument('--fkt-output-prefix', type=str, default=None, dest='fkt_output_prefix',
                        help='Prefix for F(k,t) files (default: from --output or fkt_lambda<λ>)')
    
    args = parser.parse_args()

    # If --time is set, override equil and prod
    equil_ps = args.equil
    prod_ps = args.prod
    if args.time is not None:
        equil_ps = 0.1 * args.time
        prod_ps = 0.9 * args.time

    # Box size: constant-density scaling or explicit
    if args.constant_density:
        box_size_nm = box_size_nm_at_constant_density(args.dimers)
        print(f"Constant density: {args.dimers} dimers → box {box_size_nm:.4f} nm (ref: {DIAMER_REF_N} in {DIAMER_REF_L_BOHR} Bohr)")
    else:
        box_size_nm = args.box_size

    # Coupling: --lambda sets λ directly; else λ = g/√N
    lambda_coupling = args.lambda_au if args.lambda_au is not None else (args.g / np.sqrt(args.dimers))

    try:
        success = run_test(
            num_molecules=args.dimers,
            lambda_coupling=lambda_coupling,
            temperature_K=args.temp,
            dt=args.dt,
            equilibration_time_ps=equil_ps,
            production_time_ps=prod_ps,
            cavity_freq_cm=args.cavity_freq,
            disable_dipole_output=args.no_dipole,
            box_size_nm=box_size_nm,
            fraction_OO=args.fraction_OO,
            friction=args.friction,
            minimize=not args.no_minimize,
            output_file=args.output,
            seed=args.seed,
            report_interval_steps=args.report_interval,
            pdb_file=args.pdb_file,
            pdb_interval_steps=args.pdb_interval,
            enable_fkt=args.enable_fkt,
            fkt_kmag=args.fkt_kmag,
            fkt_num_wavevectors=args.fkt_num_wavevectors,
            fkt_reference_interval_ps=args.fkt_reference_interval_ps,
            fkt_max_refs=args.fkt_max_refs,
            fkt_output_period_ps=args.fkt_output_period_ps,
            fkt_output_prefix=args.fkt_output_prefix,
        )
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
