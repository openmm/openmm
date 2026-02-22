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

F(k,t) and T=0: F(k,t) uses wrapped molecular positions (cav-hoomd parity). At nearly
0 K the structure is effectively frozen, so F(k,t) stays close to its initial value.
With setVelocitiesToTemperature(0), the first Verlet step still adds velocity from
forces, so the system is not strictly frozen unless using a special protocol (e.g.
minimizer only). With a thermostat at very low T, motion is small and F(k,t) does
not decorrelate from coordinate drift.
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
    from openmm.app import Simulation, StateDataReporter, ForceField, Topology, Element, PME, PDBFile
    print("✓ OpenMM loaded successfully")
except ImportError as e:
    print(f"Error importing OpenMM: {e}")
    print("Make sure OpenMM is installed and the Python path is correct.")
    sys.exit(1)

try:
    import gsd.hoomd
    HAS_GSD = True
except ImportError:
    HAS_GSD = False

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

# Nonbonded cutoff matching cav-hoomd: rcut = 15 Bohr for both LJ and PPPM real-space
DIAMER_RCUT_BOHR = 15
DIAMER_RCUT_NM = DIAMER_RCUT_BOHR * BOHR_TO_NM  # 0.7938 nm


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


def _wrap_positions_into_box(pos_nm, box_nm):
    """
    Wrap positions into primary image [0, L) per dimension (orthorhombic).
    Handles GSD in [0,L) or [-L/2, L/2) convention so OpenMM sees consistent primary image.
    """
    pos = np.asarray(pos_nm, dtype=np.float64)
    box = np.asarray(box_nm, dtype=np.float64)
    if box.ndim == 1:
        box = box[:3]
    # pos - floor(pos / L) * L gives [0, L); handles negative (e.g. [-L/2, L/2))
    for d in range(3):
        Ld = box[d]
        if Ld <= 0:
            raise ValueError(f"Box length L[{d}] = {Ld} must be positive")
        pos[:, d] = pos[:, d] - np.floor(pos[:, d] / Ld) * Ld
    return pos


def _load_gsd_config(path, frame_index):
    """
    Load GSD snapshot into a config dict: positions, box (Bohr), bonds_group, bonds_typeid, n_mol, has_cavity.
    Used when building OpenMM system from GSD so topology order matches the file.
    """
    if not HAS_GSD:
        raise ImportError("gsd.hoomd is required; install with: pip install gsd")
    path = Path(path)
    with gsd.hoomd.open(path, "r") as f:
        nframes = len(f)
        if frame_index < 0:
            frame_index = nframes + frame_index
        snap = f[frame_index]
    box = np.array(snap.configuration.box[:3], dtype=np.float64)
    pos = np.array(snap.particles.position, dtype=np.float64)
    bonds_group = np.array(snap.bonds.group, dtype=np.intp)
    bonds_typeid = np.array(snap.bonds.typeid, dtype=np.intp)
    n_mol = len(bonds_group)
    n_mol_particles = 2 * n_mol
    has_cavity = snap.particles.N == n_mol_particles + 1
    if snap.particles.N != n_mol_particles and not has_cavity:
        raise ValueError(
            f"GSD has {snap.particles.N} particles, expected {n_mol_particles} or {n_mol_particles}+1 (with cavity)"
        )
    return {
        "positions_bohr": pos,
        "box_bohr": box,
        "bonds_group": bonds_group,
        "bonds_typeid": bonds_typeid,
        "n_mol": n_mol,
        "has_cavity": has_cavity,
    }


def _unwrap_molecules_bohr(pos_bohr, box_bohr, bonds_group):
    """
    Unwrap positions (Bohr) so bonded atoms are close; then shift so min in [0, L).
    OpenMM HarmonicBondForce uses direct distance; bonded pairs must not span the box.
    Positions may extend past L; do NOT wrap after this (wrapping would break bond lengths).
    """
    pos = np.asarray(pos_bohr, dtype=np.float64).copy()
    box = np.asarray(box_bohr[:3], dtype=np.float64)
    bonds = np.asarray(bonds_group, dtype=np.intp)
    for b in range(len(bonds)):
        i, j = int(bonds[b][0]), int(bonds[b][1])
        dr = pos[j] - pos[i]
        shift = np.round(dr / box) * box
        pos[j] -= shift
    min_pos = np.amin(pos, axis=0)
    pos -= np.floor(min_pos / box) * box
    return pos


def _build_system_from_gsd(gsd_path, frame_index, gsd_in_nm, no_cavity, ff_dir,
                           cutoff_nm=DIAMER_RCUT_NM):
    """
    Build OpenMM system from GSD so particle order and types match the file (avoids blow-up from order mismatch).
    Returns (system, positions, topology, cavity_index or None, num_molecules, box_size_nm).
    Uses PME with cutoff_nm (default 15 Bohr = 0.794 nm, cav-hoomd parity).
    """
    if ff_dir is None:
        ff_dir = Path(__file__).parent
    ff_path = ff_dir / "diamer_forcefield.xml"
    if not ff_path.exists():
        raise FileNotFoundError(f"ForceField XML not found: {ff_path}")
    cfg = _load_gsd_config(gsd_path, frame_index)
    n_mol = cfg["n_mol"]
    has_cavity = cfg["has_cavity"] and not no_cavity
    bonds_group = cfg["bonds_group"]
    bonds_typeid = cfg["bonds_typeid"]
    box_bohr = np.asarray(cfg["box_bohr"][:3], dtype=np.float64)
    scale = 1.0 if gsd_in_nm else BOHR_TO_NM
    box_nm = box_bohr * scale
    # Use raw GSD positions directly (same convention as cav-hoomd: [-L/2, L/2]).
    # HarmonicBondForce.setUsesPeriodicBoundaryConditions(True) handles bonds that
    # span the boundary, so _unwrap_molecules_bohr() is not needed. Using raw positions
    # ensures the cavity dipole matches cav-hoomd from t=0.
    pos_nm = cfg["positions_bohr"] * scale

    topology = Topology()
    chain = topology.addChain()
    positions_list = []
    for b in range(n_mol):
        i, j = int(bonds_group[b][0]), int(bonds_group[b][1])
        is_oo = int(bonds_typeid[b]) == 0
        res_name = "OO" if is_oo else "NN"
        elem = Element.getBySymbol("O") if is_oo else Element.getBySymbol("N")
        res = topology.addResidue(res_name, chain)
        a1 = topology.addAtom("A", elem, res)
        a2 = topology.addAtom("B", elem, res)
        topology.addBond(a1, a2)
        positions_list.append(openmm.Vec3(float(pos_nm[i, 0]), float(pos_nm[i, 1]), float(pos_nm[i, 2])) * unit.nanometer)
        positions_list.append(openmm.Vec3(float(pos_nm[j, 0]), float(pos_nm[j, 1]), float(pos_nm[j, 2])) * unit.nanometer)
    cavity_index = None
    if has_cavity:
        cav_res = topology.addResidue("CAV", chain)
        topology.addAtom("Q", Element.getBySymbol("He"), cav_res)
        cavity_index = 2 * n_mol
        positions_list.append(openmm.Vec3(float(pos_nm[-1, 0]), float(pos_nm[-1, 1]), float(pos_nm[-1, 2])) * unit.nanometer)

    topology.setPeriodicBoxVectors((
        openmm.Vec3(box_nm[0], 0, 0) * unit.nanometer,
        openmm.Vec3(0, box_nm[1], 0) * unit.nanometer,
        openmm.Vec3(0, 0, box_nm[2]) * unit.nanometer,
    ))
    forcefield = ForceField(str(ff_path))
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=PME,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
    )
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_nm[0], 0, 0),
        openmm.Vec3(0, box_nm[1], 0),
        openmm.Vec3(0, 0, box_nm[2]),
    )
    # Enable PBC on HarmonicBondForce so bonds that span the periodic boundary
    # are computed using minimum-image convention (same as cav-hoomd). This lets us
    # use raw GSD positions without unwrapping.
    for f in system.getForces():
        if isinstance(f, openmm.HarmonicBondForce):
            f.setUsesPeriodicBoundaryConditions(True)

    if cavity_index is not None:
        nbf = None
        ljf = None
        for f in system.getForces():
            if isinstance(f, openmm.NonbondedForce):
                nbf = f
            if isinstance(f, openmm.CustomNonbondedForce) and f.getName() == "LennardJones":
                ljf = f
        for i in range(system.getNumParticles()):
            if i == cavity_index:
                continue
            if ljf is not None:
                ljf.addExclusion(cavity_index, i)
            if nbf is not None:
                nbf.addException(cavity_index, i, 0.0, 0.1, 0.0)
    return system, positions_list, topology, cavity_index, n_mol, float(box_nm[0])


def _load_positions_and_box_from_gsd(path, frame_index, gsd_in_nm=False, bohr_to_nm=BOHR_TO_NM):
    """
    Load positions and box from a GSD file (cav-hoomd or HOOMD format).

    Units:
      - If gsd_in_nm is False (default): assumes GSD is in atomic units (Bohr);
        multiplies positions and box by bohr_to_nm to get nm.
      - If gsd_in_nm is True: assumes GSD is already in nm (no conversion).

    Positions are always wrapped into [0, L) per dimension so OpenMM's primary
    image is consistent (handles GSD in [0,L) or [-L/2, L/2) convention).

    Returns (positions_list, box_nm) where positions_list is list of
    openmm.Vec3 * unit.nanometer, box_nm is (Lx, Ly, Lz) in nm.
    """
    if not HAS_GSD:
        raise ImportError("gsd.hoomd is required for --initial-gsd; install with: pip install gsd")
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Initial GSD not found: {path}")
    with gsd.hoomd.open(path, "r") as f:
        nframes = len(f)
        if frame_index < 0:
            frame_index = nframes + frame_index
        if frame_index < 0 or frame_index >= nframes:
            raise ValueError(f"Frame index {frame_index} out of range for GSD with {nframes} frames")
        frame = f[frame_index]
    pos = np.array(frame.particles.position, dtype=np.float64)
    box_raw = np.array(frame.configuration.box[:3], dtype=np.float64)

    if gsd_in_nm:
        pos_nm = np.asarray(pos, dtype=np.float64)
        box_nm = np.asarray(box_raw, dtype=np.float64)
    else:
        # GSD from cav-hoomd / initlattice is in Bohr
        pos_nm = pos * bohr_to_nm
        box_nm = box_raw * bohr_to_nm

    # Wrap into [0, L) so OpenMM primary image is consistent
    pos_nm = _wrap_positions_into_box(pos_nm, box_nm)

    positions_list = [openmm.Vec3(float(x), float(y), float(z)) * unit.nanometer for x, y, z in pos_nm]
    return positions_list, box_nm, box_raw


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
    nonbonded_force.setNonbondedMethod(openmm.NonbondedForce.PME)
    nonbonded_force.setCutoffDistance(DIAMER_RCUT_NM)  # 15 Bohr, cav-hoomd parity
    
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
        nonbondedMethod=PME,
        nonbondedCutoff=DIAMER_RCUT_NM * unit.nanometer,
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


def _topology_subset(topology, num_atoms):
    """
    Return a new Topology containing only the first num_atoms atoms.
    Used for PDB output when cavity is in the full topology but we write molecular atoms only.
    """
    from openmm.app import Topology
    new_top = Topology()
    new_top.setPeriodicBoxVectors(topology.getPeriodicBoxVectors())
    count = 0
    for chain in topology.chains():
        if count >= num_atoms:
            break
        new_chain = new_top.addChain(getattr(chain, "id", None))
        for residue in chain.residues():
            if count >= num_atoms:
                break
            new_res = new_top.addResidue(
                residue.name, new_chain,
                getattr(residue, "id", None),
                getattr(residue, "insertionCode", "") or "",
            )
            for atom in residue.atoms():
                if count >= num_atoms:
                    break
                new_top.addAtom(
                    atom.name, atom.element, new_res,
                    getattr(atom, "id", None),
                )
                count += 1
    return new_top


def run_test(num_molecules=250, lambda_coupling=0.001, temperature_K=100.0,
             dt=0.001, equilibration_time_ps=100.0, production_time_ps=900.0,
             cavity_freq_cm=1560.0, disable_dipole_output=False,
             box_size_nm=2.5, fraction_OO=0.8, friction=0.01, minimize=True,
             output_file=None, seed=42, report_interval_steps=1000,
             pdb_file=None, pdb_interval_steps=1000,
             enable_fkt=False, fkt_kmag=113.4, fkt_num_wavevectors=50,
             fkt_reference_interval_ps=1.0, fkt_max_refs=10,
             fkt_output_period_ps=1.0, fkt_output_prefix=None,
             nve=False, bussi_tau_ps=5.0,
             initial_gsd=None, initial_gsd_frame=-1, gsd_in_nm=False, no_cavity=False,
             cutoff_nm=DIAMER_RCUT_NM):
    """Run the cavity diamer simulation test.
    If nve=True: NVE (microcanonical) with Verlet; initial velocities thermalized at temperature_K.
    If nve=False: Verlet + Bussi (molecules only), cav-hoomd parity; default Bussi tau 5.0 ps.
    If no_cavity=True: 500 particles only (bare glassy system), no cavity particle or CavityForce.
    cutoff_nm: LJ + PME real-space cutoff (default 15 Bohr = 0.794 nm, cav-hoomd parity).
    """
    
    print("=" * 60)
    print("Cavity OpenMM Diamer Test")
    print("=" * 60)

    # Load initial structure from GSD if requested (same initial conditions as cav-hoomd).
    # When initial_gsd is set we build the OpenMM system FROM the GSD topology so particle order
    # and types (O-O vs N-N) match the file; otherwise positions-by-index would mismatch and blow up.
    build_from_gsd = bool(initial_gsd)
    gsd_positions = None
    gsd_box_nm = None
    if build_from_gsd:
        if not HAS_GSD:
            raise ImportError("gsd.hoomd is required for --initial-gsd; install with: pip install gsd")
        ff_dir = Path(__file__).parent
        system, positions, topology, cavity_index, num_molecules, box_size_nm = _build_system_from_gsd(
            initial_gsd, initial_gsd_frame, gsd_in_nm, no_cavity, ff_dir,
            cutoff_nm=cutoff_nm
        )
        n_sys = system.getNumParticles()
        print(f"\n--- Initial structure from GSD ---")
        print(f"  File: {initial_gsd}, frame: {initial_gsd_frame}")
        print(f"  Box from GSD: {box_size_nm:.4f} nm (cubic)")
        print(f"  Positions: {n_sys} particles (topology from GSD; order matches file); energy minimization will be skipped")
        if box_size_nm < 0.3:
            print(f"  *** WARNING: Box {box_size_nm:.4f} nm is very small. If GSD is in nm (not Bohr), use --gsd-in-nm. ***")
        elif box_size_nm > 50.0:
            print(f"  *** WARNING: Box {box_size_nm:.4f} nm is very large. If GSD is in Bohr, do not use --gsd-in-nm. ***")
        minimize = False
    
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
    if nve:
        print(f"  Mode: NVE (microcanonical); initial velocities thermalized at T = {temperature_K} K")
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
        print(f"    wrapped molecular positions (cav-hoomd parity; correct at 0 K)")
        print(f"    reference interval (new file every): {fkt_reference_interval_ps} ps")
        print(f"    max reference files (FIFO drop): {fkt_max_refs}")
        print(f"    F(k,t) output period: {fkt_output_period_ps} ps")
    
    # Create system: from GSD (topology matches file) or from forcefield (fixed order).
    if build_from_gsd:
        # system, positions, topology, cavity_index, num_molecules, box_size_nm already set above
        cavity_available = cavity_index is not None
        num_molecular_particles = cavity_index if cavity_index is not None else system.getNumParticles()
        print("\n--- Creating Diamer System ---")
        print(f"Created system with {system.getNumParticles()} particles ({num_molecules} molecules)")
        print(f"  Box size: {box_size_nm:.4f} nm")
        print(f"  Topology from GSD (particle order matches file)")
        if cavity_available:
            if lambda_coupling == 0:
                system.setParticleMass(cavity_index, 1e6)
                print(f"  Cavity particle at index {cavity_index}, mass set to 1e6 amu (fixed for g=0)")
            else:
                system.setParticleMass(cavity_index, photon_mass)
                print(f"  Cavity particle at index {cavity_index}, mass set to {photon_mass:.6f} amu")
        else:
            print("  No cavity (bare glassy system, 500 particles)")
    else:
        print("\n--- Creating Diamer System ---")
        include_cavity = not no_cavity
        result = create_diamer_system_from_forcefield(
            num_molecules=num_molecules,
            fraction_OO=fraction_OO,
            box_size_nm=box_size_nm,
            seed=seed,
            include_cavity=include_cavity,
        )
        cavity_available = False
        cavity_index = None
        if len(result) == 4:
            system, positions, topology, cavity_index = result
            num_molecular_particles = cavity_index
            cavity_available = True
            if lambda_coupling == 0:
                system.setParticleMass(cavity_index, 1e6)
                print("\n--- Cavity (from topology) ---")
                print(f"  Cavity particle at index {cavity_index}, mass set to 1e6 amu (fixed for g=0)")
            else:
                system.setParticleMass(cavity_index, photon_mass)
                print("\n--- Cavity (from topology) ---")
                print(f"  Cavity particle at index {cavity_index}, mass set to {photon_mass:.6f} amu")
        else:
            system, positions, topology = result
            num_molecular_particles = system.getNumParticles()
            if not no_cavity:
                print("\n--- Adding Cavity Particle ---")
                cavity_index = add_cavity_particle(system, positions, omegac_au, photon_mass)
                num_molecular_particles = cavity_index
                cavity_available = True
                if lambda_coupling == 0:
                    system.setParticleMass(cavity_index, 1e6)
                    print(f"  Cavity mass set to 1e6 amu (fixed for g=0)")
            else:
                print("  No cavity (bare glassy system, 500 particles)")

        if gsd_positions is not None:
            n_sys = system.getNumParticles()
            if len(gsd_positions) > n_sys:
                gsd_positions = gsd_positions[:n_sys]
                print(f"  Using first {n_sys} particles from GSD (no-cavity mode)")
            if len(gsd_positions) != n_sys:
                raise ValueError(
                    f"GSD has {len(gsd_positions)} particles but system has {n_sys}; "
                    "need at least system.getNumParticles() (500 for no-cavity, 501 with cavity)"
                )
            positions = gsd_positions
            topology.setPeriodicBoxVectors((
                openmm.Vec3(gsd_box_nm[0], 0, 0) * unit.nanometer,
                openmm.Vec3(0, gsd_box_nm[1], 0) * unit.nanometer,
                openmm.Vec3(0, 0, gsd_box_nm[2]) * unit.nanometer,
            ))
            system.setDefaultPeriodicBoxVectors(
                openmm.Vec3(gsd_box_nm[0], 0, 0),
                openmm.Vec3(0, gsd_box_nm[1], 0),
                openmm.Vec3(0, 0, gsd_box_nm[2]),
            )
            print("  Replaced positions and box with GSD values")
    
    # Store charges for dipole calculation
    charges = []
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for i in range(system.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(i)
                charges.append(charge)
    
    # Add CavityForce only when cavity is present (not --no-cavity)
    cavity_force = None
    if cavity_available:
        print("\n--- Adding Cavity Force ---")
        try:
            # Create CavityForce with target lambda from t=0 (omegac in atomic units)
            cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass)
            system.addForce(cavity_force)
            print(f"  CavityForce added successfully")
            print(f"  Omega_c: {omegac_au:.6f} a.u.")
            print(f"  Photon mass: {photon_mass:.6f} amu = {photon_mass * 1822.888:.1f} a.u.")
            AMU_TO_AU = 1822.888
            photon_mass_au = photon_mass * AMU_TO_AU
            K_au = photon_mass_au * omegac_au**2
            K_openmm = K_au * HARTREE_TO_KJMOL / (BOHR_TO_NM**2)
            print(f"  Expected spring constant K: {K_openmm:.0f} kJ/(mol·nm^2)")
            print(f"  Lambda coupling: {lambda_coupling} (ACTIVE from t=0)")
            displacer = openmm.CavityParticleDisplacer(cavity_index, omegac_au, photon_mass)
            displacer.setSwitchOnLambda(lambda_coupling)
            system.addForce(displacer)
            print(f"  CavityParticleDisplacer added (finite-Q mode)")
        except AttributeError as e:
            print(f"  CavityForce not available: {e}")
            print("  Skipping cavity coupling test")
            cavity_available = False
            cavity_force = None
    else:
        print("\n--- Cavity Force ---")
        print("  Skipped (--no-cavity: bare glassy system)")
    
    # Thermostat: Bussi only when not NVE (cav-hoomd parity: molecules = ConstantVolume + Bussi)
    bussi_available = False
    if not nve:
        print("\n--- Adding Bussi Thermostat ---")
        try:
            bussi = openmm.BussiThermostat(temperature_K, bussi_tau_ps)
            bussi.setApplyToAllParticles(False)
            for i in range(num_molecular_particles):
                bussi.addParticle(i)
            system.addForce(bussi)
            print(f"  BussiThermostat added for {bussi.getNumParticles()} particles")
            print(f"  Temperature: {temperature_K} K")
            print(f"  Tau: {bussi_tau_ps} ps (cav-hoomd parity)")
            bussi_available = True
        except AttributeError as e:
            print(f"  BussiThermostat not available: {e}")

    # Integrator: Verlet only (NVE) or Verlet + Bussi (no Langevin; molecules get Bussi only)
    print("\n--- Creating Integrator ---")
    if nve:
        integrator = openmm.VerletIntegrator(dt * unit.picosecond)
        print(f"  VerletIntegrator created (NVE; no thermostat)")
        print(f"  Initial velocities will be set to T = {temperature_K} K")
    else:
        integrator = openmm.VerletIntegrator(dt * unit.picosecond)
        print(f"  VerletIntegrator created")
        print(f"  Verlet + Bussi (molecules only); Bussi tau {bussi_tau_ps} ps (cav-hoomd parity)")
    
    # Create simulation - try CUDA first, fall back to Reference (CPU)
    print("\n--- Creating Simulation ---")
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        print(f"  Using CUDA platform (GPU acceleration)")
    except Exception:
        platform = openmm.Platform.getPlatformByName('Reference')
        print(f"  Using Reference platform (CPU; CUDA not available)")
    
    # Create context
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)

    # Phase 1 debug: log initial PE/KE when built from GSD (pinpoint when blow-up occurs)
    if build_from_gsd:
        state = context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        print("\n--- [DEBUG] Initial state (after setPositions) ---")
        print(f"  PE: {pe:.6g} kJ/mol  KE: {ke:.6g} kJ/mol")
        if abs(pe) > 1e15:
            print("  *** Initial PE is huge; problem is likely initial configuration or force field. ***")
    
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
    
    # Set velocities (thermalize at target T; in NVE this is the only thermalization)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    # Cavity particle has negligible mass (0.000549 amu); thermalizing it at 100 K gives it a huge
    # v_rms = sqrt(kT/m) that makes the stiff cavity spring numerically unstable with 1 fs dt.
    # Start the cavity from rest so the initial configuration (e.g. from GSD) is not destroyed.
    if cavity_available:
        state = context.getState(getVelocities=True)
        vels = list(state.getVelocities())
        vels[cavity_index] = openmm.Vec3(0, 0, 0) * (unit.nanometer / unit.picosecond)
        context.setVelocities(vels)

    # Phase 1 debug: log PE/KE after thermalization when built from GSD
    if build_from_gsd:
        state = context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        print("\n--- [DEBUG] After thermalization (before first step) ---")
        print(f"  PE: {pe:.6g} kJ/mol  KE: {ke:.6g} kJ/mol")
    
    # Calculate total steps
    total_steps = equilibration_steps + production_steps
    
    print("\n--- Running Simulation ---")
    print(f"  Equilibration: {equilibration_time_ps} ps ({equilibration_steps} steps)")
    print(f"  Production:    {production_time_ps} ps ({production_steps} steps)")
    print(f"  Total:         {equilibration_time_ps + production_time_ps} ps ({total_steps} steps)")
    if enable_fkt:
        print(f"  F(k,t) tracking starts after equilibration (production t=0)")
    print("")
    
    # Prepare dipole trajectory storage
    dipole_times = []
    dipole_trajectory = []
    step_counter = 0  # Manual step counter
    
    report_interval = report_interval_steps
    num_reports = -(-total_steps // report_interval)  # ceiling division so no steps are dropped
    
    # PDB output (molecular atoms only; exclude cavity)
    pdb_handle = None
    pdb_model_index = 0
    topology_molecular = None  # topology with num_molecular_particles atoms for PDB
    if pdb_file:
        pdb_handle = open(pdb_file, 'w')
        topology_molecular = _topology_subset(topology, num_molecular_particles)
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
        # F(k,t) uses wrapped positions (cav-hoomd parity; correct at 0 K, no unwrapping drift)
    
    first_step_pe_ke_logged = False
    for i in range(num_reports):
        steps_this_block = min(report_interval, total_steps - step_counter)
        for step in range(steps_this_block):
            integrator.step(1)
            step_counter += 1
            # Phase 1 debug: log PE/KE once after first step when built from GSD
            if build_from_gsd and not first_step_pe_ke_logged and step_counter >= 1:
                state = context.getState(getEnergy=True)
                pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                print("\n--- [DEBUG] After first integrator step ---")
                print(f"  PE: {pe:.6g} kJ/mol  KE: {ke:.6g} kJ/mol  (step_counter={step_counter})")
                first_step_pe_ke_logged = True
            if not disable_dipole_output and (step_counter % dipole_output_interval_steps) == 0:
                state = context.getState(getPositions=True)
                current_time = step_counter * dt
                dipole = compute_dipole_moment(state, charges, num_molecular_particles)
                dipole_times.append(current_time)
                dipole_trajectory.append(dipole)
            # PDB frame output (molecular atoms only; PDBFile expects angstroms)
            if pdb_handle and topology_molecular is not None and (step_counter % pdb_interval_steps) == 0:
                _state = context.getState(getPositions=True)
                _pos_all = _state.getPositions(asNumpy=True)
                _pos_mol = _pos_all[:num_molecular_particles].value_in_unit(unit.angstroms)
                PDBFile.writeModel(topology_molecular, _pos_mol, pdb_handle, modelIndex=pdb_model_index + 1)
                pdb_model_index += 1
            # F(k,t): use wrapped molecular positions (enforcePeriodicBox=True ensures
            # particles are in the primary image, matching HOOMD's snap.particles.position).
            # Without this flag, OpenMM returns unwrapped positions where boundary crossings
            # add spurious phase shifts exp(i k·n L) that decorrelate ρ_k artificially.
            # Skip during equilibration: freshly randomised velocities need ~10 tau_Bussi
            # to reach the correct position-velocity correlations; collecting F(k,t) during
            # this transient biases ref 0 and contaminates the grand average.
            if enable_fkt and fkt_tracker is not None and (step_counter % fkt_interval_steps) == 0:
                current_time_ps = step_counter * dt
                if current_time_ps >= equilibration_time_ps:
                    _st = context.getState(getPositions=True, enforcePeriodicBox=True)
                    _pos_all = _st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                    pos_mol = np.asarray(_pos_all[:num_molecular_particles], dtype=np.float64)
                    fkt_time_ps = current_time_ps - equilibration_time_ps
                    fkt_tracker.update(fkt_time_ps, pos_mol)
        
        # Report progress with timing
        current_wall_time = time.time()
        elapsed = current_wall_time - start_time
        elapsed_since_last = current_wall_time - last_report_time
        last_report_time = current_wall_time
        
        state = context.getState(getEnergy=True, getPositions=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        n_dof = 3 * system.getNumParticles() - system.getNumConstraints()
        k_B_kJ = 0.008314462618  # kJ/(mol·K)
        T_kin = (2.0 * ke / (n_dof * k_B_kJ)) if n_dof > 0 else 0.0
        positions_now = state.getPositions()
        cx = cy = cz = 0.0
        if cavity_available and cavity_index is not None:
            cavity_pos = positions_now[cavity_index]
            cx = cavity_pos[0].value_in_unit(unit.nanometer) if hasattr(cavity_pos[0], 'value_in_unit') else float(cavity_pos[0])
            cy = cavity_pos[1].value_in_unit(unit.nanometer) if hasattr(cavity_pos[1], 'value_in_unit') else float(cavity_pos[1])
            cz = cavity_pos[2].value_in_unit(unit.nanometer) if hasattr(cavity_pos[2], 'value_in_unit') else float(cavity_pos[2])
        
        sim_time_ps = step_counter * dt
        progress_pct = (step_counter / total_steps) * 100
        steps_per_sec = steps_this_block / elapsed_since_last
        
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
        if cavity_available and cavity_force is not None:
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
        
        print(f"         PE: {pe:8.2f} kJ/mol | T_kin: {T_kin:6.1f} K{cavity_energy_str}")
        if cavity_available:
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
                         'cavity_index': cavity_index if cavity_available else None,
                         'status': 'running',
                         'step': step_counter,
                         'total_steps': total_steps
                     })
    
    # Close PDB trajectory
    if pdb_handle:
        PDBFile.writeFooter(topology_molecular if topology_molecular is not None else topology, pdb_handle)
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
                    'cavity_index': cavity_index if cavity_available else None,
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

  # NVE: no thermostat; initial velocities thermalized at --temp, then constant energy
  python run_simulation.py --nve --dimers 50 --equil 0 --prod 100
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
    parser.add_argument('--nve', action='store_true',
                        help='NVE (microcanonical): Verlet integrator, no thermostat; initial velocities thermalized at --temp')
    parser.add_argument('--bussi-tau', type=float, default=5.0, dest='bussi_tau',
                        help='Bussi thermostat time constant in ps when not --nve (default: 5.0, cav-hoomd parity)')
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
    parser.add_argument('--initial-gsd', type=str, default=None, dest='initial_gsd',
                        help='Load initial positions and box from this GSD file (Bohr->nm); skip minimization')
    parser.add_argument('--initial-gsd-frame', type=int, default=-1, dest='initial_gsd_frame',
                        help='Frame index to read from --initial-gsd (default: -1, last frame)')
    parser.add_argument('--gsd-in-nm', action='store_true', dest='gsd_in_nm',
                        help='GSD positions and box are already in nm (do not convert from Bohr). Use if init-0.gsd was written in nm.')
    parser.add_argument('--no-cavity', action='store_true', dest='no_cavity',
                        help='Run bare glassy system only (500 particles, no cavity); matches cav-hoomd --no-cavity')
    parser.add_argument('--cutoff', type=float, default=DIAMER_RCUT_NM, dest='cutoff_nm',
                        help=f'LJ + PME real-space cutoff in nm (default: {DIAMER_RCUT_NM:.4f} = 15 Bohr, cav-hoomd parity)')
    
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
            nve=args.nve,
            bussi_tau_ps=args.bussi_tau,
            initial_gsd=args.initial_gsd,
            initial_gsd_frame=args.initial_gsd_frame,
            gsd_in_nm=args.gsd_in_nm,
            no_cavity=args.no_cavity,
            cutoff_nm=args.cutoff_nm,
        )
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
