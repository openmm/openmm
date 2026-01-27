#!/usr/bin/env python3
"""
Cavity Molecular Dynamics Simulation of Water with UMA
========================================================

System: Water molecules using FAIRChem UMA potential
Cavity: Optional resonant coupling with photon mode
Protocol: Equilibrate → Production with dipole recording

The UMA ML potential handles molecular interactions, while the cavity
coupling is added on top via CavityForce and CavityParticleDisplacer.
"""

import sys
import numpy as np
import time
from pathlib import Path

try:
    import openmm
    from openmm import unit
    from openmm.app import Topology, Element
    from openmmml import MLPotential
    import torch
    print("✓ OpenMM and OpenMM-ML loaded successfully")
    print(f"  PyTorch version: {torch.__version__}")
    print(f"  CUDA available: {torch.cuda.is_available()}")
except ImportError as e:
    print(f"Error importing required packages: {e}")
    sys.exit(1)

# Unit conversions
BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5
HARTREE_TO_CM = 219474.63
AMU_TO_AU = 1822.888


def wavenumber_to_hartree(omega_cm):
    """Convert wavenumber (cm⁻¹) to Hartree (atomic units)."""
    return omega_cm / HARTREE_TO_CM


def create_water_topology(num_molecules=10):
    """Create a water box topology for UMA simulation."""
    print(f"\n--- Creating Water Topology ---")
    print(f"  Model: UMA ML Potential")
    print(f"  Molecules: {num_molecules}")
    
    oh_bond = 0.09572  # nm
    hoh_angle = 104.52 * np.pi / 180.0
    
    # Calculate box size (density ~1 g/cm³)
    vol_per_mol_nm3 = 0.0299
    total_vol = num_molecules * vol_per_mol_nm3
    box_size_nm = max((total_vol ** (1/3)) * 1.2, 1.0)
    
    print(f"  Box size: {box_size_nm:.2f} nm")
    
    topology = Topology()
    positions = []
    
    side = int(np.ceil(num_molecules**(1/3)))
    spacing = box_size_nm / side
    
    mol_count = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if mol_count >= num_molecules:
                    break
                    
                chain = topology.addChain()
                residue = topology.addResidue('HOH', chain)
                topology.addAtom('O', Element.getBySymbol('O'), residue)
                topology.addAtom('H', Element.getBySymbol('H'), residue)
                topology.addAtom('H', Element.getBySymbol('H'), residue)
                
                cx = (i + 0.5) * spacing
                cy = (j + 0.5) * spacing
                cz = (k + 0.5) * spacing
                
                np.random.seed(42 + mol_count)
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                psi = np.random.rand() * 2 * np.pi
                
                h1_local = np.array([oh_bond * np.sin(hoh_angle/2), 0, oh_bond * np.cos(hoh_angle/2)])
                h2_local = np.array([-oh_bond * np.sin(hoh_angle/2), 0, oh_bond * np.cos(hoh_angle/2)])
                o_local = np.array([0, 0, 0])
                
                Rz = np.array([[np.cos(theta), -np.sin(theta), 0],
                               [np.sin(theta), np.cos(theta), 0],
                               [0, 0, 1]])
                Ry = np.array([[np.cos(phi), 0, np.sin(phi)],
                               [0, 1, 0],
                               [-np.sin(phi), 0, np.cos(phi)]])
                Rx = np.array([[1, 0, 0],
                               [0, np.cos(psi), -np.sin(psi)],
                               [0, np.sin(psi), np.cos(psi)]])
                R = Rz @ Ry @ Rx
                
                positions.append((R @ o_local + np.array([cx, cy, cz])).tolist())
                positions.append((R @ h1_local + np.array([cx, cy, cz])).tolist())
                positions.append((R @ h2_local + np.array([cx, cy, cz])).tolist())
                
                mol_count += 1
                
            if mol_count >= num_molecules:
                break
        if mol_count >= num_molecules:
            break
    
    # Get charges for dipole calculation (TIP3P-like)
    charges = []
    for atom in topology.atoms():
        if atom.element.symbol == 'O':
            charges.append(-0.8476)
        else:
            charges.append(0.4238)
    
    print(f"  Created {mol_count} water molecules ({len(positions)} atoms)")
    
    return topology, positions, box_size_nm, charges


def add_cavity_particle(system, positions, omegac_au, photon_mass_amu):
    """Add cavity photon mode as a dummy particle."""
    print(f"\n--- Adding Cavity Particle ---")
    
    cavity_index = system.addParticle(photon_mass_amu * unit.amu)
    positions.append([0.0, 0.0, 0.0])  # Cavity at origin
    
    omegac_cm = omegac_au * HARTREE_TO_CM
    print(f"  Cavity particle index: {cavity_index}")
    print(f"  Omega_c: {omegac_au:.6f} a.u. = {omegac_cm:.1f} cm⁻¹")
    print(f"  Photon mass: {photon_mass_amu:.6f} amu")
    
    return cavity_index


def add_charges_for_cavity(system, topology, cavity_index):
    """
    Add a NonbondedForce with charges for cavity dipole calculation.
    
    IMPORTANT: UMA handles all molecular interactions via PythonForce, but
    CavityForce needs atomic charges to compute the dipole moment.
    We add a NonbondedForce with:
    - Proper charges (TIP3P-like for water)
    - Zero LJ parameters (epsilon=0) so it adds no nonbonded energy
    - NoCutoff method (just for charge storage, no actual interactions)
    """
    print(f"\n--- Adding Charges for Cavity Coupling ---")
    
    # Create NonbondedForce for charges only
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    
    # Add charges for molecular atoms
    # Using TIP3P charges: O=-0.834, H=+0.417
    n_molecular = 0
    for atom in topology.atoms():
        if atom.element.symbol == 'O':
            charge = -0.834
        elif atom.element.symbol == 'H':
            charge = 0.417
        else:
            charge = 0.0
        # Zero LJ (sigma=0.1nm, epsilon=0) - no actual LJ interaction
        nb.addParticle(charge, 0.1, 0.0)
        n_molecular += 1
    
    # Add zero charge for cavity particle
    nb.addParticle(0.0, 0.1, 0.0)
    
    # Disable all pairwise interactions by adding exceptions
    # (we only want charges for dipole, not actual Coulomb interactions)
    n_total = nb.getNumParticles()
    for i in range(n_total):
        for j in range(i+1, n_total):
            nb.addException(i, j, 0.0, 0.1, 0.0)  # Zero charge product, zero LJ
    
    system.addForce(nb)
    print(f"  ✓ NonbondedForce added with {n_molecular} molecular charges")
    print(f"  Note: Charges used for dipole only (all interactions from UMA)")
    
    return nb


def setup_cavity_coupling(system, topology, cavity_index, omegac_au, lambda_coupling, photon_mass_amu):
    """Set up cavity coupling forces."""
    print(f"\n--- Setting Up Cavity Coupling ---")
    print(f"  Lambda: {lambda_coupling}")
    
    # First add charges so CavityForce can compute dipole
    add_charges_for_cavity(system, topology, cavity_index)
    
    # CavityForce
    cavity_force = openmm.CavityForce(cavity_index, omegac_au, lambda_coupling, photon_mass_amu)
    system.addForce(cavity_force)
    print(f"  ✓ CavityForce added")
    
    # CavityParticleDisplacer for finite-q correction
    displacer = openmm.CavityParticleDisplacer(cavity_index, omegac_au, photon_mass_amu)
    displacer.setSwitchOnStep(0)
    displacer.setSwitchOnLambda(lambda_coupling)
    system.addForce(displacer)
    print(f"  ✓ CavityParticleDisplacer added")
    
    return cavity_force, displacer


def compute_dipole(positions, charges, num_molecular):
    """Compute total molecular dipole moment (excluding cavity particle)."""
    dipole = np.zeros(3)
    for i in range(num_molecular):
        dipole += charges[i] * np.array(positions[i])
    return dipole


def run_simulation(lambda_coupling=0.0, cavity_freq_cm=3500.0,
                  num_molecules=10, temperature_K=300.0, 
                  dt_ps=0.0005, equilibration_ps=10.0, production_ps=100.0,
                  output_prefix="uma_water_cavity"):
    """
    Run UMA water MD simulation with optional cavity coupling.
    
    Parameters
    ----------
    lambda_coupling : float
        Coupling strength (0.0 = no cavity)
    cavity_freq_cm : float
        Cavity frequency in cm⁻¹ (default: 3500 for OH stretch)
    num_molecules : int
        Number of water molecules
    temperature_K : float
        Temperature in K
    dt_ps : float
        Timestep in ps
    equilibration_ps : float
        Equilibration time in ps
    production_ps : float
        Production time in ps
    output_prefix : str
        Prefix for output files
    """
    print("=" * 80)
    print("UMA Water Cavity Molecular Dynamics Simulation")
    print("=" * 80)
    
    # Cavity parameters
    omegac_au = wavenumber_to_hartree(cavity_freq_cm)
    photon_mass_amu = 1.0 / AMU_TO_AU
    
    mode_name = "OH stretch" if cavity_freq_cm > 2000 else "HOH bending"
    use_cavity = lambda_coupling > 0
    
    print(f"\nSimulation Parameters:")
    print(f"  Water molecules: {num_molecules}")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt_ps} ps ({dt_ps*1000:.2f} fs)")
    print(f"  Equilibration: {equilibration_ps} ps")
    print(f"  Production: {production_ps} ps")
    if use_cavity:
        print(f"  Cavity: ON")
        print(f"    Frequency: {cavity_freq_cm:.0f} cm⁻¹ ({mode_name})")
        print(f"    Lambda: {lambda_coupling}")
    else:
        print(f"  Cavity: OFF (baseline simulation)")
    
    # Calculate steps
    equilibration_steps = int(equilibration_ps / dt_ps)
    production_steps = int(production_ps / dt_ps)
    total_steps = equilibration_steps + production_steps
    
    # Create system
    np.random.seed(42)
    topology, positions, box_size_nm, charges = create_water_topology(num_molecules)
    num_molecular = len(positions)  # Before adding cavity
    
    # Create UMA potential
    print(f"\n--- Creating UMA Potential ---")
    print(f"  Loading uma-s-1p1 model...")
    potential = MLPotential('uma-s-1p1-pythonforce')
    
    # Create system
    print(f"\n--- Creating OpenMM System ---")
    system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
    print(f"  ✓ System created with {system.getNumParticles()} particles")
    
    # Add cavity if coupling > 0
    cavity_index = None
    if use_cavity:
        cavity_index = add_cavity_particle(system, positions, omegac_au, photon_mass_amu)
        cavity_force, displacer = setup_cavity_coupling(
            system, topology, cavity_index, omegac_au, lambda_coupling, photon_mass_amu
        )
    
    # Create integrator
    print(f"\n--- Creating Integrator ---")
    friction = 1.0  # ps⁻¹
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt_ps * unit.picosecond
    )
    print(f"  LangevinMiddleIntegrator: friction = {friction} ps⁻¹, dt = {dt_ps} ps")
    
    # Create context
    print(f"\n--- Creating Simulation Context ---")
    platform = openmm.Platform.getPlatformByName('CUDA')
    print(f"  Using CUDA platform")
    
    context = openmm.Context(system, integrator, platform)
    
    # Convert positions to proper format
    pos_with_units = [openmm.Vec3(p[0], p[1], p[2]) * unit.nanometer for p in positions]
    context.setPositions(pos_with_units)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    
    state = context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"  Initial PE: {pe:.1f} kJ/mol")
    
    # Storage
    dipole_times = []
    dipole_values = []
    cavity_values = []
    
    output_file = f"{output_prefix}_lambda{lambda_coupling:.4f}.npz"
    
    print(f"\n--- Running Simulation ---")
    print(f"  Total steps: {total_steps:,}")
    print(f"  Equilibration: {equilibration_steps:,} steps")
    print(f"  Production: {production_steps:,} steps")
    print(f"  Output file: {output_file}")
    
    start_time = time.time()
    report_interval = max(100, total_steps // 100)
    
    for step in range(total_steps):
        integrator.step(1)
        
        # Record data during production
        if step >= equilibration_steps:
            state = context.getState(getPositions=True)
            pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            dipole = compute_dipole(pos, charges, num_molecular)
            
            dipole_times.append((step - equilibration_steps) * dt_ps)
            dipole_values.append(dipole.copy())
            
            if cavity_index is not None:
                q_cav = pos[cavity_index]
                cavity_values.append(q_cav.copy())
        
        # Progress report
        if (step + 1) % report_interval == 0 or (step + 1) == total_steps:
            phase = "EQUIL" if step < equilibration_steps else "PROD"
            pct = 100 * (step + 1) / total_steps
            elapsed = time.time() - start_time
            rate = (step + 1) / max(elapsed, 1e-8)
            eta = (total_steps - step - 1) / max(rate, 1e-8)
            ns_per_day = (rate * dt_ps * 86400) / 1000
            
            state = context.getState(getEnergy=True, getPositions=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            dipole_now = compute_dipole(pos, charges, num_molecular)
            
            # Build status line
            status = f"[{pct:5.1f}%] {phase} | PE: {pe:10.1f} kJ/mol | "
            status += f"d_xy: ({dipole_now[0]:6.2f}, {dipole_now[1]:6.2f}) e·nm | "
            if cavity_index is not None:
                q_cav = pos[cavity_index]
                status += f"q_cav: ({q_cav[0]:6.3f}, {q_cav[1]:6.3f}) nm | "
            status += f"Speed: {ns_per_day:5.3f} ns/day | ETA: {eta/60:5.1f}m"
            print(status, flush=True)
            
            # Save during production
            if phase == "PROD" and len(dipole_times) > 0:
                save_data = {
                    'time_ps': np.array(dipole_times),
                    'dipole_nm': np.array(dipole_values),
                    'metadata': {
                        'lambda_coupling': lambda_coupling,
                        'cavity_freq_cm': cavity_freq_cm,
                        'omega_c_au': omegac_au,
                        'dt_ps': dt_ps,
                        'num_molecules': num_molecules,
                        'box_size_nm': box_size_nm,
                        'temperature_K': temperature_K,
                        'model': 'uma-s-1p1',
                        'status': 'running',
                        'step': step + 1,
                        'total_steps': total_steps
                    }
                }
                if cavity_values:
                    save_data['cavity_nm'] = np.array(cavity_values)
                np.savez(output_file, **save_data)
    
    # Final save
    elapsed = time.time() - start_time
    print(f"\n--- Simulation Complete ---")
    print(f"  Total time: {elapsed/60:.1f} min")
    print(f"  Performance: {(step+1)/elapsed:.1f} steps/s")
    
    save_data = {
        'time_ps': np.array(dipole_times),
        'dipole_nm': np.array(dipole_values),
        'metadata': {
            'lambda_coupling': lambda_coupling,
            'cavity_freq_cm': cavity_freq_cm,
            'omega_c_au': omegac_au,
            'dt_ps': dt_ps,
            'num_molecules': num_molecules,
            'box_size_nm': box_size_nm,
            'temperature_K': temperature_K,
            'model': 'uma-s-1p1',
            'status': 'complete',
            'elapsed_time_s': elapsed
        }
    }
    if cavity_values:
        save_data['cavity_nm'] = np.array(cavity_values)
    np.savez(output_file, **save_data)
    
    print(f"  ✓ Saved: {output_file}")
    print(f"  Data points: {len(dipole_times):,}")


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(description='UMA Water Cavity MD Simulation')
    parser.add_argument('--lambda', type=float, default=0.0, dest='lambda_coupling',
                       help='Coupling strength (default: 0.0 = no cavity)')
    parser.add_argument('--cavity-freq', type=float, default=3500.0, dest='cavity_freq',
                       help='Cavity frequency in cm⁻¹ (default: 3500 for OH stretch)')
    parser.add_argument('--molecules', type=int, default=10,
                       help='Number of water molecules (default: 10)')
    parser.add_argument('--temp', type=float, default=300.0,
                       help='Temperature in K (default: 300)')
    parser.add_argument('--dt', type=float, default=0.0005,
                       help='Timestep in ps (default: 0.0005 = 0.5 fs)')
    parser.add_argument('--equil', type=float, default=10.0,
                       help='Equilibration time in ps (default: 10)')
    parser.add_argument('--prod', type=float, default=100.0,
                       help='Production time in ps (default: 100)')
    parser.add_argument('--output', type=str, default='uma_water_cavity',
                       help='Output prefix (default: uma_water_cavity)')
    parser.add_argument('--test', action='store_true',
                       help='Quick test mode')
    
    args = parser.parse_args()
    
    if args.test:
        print("\n*** QUICK TEST MODE ***\n")
        args.molecules = 3
        args.equil = 0.5
        args.prod = 1.0
    
    run_simulation(
        lambda_coupling=args.lambda_coupling,
        cavity_freq_cm=args.cavity_freq,
        num_molecules=args.molecules,
        temperature_K=args.temp,
        dt_ps=args.dt,
        equilibration_ps=args.equil,
        production_ps=args.prod,
        output_prefix=args.output
    )


if __name__ == "__main__":
    main()
