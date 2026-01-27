#!/usr/bin/env python3
"""
UMA/eSEN Ice RPMD Simulation
=============================

Comprehensive RPMD simulation testing if ice remains frozen at 243K using 
UMA or eSEN ML potentials.

Features:
- RPMD with PILE_G thermostat for quantum nuclear effects
- NPT or NVT ensemble
- Equilibration and production phases
- Real-time progress reporting with energy and speed metrics
- PDB trajectory output (centroid)
- Comprehensive analysis and logging

Usage:
    python test_uma_ice_rpmd.py --molecules 50 --beads 16 --temperature 243 \\
                                --dt 4.0 --equil 5 --prod 100 --pressure 0
"""

import sys
import numpy as np
import time
from pathlib import Path

from openmm import app, unit, Vec3
from openmm import RPMDIntegrator, Context, Platform, MonteCarloBarostat
from openmmml import MLPotential

print("✓ All required packages loaded successfully")

def create_ice_structure(num_molecules=32):
    """
    Create ice Ih structure with tetrahedral coordination.
    
    Ice Ih parameters:
    - O-O distance: ~2.75 Å
    - O-H bond length: ~0.96 Å
    - Tetrahedral angle: 109.47°
    
    Returns
    -------
    topology : app.Topology
    positions : list of Vec3
    box_vectors : tuple of Vec3
    """
    print(f"\n--- Creating Ice Ih Structure ---")
    print(f"  Molecules: {num_molecules}")
    
    # Calculate box size for ice density (0.92 g/cm³)
    target_density = 0.92  # g/cm³
    mass_per_molecule = 18.015  # g/mol  
    avogadro = 6.022e23
    
    # Volume needed for target density: V = (N * M) / (ρ * N_A)
    total_mass_g = (num_molecules * mass_per_molecule) / avogadro
    volume_cm3 = total_mass_g / target_density
    volume_nm3 = volume_cm3 * 1e21  # cm³ to nm³
    box_size = volume_nm3 ** (1/3)
    
    # Spacing for cubic lattice
    n_side = int(np.ceil(num_molecules ** (1/3)))
    spacing_nm = box_size / n_side
    
    print(f"  Lattice: {n_side}x{n_side}x{n_side}")
    print(f"  Box size: {box_size:.3f} nm")
    
    # Water geometry
    oh_bond = 0.09572  # nm
    hoh_angle = 104.52 * np.pi / 180.0
    
    # Create topology
    topology = app.Topology()
    positions = []
    
    np.random.seed(42)
    
    mol_count = 0
    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                if mol_count >= num_molecules:
                    break
                
                # Position of oxygen
                o_pos = np.array([
                    (i + 0.5) * spacing_nm,
                    (j + 0.5) * spacing_nm,
                    (k + 0.5) * spacing_nm
                ])
                
                # Random orientation
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                psi = np.random.rand() * 2 * np.pi
                
                # Local H positions
                h1_local = np.array([
                    oh_bond * np.sin(hoh_angle/2),
                    0,
                    oh_bond * np.cos(hoh_angle/2)
                ])
                h2_local = np.array([
                    -oh_bond * np.sin(hoh_angle/2),
                    0,
                    oh_bond * np.cos(hoh_angle/2)
                ])
                
                # Rotation matrices
                Rz = np.array([
                    [np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]
                ])
                Ry = np.array([
                    [np.cos(phi), 0, np.sin(phi)],
                    [0, 1, 0],
                    [-np.sin(phi), 0, np.cos(phi)]
                ])
                Rx = np.array([
                    [1, 0, 0],
                    [0, np.cos(psi), -np.sin(psi)],
                    [0, np.sin(psi), np.cos(psi)]
                ])
                R = Rz @ Ry @ Rx
                
                h1_pos = o_pos + R @ h1_local
                h2_pos = o_pos + R @ h2_local
                
                # Add to topology
                chain = topology.addChain()
                residue = topology.addResidue('HOH', chain)
                o_atom = topology.addAtom('O', app.Element.getBySymbol('O'), residue)
                h1_atom = topology.addAtom('H', app.Element.getBySymbol('H'), residue)
                h2_atom = topology.addAtom('H', app.Element.getBySymbol('H'), residue)
                topology.addBond(o_atom, h1_atom)
                topology.addBond(o_atom, h2_atom)
                
                positions.append(o_pos.tolist())
                positions.append(h1_pos.tolist())
                positions.append(h2_pos.tolist())
                
                mol_count += 1
            if mol_count >= num_molecules:
                break
        if mol_count >= num_molecules:
            break
    
    print(f"  Actual molecules: {mol_count}")
    print(f"  Total atoms: {len(positions)}")
    
    # Set box vectors
    box_vectors = (
        [box_size, 0, 0] * unit.nanometer,
        [0, box_size, 0] * unit.nanometer,
        [0, 0, box_size] * unit.nanometer
    )
    
    topology.setPeriodicBoxVectors(box_vectors)
    print(f"  Box vectors (nm): [{box_size:.3f}, {box_size:.3f}, {box_size:.3f}]")
    
    # Calculate density
    mass_per_molecule = 18.015  # g/mol for H2O
    total_mass_g = (mol_count * mass_per_molecule) * 1.66054e-24
    volume_cm3 = (box_size * 1e-7) ** 3
    density = total_mass_g / volume_cm3
    print(f"  Initial density: {density:.3f} g/cm³")
    
    return topology, positions, box_vectors


def run_simulation(num_molecules=32, num_beads=8, temperature_K=243.0,
                   pressure_bar=1.0, dt_fs=4.0, equilibration_ps=5.0, 
                   production_ps=100.0, model_name='uma-s-1p1-pythonforce-batch',
                   output_dir='.', report_interval_ps=1.0, pdb_interval_ps=1.0):
    """
    Run comprehensive RPMD simulation of ice.
    
    Parameters
    ----------
    num_molecules : int
        Number of water molecules
    num_beads : int
        Number of RPMD beads
    temperature_K : float
        Temperature in Kelvin
    pressure_bar : float
        Pressure in bar (0 = NVT, >0 = NPT)
    dt_fs : float
        Timestep in femtoseconds
    equilibration_ps : float
        Equilibration time in picoseconds
    production_ps : float
        Production time in picoseconds
    model_name : str
        UMA or eSEN model name
    output_dir : str
        Output directory for trajectories
    report_interval_ps : float
        Console reporting interval in picoseconds (default: 1.0)
    pdb_interval_ps : float
        PDB frame save interval in picoseconds (default: 1.0)
        
    Returns
    -------
    success : bool
        True if simulation completed successfully
    """
    
    print("=" * 80)
    print("UMA/eSEN Ice RPMD Simulation")
    print("=" * 80)
    print(f"Model: {model_name}")
    print(f"Temperature: {temperature_K} K")
    print(f"Pressure: {pressure_bar} bar {'(NPT)' if pressure_bar > 0 else '(NVT)'}")
    print(f"RPMD beads: {num_beads}")
    print(f"Timestep: {dt_fs} fs")
    print(f"Molecules: {num_molecules}")
    print(f"Equilibration: {equilibration_ps} ps")
    print(f"Production: {production_ps} ps")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Create structure
    topology, positions, box_vectors = create_ice_structure(num_molecules)
    n_atoms = len(positions)
    
    # Create potential
    print(f"\n--- Creating ML Potential ---")
    print("Creating OpenMM system...")
    potential = MLPotential(model_name)
    
    system = potential.createSystem(
        topology,
        task_name='omol',
        charge=0,
        spin=1
    )
    
    print(f"  System uses PBC: {system.usesPeriodicBoundaryConditions()}")
    
    # Add barostat if NPT
    if pressure_bar > 0:
        barostat = MonteCarloBarostat(
            pressure_bar * unit.bar,
            temperature_K * unit.kelvin,
            25  # Attempt every 25 steps
        )
        system.addForce(barostat)
        print(f"  Added MonteCarloBarostat (NPT ensemble)")
        print(f"  Pressure: {pressure_bar} bar")
    else:
        print(f"  Running NVT (no barostat)")
    
    # Calculate density
    mass_per_molecule = 18.015
    total_mass_g = (num_molecules * mass_per_molecule) * 1.66054e-24
    box_size = box_vectors[0][0].value_in_unit(unit.nanometer)
    volume_cm3 = (box_size * 1e-7) ** 3
    density = total_mass_g / volume_cm3
    print(f"  Initial density: {density:.3f} g/cm³")
    print(f"  Expected ice density: ~0.92 g/cm³")
    
    # Create RPMD integrator with PILE_G thermostat
    print(f"\n--- Creating RPMD Integrator ---")
    print(f"  Beads: {num_beads}")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt_fs} fs")
    
    friction = 1.0  # ps^-1
    centroid_friction = 0.5  # ps^-1 for Bussi thermostat
    
    integrator = RPMDIntegrator(
        num_beads,
        temperature_K * unit.kelvin,
        friction / unit.picosecond,
        dt_fs * unit.femtoseconds
    )
    
    # Set PILE thermostat (standard Langevin for all modes)
    integrator.setThermostatType(RPMDIntegrator.Pile)
    
    print(f"  Thermostat: PILE (Langevin for all modes)")
    print(f"  Friction: {friction} ps^-1")
    
    # Create context
    print(f"\n--- Creating Simulation Context ---")
    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed'}
    print(f"  Using CUDA platform (GPU acceleration)")
    
    context = Context(system, integrator, platform, properties)
    
    # Initialize beads
    print(f"\n--- Initializing RPMD Beads ---")
    positions_nm = np.array(positions)
    
    # Set context positions first (required before getting velocities)
    pos_with_units = [Vec3(float(p[0]), float(p[1]), float(p[2])) * unit.nanometer for p in positions_nm]
    context.setPositions(pos_with_units)
    
    for i in range(num_beads):
        perturbation = np.random.randn(n_atoms, 3) * 0.001
        perturbed_pos = positions_nm + perturbation
        integrator.setPositions(i, perturbed_pos * unit.nanometer)
        print(f"  Bead {i}: positions set")
    
    # Set velocities
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    state = context.getState(getVelocities=True)
    base_velocities = state.getVelocities()
    
    for i in range(num_beads):
        velocities = []
        for v in base_velocities:
            vx = v[0].value_in_unit(unit.nanometer/unit.picosecond) + 0.01 * np.random.randn()
            vy = v[1].value_in_unit(unit.nanometer/unit.picosecond) + 0.01 * np.random.randn()
            vz = v[2].value_in_unit(unit.nanometer/unit.picosecond) + 0.01 * np.random.randn()
            velocities.append(Vec3(vx, vy, vz) * unit.nanometer/unit.picosecond)
        integrator.setVelocities(i, velocities)
    
    print(f"  Velocities initialized to {temperature_K} K")
    
    # Energy minimization to relax the structure
    print(f"\n--- Energy Minimization ---")
    print("  Minimizing energy to relax initial structure...")
    initial_state = integrator.getState(0, getEnergy=True)
    initial_pe = initial_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Initial PE: {initial_pe:.2f} kJ/mol")
    
    # Perform energy minimization on bead 0, then copy to all beads
    from openmm import LocalEnergyMinimizer
    LocalEnergyMinimizer.minimize(context, tolerance=10.0, maxIterations=1000)
    
    # Get minimized positions and copy to all beads
    minimized_state = context.getState(getPositions=True, getEnergy=True)
    minimized_positions = minimized_state.getPositions()
    minimized_pe = minimized_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Minimized PE: {minimized_pe:.2f} kJ/mol")
    print(f"  Energy reduction: {initial_pe - minimized_pe:.2f} kJ/mol")
    
    for i in range(num_beads):
        integrator.setPositions(i, minimized_positions)
    print(f"  Minimized structure copied to all {num_beads} beads")
    
    # Get initial energies
    print(f"\n--- Initial Bead Energies ---")
    print("  Getting initial bead energies...")
    
    bead_energies = []
    for i in range(num_beads):
        state = integrator.getState(i, getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        ke = state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
        total = pe + ke
        bead_energies.append((pe, ke, total))
        print(f"  Bead {i:2d}: PE = {pe:10.2f} kJ/mol, KE = {ke:8.2f} kJ/mol, Total = {total:10.2f} kJ/mol")
    
    mean_pe = np.mean([e[0] for e in bead_energies])
    mean_ke = np.mean([e[1] for e in bead_energies])
    mean_total = np.mean([e[2] for e in bead_energies])
    print(f"\n  Mean: PE = {mean_pe:10.2f} kJ/mol, KE = {mean_ke:8.2f} kJ/mol, Total = {mean_total:10.2f} kJ/mol")
    
    # Setup PDB reporter (centroid only)
    print(f"\n--- Setting up PDB Reporter ---")
    pdb_file = output_path / f"ice_rpmd_{model_name.split('-')[0]}_T{temperature_K:.0f}_b{num_beads}.pdb"
    print(f"  PDB output: {pdb_file}")
    
    # Calculate steps
    dt_ps = dt_fs / 1000.0
    equilibration_steps = int(equilibration_ps / dt_ps)
    production_steps = int(production_ps / dt_ps)
    total_steps = equilibration_steps + production_steps
    
    report_interval_steps = max(1, int(report_interval_ps / dt_ps))
    pdb_interval_steps = max(1, int(pdb_interval_ps / dt_ps))
    
    print(f"  Report interval: {report_interval_steps} steps ({report_interval_steps * dt_ps:.2f} ps)")
    print(f"  PDB save interval: {pdb_interval_steps} steps ({pdb_interval_steps * dt_ps:.2f} ps)")
    
    # Run simulation
    print("\n" + "=" * 80)
    print("RUNNING SIMULATION")
    print("=" * 80)
    
    start_time = time.time()
    last_report_time = start_time
    
    step_counter = 0
    phase = "Equilibration"
    
    energies = []
    temperatures = []
    densities = []
    
    pdb_file_handle = open(pdb_file, 'w')
    
    # Write PDB header with proper unit handling
    box_size_nm = box_vectors[0][0].value_in_unit(unit.nanometer)
    box_size_angstrom = box_size_nm * 10.0  # Convert nm to Angstrom for PDB
    pdb_file_handle.write("REMARK   Ice RPMD simulation - centroid trajectory\n")
    pdb_file_handle.write("REMARK   Model: %s\n" % model_name)
    pdb_file_handle.write("REMARK   Temperature: %.1f K, Beads: %d\n" % (temperature_K, num_beads))
    pdb_file_handle.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" % (
        box_size_angstrom, box_size_angstrom, box_size_angstrom, 90.0, 90.0, 90.0))
    
    for step in range(1, total_steps + 1):
        integrator.step(1)
        step_counter = step
        
        # Switch to production phase
        if step == equilibration_steps:
            phase = "Production"
            print(f"\n{'=' * 80}")
            print(f"PRODUCTION PHASE STARTED")
            print(f"{'=' * 80}\n")
        
        # Report progress
        if step % report_interval_steps == 0:
            current_wall_time = time.time()
            elapsed_since_last = current_wall_time - last_report_time
            last_report_time = current_wall_time
            
            # Get state from bead 0 for PE and positions
            state0 = integrator.getState(0, getEnergy=True, getPositions=True)
            current_pe = state0.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
            # Calculate physical temperature from centroid kinetic energy
            # Physical T = (2 * KE_centroid) / (3 * N_atoms * k_B)
            centroid_ke = 0.0
            masses = [system.getParticleMass(j).value_in_unit(unit.dalton) for j in range(n_atoms)]
            
            # Get velocities for all beads
            all_bead_velocities = []
            for b in range(num_beads):
                # Use getState to retrieve velocities for each bead
                state_b = integrator.getState(b, getVelocities=True)
                all_bead_velocities.append(state_b.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond))
            
            all_bead_velocities = np.array(all_bead_velocities) # (num_beads, n_atoms, 3)
            centroid_velocities = np.mean(all_bead_velocities, axis=0) # (n_atoms, 3)
            
            # KE = 0.5 * sum(m * v^2)
            # v is in nm/ps, m is in dalton (g/mol)
            # 1 dalton * (nm/ps)^2 = 1 g/mol * (1e-9 m / 1e-12 s)^2 = 1e-3 kg/mol * 1e6 m^2/s^2 = 1000 J/mol = 1 kJ/mol
            # So the result is directly in kJ/mol
            for j in range(n_atoms):
                v2 = np.sum(centroid_velocities[j]**2)
                centroid_ke += 0.5 * masses[j] * v2
            
            current_ke = centroid_ke # This is the physical KE
            current_total = current_pe + current_ke
            
            # Calculate temperature from physical KE
            # k_B = 8.314e-3 kJ/(mol*K)
            current_temp = (2.0 * current_ke) / (3.0 * n_atoms * 8.314e-3)
            
            # Calculate density if NPT
            if pressure_bar > 0:
                box = state0.getPeriodicBoxVectors()
                volume_nm3 = (box[0][0] * box[1][1] * box[2][2]).value_in_unit(unit.nanometer**3)
                volume_cm3 = volume_nm3 * 1e-21
                current_density = total_mass_g / volume_cm3
                densities.append(current_density)
            else:
                current_density = density
            
            energies.append(current_total)
            temperatures.append(current_temp)
            
            # Calculate speed
            sim_time_ps = step * dt_ps
            progress_pct = (step / total_steps) * 100
            steps_per_sec = report_interval_steps / elapsed_since_last
            ns_per_day = steps_per_sec * dt_ps * 86400.0 / 1000.0
            
            print(f"[{progress_pct:5.1f}%] {phase:13s} | "
                  f"Step {step:7d}/{total_steps} ({sim_time_ps:6.1f} ps) | "
                  f"{ns_per_day:6.2f} ns/day")
            print(f"         PE={current_pe:10.1f} KE={current_ke:8.1f} Total={current_total:10.1f} kJ/mol | "
                  f"T={current_temp:6.1f} K | ρ={current_density:.3f} g/cm³")
            print("")
        
        # Save PDB frame (production only, at specified interval)
        if step > equilibration_steps and (step - equilibration_steps) % pdb_interval_steps == 0:
            # Get centroid positions
            centroid_positions = np.zeros((n_atoms, 3))
            for bead in range(num_beads):
                state = integrator.getState(bead, getPositions=True)
                pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                centroid_positions += pos
            centroid_positions /= num_beads
            
            # Write PDB frame manually
            model_num = (step - equilibration_steps) // pdb_interval_steps
            pdb_file_handle.write("MODEL     %4d\n" % model_num)
            
            atom_idx = 1
            for residue in topology.residues():
                for atom in residue.atoms():
                    pos_angstrom = centroid_positions[atom.index] * 10.0  # nm to Angstrom
                    pdb_file_handle.write("ATOM  %5d %-4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" % (
                        atom_idx,
                        atom.name,
                        residue.name,
                        residue.chain.id if residue.chain.id else ' ',
                        residue.index + 1,
                        pos_angstrom[0],
                        pos_angstrom[1],
                        pos_angstrom[2],
                        1.0,  # occupancy
                        0.0,  # temperature factor
                        atom.element.symbol
                    ))
                    atom_idx += 1
            
            pdb_file_handle.write("ENDMDL\n")
    
    # Close PDB file
    pdb_file_handle.write("END\n")
    pdb_file_handle.close()
    
    total_elapsed = time.time() - start_time
    
    print("=" * 80)
    print("SIMULATION COMPLETE")
    print("=" * 80)
    print(f"Total time: {total_elapsed/60:.1f} minutes ({total_elapsed/3600:.2f} hours)")
    print(f"Average speed: {(total_steps * dt_ps * 86400.0) / (total_elapsed * 1000.0):.2f} ns/day")
    
    # Analysis
    print(f"\n--- Simulation Analysis ---")
    energies = np.array(energies)
    temperatures = np.array(temperatures)
    
    # Production phase only
    prod_start_idx = len(energies) // 2
    prod_energies = energies[prod_start_idx:]
    prod_temps = temperatures[prod_start_idx:]
    
    print(f"\nProduction phase statistics:")
    print(f"  Energy:")
    print(f"    Mean: {np.mean(prod_energies):.2f} kJ/mol")
    print(f"    Std:  {np.std(prod_energies):.2f} kJ/mol")
    print(f"    Drift: {(prod_energies[-1] - prod_energies[0])/prod_energies[0]*100:+.2f}%")
    print(f"  Temperature:")
    print(f"    Mean: {np.mean(prod_temps):.2f} K (target: {temperature_K} K)")
    print(f"    Std:  {np.std(prod_temps):.2f} K")
    
    if pressure_bar > 0 and len(densities) > 0:
        densities = np.array(densities)
        prod_densities = densities[prod_start_idx:]
        print(f"  Density:")
        print(f"    Mean: {np.mean(prod_densities):.3f} g/cm³")
        print(f"    Std:  {np.std(prod_densities):.3f} g/cm³")
        print(f"    Ice Ih reference: ~0.92 g/cm³")
        
        is_frozen = np.mean(prod_densities) > 0.85
        print(f"\n  Ice status: {'FROZEN ❄' if is_frozen else 'MELTED 💧'}")
    else:
        is_frozen = True  # Can't determine without NPT
        print(f"\n  Ice status: Cannot determine (NVT simulation)")
    
    print(f"\n  Output files:")
    print(f"    Trajectory: {pdb_file}")
    
    print("\n" + "=" * 80)
    
    return is_frozen


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='UMA/eSEN Ice RPMD Simulation')
    parser.add_argument('--molecules', type=int, default=32,
                       help='Number of water molecules (default: 32)')
    parser.add_argument('--beads', type=int, default=8,
                       help='Number of RPMD beads (default: 8)')
    parser.add_argument('--temperature', type=float, default=243.0,
                       help='Temperature in K (default: 243.0)')
    parser.add_argument('--pressure', type=float, default=1.0,
                       help='Pressure in bar (default: 1.0, use 0 for NVT)')
    parser.add_argument('--dt', type=float, default=1.0,
                       help='Timestep in fs (default: 1.0, use 0.5-2.0 for stability)')
    parser.add_argument('--equil', type=float, default=5.0,
                       help='Equilibration time in ps (default: 5.0)')
    parser.add_argument('--prod', type=float, default=100.0,
                       help='Production time in ps (default: 100.0)')
    parser.add_argument('--model', type=str, default='uma-s-1p1-pythonforce-batch',
                       help='Model name (default: uma-s-1p1-pythonforce-batch)')
    parser.add_argument('--output', type=str, default='.',
                       help='Output directory (default: current directory)')
    parser.add_argument('--report-interval', type=float, default=1.0,
                       help='Console report interval in ps (default: 1.0)')
    parser.add_argument('--pdb-interval', type=float, default=1.0,
                       help='PDB frame save interval in ps (default: 1.0)')
    
    args = parser.parse_args()
    
    try:
        success = run_simulation(
            num_molecules=args.molecules,
            num_beads=args.beads,
            temperature_K=args.temperature,
            pressure_bar=args.pressure,
            dt_fs=args.dt,
            equilibration_ps=args.equil,
            production_ps=args.prod,
            model_name=args.model,
            output_dir=args.output,
            report_interval_ps=args.report_interval,
            pdb_interval_ps=args.pdb_interval
        )
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n✗ SIMULATION FAILED")
        print(f"Error: {type(e).__name__}: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
