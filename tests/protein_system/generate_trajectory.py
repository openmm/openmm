#!/usr/bin/env python3
"""
Generate PDB Trajectory from Cavity MD Simulation
==================================================

This script runs a cavity MD simulation and outputs PDB frames at specified
intervals (default: every 1 ps) for visualization in VMD, PyMOL, etc.

Usage:
    python generate_trajectory.py --pdb-path ./3UTL.pdb --lambda 0.01 --output traj
    python generate_trajectory.py --help
"""

import sys
import time
import argparse
from pathlib import Path

import numpy as np

try:
    from openmm import openmm
    from openmm import unit
    from openmm import app
    from pdbfixer import PDBFixer
    print("OpenMM and PDBFixer loaded successfully")
except ImportError as exc:
    print(f"Error importing required packages: {exc}")
    print("Install openmm, pdbfixer, and openmmforcefields.")
    sys.exit(1)

# Unit conversions
BOHR_TO_NM = 0.0529177
HARTREE_TO_KJMOL = 2625.5
HARTREE_TO_CM = 219474.63
AMU_TO_AU = 1822.888


def wavenumber_to_hartree(omega_cm):
    """Convert wavenumber (cm^-1) to Hartree."""
    return omega_cm / HARTREE_TO_CM


def prepare_system(pdb_path, padding_nm=1.0, ionic_strength=0.15, ph=7.0, verbose=False):
    """Prepare protein + solvent system."""
    print(f"\n--- Loading PDB: {pdb_path} ---")
    fixer = PDBFixer(filename=str(pdb_path))
    
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    print(f"  Prepared structure at pH {ph}")
    
    modeller = app.Modeller(fixer.topology, fixer.positions)
    forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip4pew.xml")
    
    print("\n--- Adding Solvent ---")
    modeller.addSolvent(
        forcefield,
        model="tip4pew",
        padding=padding_nm * unit.nanometer,
        ionicStrength=ionic_strength * unit.molar,
        neutralize=True,
    )
    print(f"  Padding: {padding_nm} nm, Ionic strength: {ionic_strength} M")
    
    print("\n--- Creating System ---")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=None,
        rigidWater=False,
        ewaldErrorTolerance=0.0005,
    )
    
    real_indices = np.array(
        [i for i in range(system.getNumParticles()) if not system.isVirtualSite(i)],
        dtype=int,
    )
    
    print(f"  Total particles: {system.getNumParticles()}")
    print(f"  Real atoms: {len(real_indices)}")
    
    return system, modeller.topology, modeller.getPositions(), real_indices


def add_cavity(system, positions, omega_cm, lambda_coupling):
    """Add cavity particle and coupling forces."""
    omega_au = wavenumber_to_hartree(omega_cm)
    photon_mass_amu = 1.0 / AMU_TO_AU
    
    cavity_index = system.addParticle(photon_mass_amu * unit.amu)
    positions.append(openmm.Vec3(0, 0, 0) * unit.nanometer)
    
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.addParticle(0.0, 0.1 * unit.nanometer, 0.0)
            break
    
    cavity_force = openmm.CavityForce(cavity_index, omega_au, lambda_coupling, photon_mass_amu)
    system.addForce(cavity_force)
    
    displacer = openmm.CavityParticleDisplacer(cavity_index, omega_au, photon_mass_amu)
    displacer.setSwitchOnStep(0)
    displacer.setSwitchOnLambda(lambda_coupling)
    system.addForce(displacer)
    
    print(f"\n--- Cavity Added ---")
    print(f"  Frequency: {omega_cm} cm^-1")
    print(f"  Lambda: {lambda_coupling}")
    
    return cavity_index


def run_trajectory_simulation(
    pdb_path,
    lambda_coupling=0.01,
    cavity_freq_cm=3663.0,
    temperature_K=300.0,
    dt_ps=0.0005,
    equilibration_ps=10.0,
    production_ps=100.0,
    frame_interval_ps=1.0,
    output_prefix="trajectory",
    padding_nm=1.0,
    ionic_strength=0.15,
    enable_cavity=True,
    verbose=False,
):
    """Run simulation and output PDB trajectory frames."""
    print("=" * 70)
    print("Cavity MD Simulation with PDB Trajectory Output")
    print("=" * 70)
    
    print(f"\nParameters:")
    print(f"  PDB: {pdb_path}")
    print(f"  Lambda: {lambda_coupling}")
    print(f"  Cavity freq: {cavity_freq_cm} cm^-1")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt_ps} ps ({dt_ps*1000:.1f} fs)")
    print(f"  Equilibration: {equilibration_ps} ps")
    print(f"  Production: {production_ps} ps")
    print(f"  Frame interval: {frame_interval_ps} ps")
    
    # Prepare system
    system, topology, positions, real_indices = prepare_system(
        pdb_path, padding_nm, ionic_strength, verbose=verbose
    )
    
    # Add cavity if enabled
    cavity_index = None
    if enable_cavity:
        cavity_index = add_cavity(system, positions, cavity_freq_cm, lambda_coupling)
    
    # Add thermostat
    try:
        bussi = openmm.BussiThermostat(temperature_K, 1.0)
        bussi.setApplyToAllParticles(False)
        for idx in real_indices:
            bussi.addParticle(int(idx))
        system.addForce(bussi)
        print(f"\n--- Bussi Thermostat: {temperature_K} K ---")
    except AttributeError:
        print("\n--- Using Langevin thermostat only ---")
    
    # Create integrator and context
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        0.5 / unit.picosecond,
        dt_ps * unit.picosecond,
    )
    
    platform = openmm.Platform.getPlatformByName("CUDA")
    properties = {"Precision": "mixed"}
    context = openmm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    
    print("\n--- Energy Minimization ---")
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=1000)
    state = context.getState(getEnergy=True)
    print(f"  Initial PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.1f} kJ/mol")
    
    # Calculate steps
    equilibration_steps = int(equilibration_ps / dt_ps)
    production_steps = int(production_ps / dt_ps)
    frame_interval_steps = int(frame_interval_ps / dt_ps)
    total_steps = equilibration_steps + production_steps
    n_frames = production_steps // frame_interval_steps
    
    print(f"\n--- Running Simulation ---")
    print(f"  Total steps: {total_steps:,}")
    print(f"  Equilibration: {equilibration_steps:,} steps")
    print(f"  Production: {production_steps:,} steps")
    print(f"  Frame interval: {frame_interval_steps:,} steps ({frame_interval_ps} ps)")
    print(f"  Expected frames: {n_frames}")
    
    # Create output directory
    output_dir = Path(output_prefix)
    output_dir.mkdir(exist_ok=True)
    
    # Number of atoms in topology (excludes cavity particle if present)
    n_topology_atoms = topology.getNumAtoms()
    
    # Save topology PDB
    topo_file = output_dir / "topology.pdb"
    state = context.getState(getPositions=True)
    pos = state.getPositions()
    # Exclude cavity particle (last particle) if present
    pos_for_pdb = pos[:n_topology_atoms]
    with open(topo_file, "w") as f:
        app.PDBFile.writeFile(topology, pos_for_pdb, f)
    print(f"\n  Saved topology: {topo_file}")
    
    start_time = time.time()
    frame_count = 0
    last_report_time = start_time
    last_report_step = 0
    
    for step in range(total_steps):
        integrator.step(1)
        
        # Progress report every 1000 steps
        if (step + 1) % 1000 == 0:
            phase = "EQUIL" if step < equilibration_steps else "PROD"
            sim_time_ps = (step + 1) * dt_ps
            pct = 100 * (step + 1) / total_steps
            
            current_time = time.time()
            dt_report = current_time - last_report_time
            steps_since_last = (step + 1) - last_report_step
            inst_rate = steps_since_last / max(dt_report, 1e-8)
            ns_per_day = (inst_rate * dt_ps * 86400) / 1000
            eta = (total_steps - step - 1) / max(inst_rate, 1e-8)
            
            last_report_time = current_time
            last_report_step = step + 1
            
            state = context.getState(getEnergy=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            
            print(
                f"[{pct:5.1f}%] {phase} | t={sim_time_ps:7.1f} ps | "
                f"PE: {pe:10.1f} kJ/mol | Frames: {frame_count} | "
                f"Speed: {ns_per_day:5.2f} ns/day | ETA: {eta/60:5.1f}m",
                flush=True,
            )
        
        # Save frame during production
        if step >= equilibration_steps:
            prod_step = step - equilibration_steps
            if prod_step % frame_interval_steps == 0:
                state = context.getState(getPositions=True)
                pos = state.getPositions()
                # Exclude cavity particle (last particle) if present
                pos_for_pdb = pos[:n_topology_atoms]
                
                frame_time_ps = prod_step * dt_ps
                frame_file = output_dir / f"frame_{frame_count:06d}_t{frame_time_ps:.1f}ps.pdb"
                
                with open(frame_file, "w") as f:
                    app.PDBFile.writeFile(topology, pos_for_pdb, f)
                
                frame_count += 1
    
    elapsed = time.time() - start_time
    print(f"\n--- Simulation Complete ---")
    print(f"  Wall time: {elapsed/60:.1f} min")
    print(f"  Frames saved: {frame_count}")
    print(f"  Output directory: {output_dir}")
    
    # Create a summary file
    summary_file = output_dir / "trajectory_info.txt"
    with open(summary_file, "w") as f:
        f.write(f"Trajectory Summary\n")
        f.write(f"==================\n\n")
        f.write(f"PDB source: {pdb_path}\n")
        f.write(f"Lambda coupling: {lambda_coupling}\n")
        f.write(f"Cavity frequency: {cavity_freq_cm} cm^-1\n")
        f.write(f"Temperature: {temperature_K} K\n")
        f.write(f"Timestep: {dt_ps} ps\n")
        f.write(f"Equilibration: {equilibration_ps} ps\n")
        f.write(f"Production: {production_ps} ps\n")
        f.write(f"Frame interval: {frame_interval_ps} ps\n")
        f.write(f"Total frames: {frame_count}\n")
        f.write(f"Total atoms: {system.getNumParticles()}\n")
        f.write(f"Real atoms: {len(real_indices)}\n")
    
    print(f"  Summary: {summary_file}")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="Generate PDB trajectory from cavity MD simulation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--pdb-path", type=str, required=True,
                        help="Path to input PDB file")
    parser.add_argument("--lambda", type=float, default=0.01, dest="lambda_coupling",
                        help="Cavity coupling strength")
    parser.add_argument("--cavity-freq", type=float, default=3663.0,
                        help="Cavity frequency in cm^-1")
    parser.add_argument("--temp", type=float, default=300.0,
                        help="Temperature in K")
    parser.add_argument("--dt", type=float, default=0.0005,
                        help="Timestep in ps")
    parser.add_argument("--equil", type=float, default=10.0,
                        help="Equilibration time in ps")
    parser.add_argument("--prod", type=float, default=100.0,
                        help="Production time in ps")
    parser.add_argument("--frame-interval", type=float, default=1.0,
                        help="Frame output interval in ps")
    parser.add_argument("--output", type=str, default="trajectory",
                        help="Output directory prefix")
    parser.add_argument("--padding", type=float, default=1.0,
                        help="Solvent padding in nm")
    parser.add_argument("--ionic-strength", type=float, default=0.15,
                        help="Ionic strength in M")
    parser.add_argument("--no-cavity", action="store_true",
                        help="Disable cavity coupling")
    parser.add_argument("--verbose", action="store_true",
                        help="Verbose output")
    
    args = parser.parse_args()
    
    run_trajectory_simulation(
        pdb_path=args.pdb_path,
        lambda_coupling=args.lambda_coupling,
        cavity_freq_cm=args.cavity_freq,
        temperature_K=args.temp,
        dt_ps=args.dt,
        equilibration_ps=args.equil,
        production_ps=args.prod,
        frame_interval_ps=args.frame_interval,
        output_prefix=args.output,
        padding_nm=args.padding,
        ionic_strength=args.ionic_strength,
        enable_cavity=not args.no_cavity,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
