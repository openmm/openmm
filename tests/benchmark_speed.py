#!/usr/bin/env python3
"""
Speed benchmark for Rabi splitting simulation.
No file I/O, just pure simulation speed measurement.
"""

import argparse
import openmm
import openmm.unit as unit
import numpy as np
import time

# Unit conversions
HARTREE_TO_CM = 219474.63
BOHR_TO_NM = 0.0529177

# System parameters (defaults)
NUM_MOLECULES = 250
BOX_SIZE_BOHR = 40.0
BOX_SIZE_NM = BOX_SIZE_BOHR * BOHR_TO_NM

# Cavity parameters
CAVITY_FREQ_CM = 1560
OMEGAC_AU = CAVITY_FREQ_CM / HARTREE_TO_CM
LAMBDA_COUPLING = 0.07
AMU_TO_AU = 1822.888
PHOTON_MASS = 1.0 / AMU_TO_AU

# Molecular parameters
MASS_O = 16.0
CHARGE_MAGNITUDE = 0.3
TEMPERATURE_K = 100.0

# Force constant for O-O (1560 cm⁻¹)
K_OO_AU = 0.73204
K_OO = K_OO_AU * 2625.5 / (BOHR_TO_NM ** 2)
BOND_LENGTH_NM = 0.121

# Timing
DT_PS = 0.001
FRICTION = 0.01


def compute_box_bohr(num_molecules, base_num=NUM_MOLECULES, base_box_bohr=BOX_SIZE_BOHR):
    """Scale box length to keep density constant."""
    scale = (num_molecules / base_num) ** (1.0 / 3.0)
    return base_box_bohr * scale


def create_system(num_molecules=NUM_MOLECULES, box_bohr=BOX_SIZE_BOHR, include_cavity=True):
    """Create the molecular system with cavity."""
    box_size_nm = box_bohr * BOHR_TO_NM
    system = openmm.System()
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_size_nm, 0, 0),
        openmm.Vec3(0, box_size_nm, 0),
        openmm.Vec3(0, 0, box_size_nm)
    )
    
    positions = []
    charges = []
    
    # Lattice placement
    n_per_side = int(np.ceil(num_molecules ** (1/3)))
    spacing = box_size_nm / n_per_side
    
    mol_count = 0
    for ix in range(n_per_side):
        for iy in range(n_per_side):
            for iz in range(n_per_side):
                if mol_count >= num_molecules:
                    break
                
                cx = (ix + 0.5) * spacing
                cy = (iy + 0.5) * spacing
                cz = (iz + 0.5) * spacing
                
                # Random orientation
                theta = np.random.uniform(0, np.pi)
                phi = np.random.uniform(0, 2 * np.pi)
                dx = BOND_LENGTH_NM / 2 * np.sin(theta) * np.cos(phi)
                dy = BOND_LENGTH_NM / 2 * np.sin(theta) * np.sin(phi)
                dz = BOND_LENGTH_NM / 2 * np.cos(theta)
                
                system.addParticle(MASS_O)
                system.addParticle(MASS_O)
                positions.append(openmm.Vec3(cx - dx, cy - dy, cz - dz))
                positions.append(openmm.Vec3(cx + dx, cy + dy, cz + dz))
                charges.extend([-CHARGE_MAGNITUDE, CHARGE_MAGNITUDE])
                
                mol_count += 1
            if mol_count >= num_molecules:
                break
        if mol_count >= num_molecules:
            break
    
    num_molecular = system.getNumParticles()
    
    # Harmonic bonds
    bond_force = openmm.HarmonicBondForce()
    for i in range(num_molecules):
        bond_force.addBond(2*i, 2*i+1, BOND_LENGTH_NM, K_OO)
    system.addForce(bond_force)
    
    # Nonbonded (needed for CavityForce to read charges)
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(1.0 * unit.nanometer)
    for i, q in enumerate(charges):
        nb.addParticle(q, 0.3, 0.0)
    for i in range(num_molecules):
        nb.addException(2*i, 2*i+1, 0.0, 0.3, 0.0)
    system.addForce(nb)
    
    cavity_index = None
    if include_cavity:
        # Cavity particle
        cavity_index = system.addParticle(PHOTON_MASS)
        positions.append(openmm.Vec3(0, 0, 0))
        nb.addParticle(0.0, 0.1, 0.0)  # Cavity has no charge/LJ
        
        # CavityForce
        cavity_force = openmm.CavityForce(cavity_index, OMEGAC_AU, LAMBDA_COUPLING, PHOTON_MASS)
        system.addForce(cavity_force)
        
        # CavityParticleDisplacer
        displacer = openmm.CavityParticleDisplacer(cavity_index, OMEGAC_AU, PHOTON_MASS)
        displacer.setSwitchOnLambda(LAMBDA_COUPLING)
        system.addForce(displacer)
    
    return system, positions, num_molecular, cavity_index, box_size_nm


def benchmark(steps, batch_size=100, platform_name='CUDA', num_molecules=NUM_MOLECULES, box_bohr=BOX_SIZE_BOHR, include_cavity=True):
    """Run benchmark and return ns/day."""
    print("=" * 60)
    print(f"Speed Benchmark: {steps} steps, batch={batch_size}")
    print("=" * 60)
    
    np.random.seed(42)
    system, positions, num_molecular, cavity_index, box_size_nm = create_system(
        num_molecules=num_molecules,
        box_bohr=box_bohr,
        include_cavity=include_cavity
    )
    
    cavity_label = " + 1 cavity particle" if include_cavity else ""
    print(f"System: {num_molecular} molecular atoms{cavity_label}")
    print(f"Box: {box_bohr:.1f} Bohr = {box_size_nm:.3f} nm")
    print(f"Platform: {platform_name}")
    
    integrator = openmm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        FRICTION / unit.picosecond,
        DT_PS * unit.picoseconds
    )
    
    platform = openmm.Platform.getPlatformByName(platform_name)
    properties = {'Precision': 'mixed'} if platform_name == 'CUDA' else {}
    
    context = openmm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin)
    
    # Warmup
    print("\nWarmup (1000 steps)...")
    integrator.step(1000)
    
    # Benchmark
    print(f"\nBenchmarking {steps} steps...")
    
    start = time.perf_counter()
    
    remaining = steps
    while remaining > 0:
        chunk = min(batch_size, remaining)
        integrator.step(chunk)
        remaining -= chunk
    
    elapsed = time.perf_counter() - start
    
    # Calculate performance
    time_simulated_ps = steps * DT_PS
    time_simulated_ns = time_simulated_ps / 1000
    ns_per_day = time_simulated_ns / elapsed * 86400
    
    print(f"\n{'='*60}")
    print(f"RESULTS:")
    print(f"  Steps: {steps}")
    print(f"  Simulated time: {time_simulated_ps:.1f} ps = {time_simulated_ns:.3f} ns")
    print(f"  Wall time: {elapsed:.2f} s")
    print(f"  Speed: {ns_per_day:.1f} ns/day")
    print(f"  Steps/sec: {steps/elapsed:.0f}")
    print(f"{'='*60}")
    
    return ns_per_day


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Speed benchmark for cavity MD.")
    parser.add_argument("--molecules", type=int, default=NUM_MOLECULES, help="Number of dimers.")
    parser.add_argument("--box-bohr", type=float, default=BOX_SIZE_BOHR, help="Box length in Bohr.")
    parser.add_argument("--scale-density", action="store_true",
                        help="Scale box to keep density constant vs default.")
    parser.add_argument("--steps", type=int, default=10000, help="Total steps for benchmark.")
    parser.add_argument("--batch", type=int, default=1000, help="Batch size for integrator.step().")
    parser.add_argument("--no-cavity", action="store_true",
                        help="Disable cavity particle and cavity forces.")
    parser.add_argument("--full-suite", action="store_true",
                        help="Run default batch-size comparison and long benchmark.")
    args = parser.parse_args()

    box_bohr = args.box_bohr
    if args.scale_density:
        box_bohr = compute_box_bohr(args.molecules)

    if args.full_suite:
        print("\n" + "="*60)
        print("BATCH SIZE COMPARISON")
        print("="*60)

        for batch in [1, 10, 100, 1000]:
            benchmark(args.steps, batch_size=batch, num_molecules=args.molecules, box_bohr=box_bohr, include_cavity=not args.no_cavity)
            print()

        print("\n" + "="*60)
        print("FULL BENCHMARK (50000 steps)")
        print("="*60)
        benchmark(50000, batch_size=1000, num_molecules=args.molecules, box_bohr=box_bohr, include_cavity=not args.no_cavity)
    else:
        benchmark(args.steps, batch_size=args.batch, num_molecules=args.molecules, box_bohr=box_bohr, include_cavity=not args.no_cavity)
