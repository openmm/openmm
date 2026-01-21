#!/usr/bin/env python3
"""
Speed benchmark for Rabi splitting simulation.
No file I/O, just pure simulation speed measurement.
"""

import openmm
import openmm.unit as unit
import numpy as np
import time

# Unit conversions
HARTREE_TO_CM = 219474.63
BOHR_TO_NM = 0.0529177

# System parameters
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


def create_system():
    """Create the molecular system with cavity."""
    system = openmm.System()
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(BOX_SIZE_NM, 0, 0),
        openmm.Vec3(0, BOX_SIZE_NM, 0),
        openmm.Vec3(0, 0, BOX_SIZE_NM)
    )
    
    positions = []
    charges = []
    
    # Lattice placement
    n_per_side = int(np.ceil(NUM_MOLECULES ** (1/3)))
    spacing = BOX_SIZE_NM / n_per_side
    
    mol_count = 0
    for ix in range(n_per_side):
        for iy in range(n_per_side):
            for iz in range(n_per_side):
                if mol_count >= NUM_MOLECULES:
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
            if mol_count >= NUM_MOLECULES:
                break
        if mol_count >= NUM_MOLECULES:
            break
    
    num_molecular = system.getNumParticles()
    
    # Harmonic bonds
    bond_force = openmm.HarmonicBondForce()
    for i in range(NUM_MOLECULES):
        bond_force.addBond(2*i, 2*i+1, BOND_LENGTH_NM, K_OO)
    system.addForce(bond_force)
    
    # Nonbonded (needed for CavityForce to read charges)
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(1.0 * unit.nanometer)
    for i, q in enumerate(charges):
        nb.addParticle(q, 0.3, 0.0)
    for i in range(NUM_MOLECULES):
        nb.addException(2*i, 2*i+1, 0.0, 0.3, 0.0)
    system.addForce(nb)
    
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
    
    return system, positions, num_molecular, cavity_index


def benchmark(steps, batch_size=100, platform_name='CUDA'):
    """Run benchmark and return ns/day."""
    print("=" * 60)
    print(f"Speed Benchmark: {steps} steps, batch={batch_size}")
    print("=" * 60)
    
    np.random.seed(42)
    system, positions, num_molecular, cavity_index = create_system()
    
    print(f"System: {num_molecular} molecular atoms + 1 cavity particle")
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
    # Test different batch sizes
    print("\n" + "="*60)
    print("BATCH SIZE COMPARISON")
    print("="*60)
    
    test_steps = 10000
    
    for batch in [1, 10, 100, 1000]:
        speed = benchmark(test_steps, batch_size=batch)
        print()
    
    # Final long benchmark
    print("\n" + "="*60)
    print("FULL BENCHMARK (50000 steps)")
    print("="*60)
    benchmark(50000, batch_size=1000)
