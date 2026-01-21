#!/usr/bin/env python3
"""
Compare speed with and without CavityForce to identify bottleneck.
"""

import openmm
import openmm.unit as unit
import numpy as np
import time

# Parameters
NUM_MOLECULES = 250
BOX_SIZE_NM = 2.117
MASS_O = 16.0
CHARGE_MAGNITUDE = 0.3
TEMPERATURE_K = 100.0
K_OO = 686349.0
BOND_LENGTH_NM = 0.121
DT_PS = 0.001
FRICTION = 0.01

# Cavity
OMEGAC_AU = 1560 / 219474.63
LAMBDA_COUPLING = 0.07
PHOTON_MASS = 1.0 / 1822.888


def create_system(with_cavity=True):
    system = openmm.System()
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(BOX_SIZE_NM, 0, 0),
        openmm.Vec3(0, BOX_SIZE_NM, 0),
        openmm.Vec3(0, 0, BOX_SIZE_NM)
    )
    
    positions = []
    charges = []
    
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
    
    # Harmonic bonds
    bond_force = openmm.HarmonicBondForce()
    for i in range(NUM_MOLECULES):
        bond_force.addBond(2*i, 2*i+1, BOND_LENGTH_NM, K_OO)
    system.addForce(bond_force)
    
    # Nonbonded
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(1.0 * unit.nanometer)
    for q in charges:
        nb.addParticle(q, 0.3, 0.0)
    for i in range(NUM_MOLECULES):
        nb.addException(2*i, 2*i+1, 0.0, 0.3, 0.0)
    
    if with_cavity:
        # Add cavity particle
        cavity_index = system.addParticle(PHOTON_MASS)
        positions.append(openmm.Vec3(0, 0, 0))
        nb.addParticle(0.0, 0.1, 0.0)
        
        # CavityForce
        cavity_force = openmm.CavityForce(cavity_index, OMEGAC_AU, LAMBDA_COUPLING, PHOTON_MASS)
        system.addForce(cavity_force)
        
        # Displacer
        displacer = openmm.CavityParticleDisplacer(cavity_index, OMEGAC_AU, PHOTON_MASS)
        displacer.setSwitchOnLambda(LAMBDA_COUPLING)
        system.addForce(displacer)
    
    system.addForce(nb)
    
    return system, positions


def benchmark(with_cavity, steps=50000, platform_name='CUDA'):
    label = "WITH CavityForce" if with_cavity else "WITHOUT CavityForce"
    print(f"\n{'='*60}")
    print(f"Benchmark: {label}")
    print(f"{'='*60}")
    
    np.random.seed(42)
    system, positions = create_system(with_cavity=with_cavity)
    
    print(f"Particles: {system.getNumParticles()}")
    print(f"Forces: {[system.getForce(i).__class__.__name__ for i in range(system.getNumForces())]}")
    
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
    integrator.step(1000)
    
    # Benchmark
    start = time.perf_counter()
    integrator.step(steps)
    elapsed = time.perf_counter() - start
    
    ns_per_day = (steps * DT_PS / 1000) / elapsed * 86400
    
    print(f"\nResults:")
    print(f"  Steps: {steps}")
    print(f"  Wall time: {elapsed:.2f} s")
    print(f"  Speed: {ns_per_day:.1f} ns/day")
    print(f"  Steps/sec: {steps/elapsed:.0f}")
    
    return ns_per_day


if __name__ == "__main__":
    print("="*60)
    print("COMPARING SIMULATION SPEED WITH/WITHOUT CAVITYFORCE")
    print("="*60)
    
    speed_no_cavity = benchmark(with_cavity=False)
    speed_with_cavity = benchmark(with_cavity=True)
    
    print(f"\n{'='*60}")
    print("SUMMARY:")
    print(f"  Without CavityForce: {speed_no_cavity:.1f} ns/day")
    print(f"  With CavityForce:    {speed_with_cavity:.1f} ns/day")
    print(f"  Slowdown factor:     {speed_no_cavity/speed_with_cavity:.2f}x")
    print(f"{'='*60}")
