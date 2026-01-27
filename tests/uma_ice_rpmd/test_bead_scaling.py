#!/usr/bin/env python3
"""Test RPMD bead scaling to verify parallelization efficiency."""

import sys
import time
import numpy as np
sys.path.insert(0, '/media/extradrive/Trajectories/openmm/fairchem/src')
sys.path.insert(0, '/media/extradrive/Trajectories/openmm/fairchem/openmm-ml')

import openmm
from openmm import unit
from openmm.app import PDBFile, Topology, Element
from openmmml.models import MLPotential
import torch

print("=" * 80)
print("UMA RPMD Bead Scaling Test")
print("=" * 80)

# Create simple water system (10 molecules = 30 atoms)
topology = Topology()
chain = topology.addChain()
for i in range(10):
    residue = topology.addResidue(f'HOH', chain)
    O = topology.addAtom('O', Element.getBySymbol('O'), residue)
    H1 = topology.addAtom('H', Element.getBySymbol('H'), residue)
    H2 = topology.addAtom('H', Element.getBySymbol('H'), residue)
    topology.addBond(O, H1)
    topology.addBond(O, H2)

# Initial positions (ice-like)
positions = []
for i in range(10):
    x = (i % 3) * 0.275
    y = ((i // 3) % 3) * 0.275  
    z = (i // 9) * 0.275
    positions.append([x, y, z])  # O
    positions.append([x+0.0957, y, z])  # H1
    positions.append([x, y+0.0957, z])  # H2
positions = positions * unit.nanometer

# Box
box_size = 0.825 * unit.nanometer
topology.setPeriodicBoxVectors([[box_size, 0, 0], [0, box_size, 0], [0, 0, box_size]])

# Create system
print("\nCreating UMA potential...")
potential = MLPotential('uma-s-1p1-pythonforce-batch')
system = potential.createSystem(topology)

# Test different bead counts
bead_counts = [1, 2, 4, 8, 16]
results = {}

for num_beads in bead_counts:
    print(f"\n{'=' * 80}")
    print(f"Testing {num_beads} bead(s)")
    print(f"{'=' * 80}")
    
    # Create integrator
    from openmm import RPMDIntegrator
    integrator = RPMDIntegrator(
        numCopies=num_beads,
        temperature=243.0 * unit.kelvin,
        frictionCoeff=1.0 / unit.picoseconds,
        stepSize=1.0 * unit.femtoseconds
    )
    
    # Create simulation
    simulation = openmm.app.Simulation(topology, system, integrator, openmm.Platform.getPlatformByName('CUDA'))
    
    # Initialize positions for all beads
    for bead in range(num_beads):
        simulation.context.setPositions(positions, bead)
    simulation.context.setVelocitiesToTemperature(243.0 * unit.kelvin)
    
    # Warmup (3 steps)
    print(f"Warming up...")
    for _ in range(3):
        integrator.step(1)
    
    # Benchmark (20 steps)
    print(f"Benchmarking 20 steps...")
    times = []
    for i in range(20):
        start = time.time()
        integrator.step(1)
        end = time.time()
        times.append((end - start) * 1000)  # ms
    
    avg_time = np.mean(times)
    std_time = np.std(times)
    
    results[num_beads] = {
        'avg_ms': avg_time,
        'std_ms': std_time,
        'times': times
    }
    
    print(f"  Average: {avg_time:.1f} ± {std_time:.1f} ms/step")
    print(f"  Min: {np.min(times):.1f} ms, Max: {np.max(times):.1f} ms")
    
    del simulation
    del integrator
    torch.cuda.empty_cache()

# Analysis
print(f"\n{'=' * 80}")
print("SCALING ANALYSIS")
print(f"{'=' * 80}")

print(f"\n{'Beads':<8} {'Time (ms)':<12} {'Per-bead (ms)':<15} {'Speedup vs 1':<15} {'Efficiency':<12}")
print("-" * 70)

baseline = results[1]['avg_ms']
for num_beads in bead_counts:
    avg = results[num_beads]['avg_ms']
    per_bead = avg / num_beads
    speedup = (baseline * num_beads) / avg
    efficiency = speedup / num_beads * 100
    
    print(f"{num_beads:<8} {avg:<12.1f} {per_bead:<15.1f} {speedup:<15.2f}x {efficiency:<12.1f}%")

print(f"\n{'=' * 80}")
print("INTERPRETATION")
print(f"{'=' * 80}")

# Check if we're getting parallelization
efficiency_8 = (baseline * 8) / results[8]['avg_ms'] / 8 * 100

if efficiency_8 > 80:
    print("✅ EXCELLENT: Near-perfect parallel scaling (>80% efficiency)")
elif efficiency_8 > 50:
    print("✅ GOOD: Strong parallel scaling (50-80% efficiency)")
    print("   This is typical for GPU batch processing with GNNs")
elif efficiency_8 > 25:
    print("⚠️  MODERATE: Some parallelization but with overhead (25-50%)")
else:
    print("❌ POOR: Mostly sequential processing (<25% efficiency)")

print(f"\nParallel efficiency at 8 beads: {efficiency_8:.1f}%")
print(f"\nIf sequential (no batching): 8 beads would take {baseline * 8:.1f} ms")
print(f"Actual with batching: {results[8]['avg_ms']:.1f} ms")
print(f"Speedup from batching: {(baseline * 8) / results[8]['avg_ms']:.2f}x")
