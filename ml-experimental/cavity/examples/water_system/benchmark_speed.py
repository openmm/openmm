#!/usr/bin/env python3
"""
Benchmark speed for cavity-coupled water system.
Tests performance on different platforms (Reference, CPU, CUDA).
"""

import time
import numpy as np
import openmm
import openmm.unit as unit
from openmm import app

def create_test_system(num_molecules=500):
    """Create a small water system for benchmarking."""
    from run_simulation import create_water_box, add_cavity_particle, setup_cavity_coupling
    
    # System parameters
    box_size_nm = 2.5
    temperature_K = 300.0
    cavity_freq_cm = 1600
    lambda_coupling = 0.01
    
    print(f"  Creating {num_molecules} water molecules...")
    system, topology, positions, charges, num_molecular = create_water_box(
        num_molecules=num_molecules,
        box_size_nm=box_size_nm,
        temperature_K=temperature_K
    )
    
    print("  Adding cavity particle...")
    cavity_index = add_cavity_particle(
        system, topology, positions, charges,
        cavity_freq_cm=cavity_freq_cm,
        box_size_nm=box_size_nm
    )
    
    print("  Setting up cavity coupling...")
    from run_simulation import wavenumber_to_hartree, AMU_TO_AU
    omegac_au = wavenumber_to_hartree(cavity_freq_cm)
    photon_mass_amu = 1.0 / AMU_TO_AU
    
    cavity_force, displacer = setup_cavity_coupling(
        system, cavity_index, omegac_au, 
        lambda_coupling, photon_mass_amu
    )
    
    return system, topology, positions

def benchmark_platform(platform_name, num_molecules=500, num_steps=1000):
    """Benchmark a specific platform."""
    print(f"\n{'='*70}")
    print(f"Benchmarking {platform_name} Platform")
    print(f"{'='*70}")
    
    try:
        # Create system
        print("Creating system...")
        system, topology, positions = create_test_system(num_molecules)
        
        # Create integrator
        integrator = openmm.LangevinMiddleIntegrator(
            300 * unit.kelvin,
            0.5 / unit.picosecond,
            0.0005 * unit.picosecond  # 0.5 fs timestep
        )
        
        # Create context
        print(f"Creating context on {platform_name}...")
        platform = openmm.Platform.getPlatformByName(platform_name)
        
        # Set platform properties for CUDA/OpenCL
        properties = {}
        if platform_name in ['CUDA', 'OpenCL']:
            properties['Precision'] = 'mixed'
        
        context = openmm.Context(system, integrator, platform, properties)
        context.setPositions(positions)
        
        # Minimize
        print("Minimizing energy...")
        openmm.LocalEnergyMinimizer.minimize(context, maxIterations=100)
        
        # Warmup
        print("Warming up...")
        integrator.step(100)
        
        # Benchmark
        print(f"Running {num_steps} steps...")
        start = time.time()
        integrator.step(num_steps)
        elapsed = time.time() - start
        
        # Calculate performance
        steps_per_sec = num_steps / elapsed
        ns_per_day = (steps_per_sec * 0.0005 * 86400) / 1000  # dt=0.5fs
        
        print(f"\n{'='*70}")
        print(f"Results for {platform_name}:")
        print(f"  Time: {elapsed:.2f} seconds")
        print(f"  Speed: {steps_per_sec:.1f} steps/s")
        print(f"  Performance: {ns_per_day:.2f} ns/day")
        print(f"{'='*70}")
        
        return ns_per_day
        
    except Exception as e:
        print(f"\nError with {platform_name}: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    print("""
    ╔══════════════════════════════════════════════════════════════════╗
    ║      Cavity-Coupled Water System - Speed Benchmark              ║
    ╚══════════════════════════════════════════════════════════════════╝
    """)
    
    num_molecules = 500
    num_steps = 1000
    
    print(f"System size: {num_molecules} water molecules ({num_molecules*3} atoms + 1 cavity)")
    print(f"Benchmark: {num_steps} MD steps (0.5 fs timestep)")
    
    # Test all available platforms
    results = {}
    for platform_name in ['Reference', 'CPU', 'CUDA', 'OpenCL']:
        try:
            openmm.Platform.getPlatformByName(platform_name)
            results[platform_name] = benchmark_platform(platform_name, num_molecules, num_steps)
        except:
            print(f"\n{platform_name} platform not available, skipping...")
    
    # Summary
    print(f"\n\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"{'Platform':<15} {'Performance (ns/day)':<20} {'Speedup':<10}")
    print(f"{'-'*70}")
    
    reference_speed = results.get('Reference', None)
    for platform, speed in sorted(results.items(), key=lambda x: x[1] or 0, reverse=True):
        if speed:
            speedup = f"{speed/reference_speed:.1f}x" if reference_speed else "N/A"
            print(f"{platform:<15} {speed:>8.2f}             {speedup:<10}")
    
    print(f"{'='*70}\n")
    
    # Recommendation
    best_platform = max(results.items(), key=lambda x: x[1] or 0)
    print(f"✓ Recommended platform: {best_platform[0]} ({best_platform[1]:.2f} ns/day)\n")

if __name__ == '__main__':
    main()
