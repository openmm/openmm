#!/usr/bin/env python3
"""
Benchmark speed for UMA ML potential water system.
Tests performance with different system sizes.
"""

import time
import numpy as np
import openmm
import openmm.unit as unit
from openmm.app import Topology, Element


def create_water_topology(num_molecules):
    """Create a water box topology."""
    oh_bond = 0.09572  # nm
    hoh_angle = 104.52 * np.pi / 180.0
    
    # Calculate box size
    vol_per_mol_nm3 = 0.0299
    total_vol = num_molecules * vol_per_mol_nm3
    box_size_nm = max((total_vol ** (1/3)) * 1.2, 1.0)
    
    topology = Topology()
    positions = []
    
    side = int(np.ceil(num_molecules**(1/3)))
    spacing = box_size_nm / side
    
    mol_count = 0
    np.random.seed(42)
    
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
    
    return topology, positions


def benchmark_uma(num_molecules=10, num_steps=100):
    """Benchmark UMA potential."""
    print(f"\n{'='*70}")
    print(f"Benchmarking UMA with {num_molecules} water molecules")
    print(f"{'='*70}")
    
    try:
        from openmmml import MLPotential
        import torch
        
        print(f"PyTorch: {torch.__version__}")
        print(f"CUDA available: {torch.cuda.is_available()}")
        if torch.cuda.is_available():
            print(f"CUDA device: {torch.cuda.get_device_name(0)}")
        
        # Create topology
        print(f"\nCreating {num_molecules} water molecules...")
        topology, positions = create_water_topology(num_molecules)
        print(f"  {len(positions)} atoms")
        
        # Create UMA potential
        print(f"\nLoading UMA model...")
        load_start = time.time()
        potential = MLPotential('uma-s-1p1-pythonforce')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
        load_time = time.time() - load_start
        print(f"  Model loaded in {load_time:.2f}s")
        
        # Create integrator and context
        print(f"\nCreating simulation context...")
        integrator = openmm.LangevinMiddleIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picosecond,
            0.0005 * unit.picosecond  # 0.5 fs
        )
        
        platform = openmm.Platform.getPlatformByName('CUDA')
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions * unit.nanometer)
        
        # Warmup
        print(f"Warming up (10 steps)...")
        integrator.step(10)
        
        # Benchmark
        print(f"Running {num_steps} steps...")
        start = time.time()
        integrator.step(num_steps)
        elapsed = time.time() - start
        
        # Calculate performance
        steps_per_sec = num_steps / elapsed
        ns_per_day = (steps_per_sec * 0.0005 * 86400) / 1000  # dt=0.5fs
        
        # Get energy to verify
        state = context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        
        print(f"\n{'='*70}")
        print(f"Results for UMA ({num_molecules} waters, {num_molecules*3} atoms):")
        print(f"  Time: {elapsed:.2f} seconds")
        print(f"  Speed: {steps_per_sec:.1f} steps/s")
        print(f"  Performance: {ns_per_day:.4f} ns/day")
        print(f"  Final PE: {pe:.1f} kJ/mol")
        print(f"{'='*70}")
        
        return {
            'num_molecules': num_molecules,
            'num_atoms': num_molecules * 3,
            'steps_per_sec': steps_per_sec,
            'ns_per_day': ns_per_day,
            'elapsed': elapsed
        }
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    print("""
    ╔══════════════════════════════════════════════════════════════════╗
    ║        UMA ML Potential Water System - Speed Benchmark           ║
    ╚══════════════════════════════════════════════════════════════════╝
    """)
    
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark UMA water system')
    parser.add_argument('--molecules', type=int, nargs='+', default=[3, 10, 27, 64],
                       help='Number of molecules to test (default: 3 10 27 64)')
    parser.add_argument('--steps', type=int, default=100,
                       help='Number of MD steps (default: 100)')
    args = parser.parse_args()
    
    results = []
    for n in args.molecules:
        result = benchmark_uma(n, args.steps)
        if result:
            results.append(result)
    
    # Summary
    if results:
        print(f"\n\n{'='*70}")
        print("SUMMARY")
        print(f"{'='*70}")
        print(f"{'Waters':<10} {'Atoms':<10} {'Steps/s':<15} {'ns/day':<15}")
        print(f"{'-'*70}")
        
        for r in results:
            print(f"{r['num_molecules']:<10} {r['num_atoms']:<10} "
                  f"{r['steps_per_sec']:<15.1f} {r['ns_per_day']:<15.4f}")
        
        print(f"{'='*70}\n")
        
        # Scaling analysis
        if len(results) >= 2:
            print("Scaling analysis:")
            for i in range(1, len(results)):
                r0, r1 = results[0], results[i]
                atom_ratio = r1['num_atoms'] / r0['num_atoms']
                speed_ratio = r0['steps_per_sec'] / r1['steps_per_sec']
                print(f"  {r0['num_atoms']} → {r1['num_atoms']} atoms: "
                      f"{atom_ratio:.1f}x atoms, {speed_ratio:.1f}x slower")


if __name__ == '__main__':
    main()
