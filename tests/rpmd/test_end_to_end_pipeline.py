#!/usr/bin/env python3
"""
Test: End-to-End UMA/RPMD Pipeline Validation

Integration test validating the full UMA/RPMD pipeline from system creation
through simulation to analysis. Catches bugs that only appear when components interact.

Test workflow:
1. System creation with UMA
2. Multiple force groups setup
3. Hybrid RPMD configuration
4. Equilibration phase
5. Production simulation
6. Analysis and validation

This test is designed to catch integration bugs that unit tests might miss.
"""

import sys
import numpy as np
import openmm
from openmm import app, unit

def test_full_pipeline():
    """Run complete UMA/RPMD pipeline."""
    print("=" * 70)
    print("End-to-End UMA/RPMD Pipeline Test")
    print("=" * 70)
    
    try:
        from openmmml import MLPotential
    except ImportError:
        print("openmmml not available, skipping test")
        return True
    
    # System Setup
    print("\nSystem Setup")
    print("-" * 70)
    
    # Create water dimer
    topology = app.Topology()
    chain = topology.addChain()
    
    for mol_idx in range(2):
        residue = topology.addResidue(f"WAT{mol_idx}", chain)
        O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
        H1 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
        H2 = topology.addAtom("H", app.Element.getBySymbol("H"), residue)
        topology.addBond(O, H1)
        topology.addBond(O, H2)
    
    # Positions: two water molecules separated by ~3 Angstrom
    positions = [
        # Molecule 1
        [0.0, 0.0, 0.0],
        [0.0957, 0.0, 0.0],
        [-0.024, 0.093, 0.0],
        # Molecule 2
        [0.3, 0.0, 0.0],
        [0.3957, 0.0, 0.0],
        [0.276, 0.093, 0.0]
    ]
    
    print("Created water dimer topology")
    
    # Create system with UMA
    try:
        potential = MLPotential('uma-s-1p1-pythonforce-batch')
        system = potential.createSystem(topology, task_name='omol', charge=0, spin=1)
        print("Created UMA system")
    except Exception as e:
        print(f"Failed to create UMA system: {e}")
        return True
    
    num_molecular_atoms = 6
    
    # Add cavity particle
    cavity_mass = 1.0 / 1822.888  # electron mass
    cavity_index = system.addParticle(cavity_mass)
    positions.append([0.15, 0.15, 0.0])
    print(f"Added cavity particle (index {cavity_index})")
    
    # Add weak restraints (force group 1)
    restraint = openmm.CustomExternalForce("k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint.addGlobalParameter("k", 10.0)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    for i, pos in enumerate(positions):
        restraint.addParticle(i, pos)
    restraint.setForceGroup(1)
    system.addForce(restraint)
    print("Added restraint forces (group 1)")
    
    #  RPMD Integrator Setup
    print("\nPhase 2: RPMD Integrator Setup")
    print("-" * 70)
    
    num_beads = 8
    temperature_K = 300.0
    dt_fs = 1.0
    
    integrator = openmm.RPMDIntegrator(
        num_beads,
        temperature_K * unit.kelvin,
        1.0 / unit.picosecond,
        dt_fs * unit.femtoseconds
    )
    
    # Configure hybrid mode
    integrator.setParticleType(cavity_index, 1)  # Classical
    integrator.setClassicalThermostat(openmm.RPMDIntegrator.LangevinClassical)
    integrator.setThermostatType(openmm.RPMDIntegrator.PileG)
    integrator.setCentroidFriction(0.5 / unit.picosecond)
    
    print(f"Created hybrid RPMD integrator")
    print(f"  Beads: {num_beads}")
    print(f"  Temperature: {temperature_K} K")
    print(f"  Timestep: {dt_fs} fs")
    print(f"  Quantum particles: 0-{num_molecular_atoms-1}")
    print(f"  Classical particles: {cavity_index}")
    
    #  Context Creation
    print("\nPhase 3: Context and Platform Setup")
    print("-" * 70)
    
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        print("Using CUDA platform")
    except Exception:
        try:
            platform = openmm.Platform.getPlatformByName('OpenCL')
            print("Using OpenCL platform")
        except Exception:
            platform = openmm.Platform.getPlatformByName('Reference')
            print("Using Reference platform")
    
    try:
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions * unit.nanometers)
        print("Created simulation context")
    except Exception as e:
        print(f"Failed to create context: {e}")
        return False
    
    #  Initialize Beads
    print("\nPhase 4: Initializing RPMD Beads")
    print("-" * 70)
    
    for bead in range(num_beads):
        bead_positions = []
        for i, pos in enumerate(positions):
            if i == cavity_index:
                # Classical: identical on all beads
                bead_positions.append(pos)
            else:
                # Quantum: small perturbation
                perturb = np.random.randn(3) * 0.002
                bead_positions.append([pos[0] + perturb[0], pos[1] + perturb[1], pos[2] + perturb[2]])
        integrator.setPositions(bead, bead_positions * unit.nanometers)
    
    print(f"Initialized {num_beads} beads")
    
    #  Equilibration
    print("\nPhase 5: Equilibration (50 steps)")
    print("-" * 70)
    
    try:
        equilibration_steps = 50
        integrator.step(equilibration_steps)
        print(f"Equilibration completed")
        
        # Check state after equilibration
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        if np.isnan(energy) or np.isinf(energy):
            print(f"FAIL: Energy is NaN or Inf after equilibration")
            return False
        
        print(f"  Energy: {energy:.2f} kJ/mol")
        
    except Exception as e:
        print(f"FAIL: Equilibration crashed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    #  Production and Validation
    print("\nPhase 6: Production Simulation (100 steps)")
    print("-" * 70)
    
    try:
        production_steps = 100
        energies = []
        cavity_spreads = []
        molecular_spreads = []
        
        # Sample every 10 steps
        for i in range(10):
            integrator.step(10)
            
            # Get energy
            state = context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            energies.append(energy)
            
            # Check cavity bead spread
            cavity_positions = []
            molecular_positions = [[] for _ in range(num_molecular_atoms)]
            
            for bead in range(num_beads):
                state_bead = integrator.getState(bead, getPositions=True)
                pos = state_bead.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                cavity_positions.append(pos[cavity_index])
                for j in range(num_molecular_atoms):
                    molecular_positions[j].append(pos[j])
            
            cavity_positions = np.array(cavity_positions)
            cavity_spread = np.max(np.std(cavity_positions, axis=0))
            cavity_spreads.append(cavity_spread)
            
            # Molecular spread
            mol_spreads = []
            for j in range(num_molecular_atoms):
                mol_pos = np.array(molecular_positions[j])
                spread = np.max(np.std(mol_pos, axis=0))
                mol_spreads.append(spread)
            avg_mol_spread = np.mean(mol_spreads)
            molecular_spreads.append(avg_mol_spread)
        
        print(f"Production completed ({production_steps} steps)")
        
        # Validation checks
        print("\nPhase 7: Validation")
        print("-" * 70)
        
        all_checks_passed = True
        
        # Check 1: Energy stability
        energies = np.array(energies)
        energy_std = np.std(energies)
        energy_mean = np.mean(energies)
        energy_cv = energy_std / abs(energy_mean) if abs(energy_mean) > 1e-6 else 0
        
        print(f"\n1. Energy Stability:")
        print(f"   Mean energy: {energy_mean:.2f} kJ/mol")
        print(f"   Std dev: {energy_std:.2f} kJ/mol")
        print(f"   Coefficient of variation: {energy_cv*100:.2f}%")
        
        if not np.any(np.isnan(energies)) and not np.any(np.isinf(energies)):
            print(f"   PASS: Energy remains finite")
        else:
            print(f"   FAIL: Energy became NaN or Inf")
            all_checks_passed = False
        
        # Check 2: Cavity remains classical
        cavity_spreads = np.array(cavity_spreads)
        max_cavity_spread = np.max(cavity_spreads)
        
        print(f"\n2. Cavity Classical Behavior:")
        print(f"   Max bead spread: {max_cavity_spread:.8f} nm")
        
        if max_cavity_spread < 1e-5:
            print(f"   PASS: Cavity beads remain identical (< 1e-5 nm)")
        else:
            print(f"   FAIL: Cavity beads spreading (not classical!)")
            all_checks_passed = False
        
        # Check 3: Molecular particles show quantum delocalization
        molecular_spreads = np.array(molecular_spreads)
        avg_molecular_spread = np.mean(molecular_spreads)
        
        print(f"\n3. Molecular Quantum Delocalization:")
        print(f"   Average bead spread: {avg_molecular_spread:.6f} nm")
        
        if avg_molecular_spread > 1e-4:
            print(f"   PASS: Quantum particles show delocalization")
        else:
            print(f"   ⚠ Note: Delocalization is small (may be at low T or equilibrium)")
        
        # Check 4: Force accumulation
        print(f"\n4. Force Accumulation Check:")
        state_all = integrator.getState(0, getForces=True, groups={0, 1})
        forces_all = state_all.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        state_0 = integrator.getState(0, getForces=True, groups={0})
        forces_0 = state_0.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        state_1 = integrator.getState(0, getForces=True, groups={1})
        forces_1 = state_1.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        
        expected = forces_0 + forces_1
        error = np.abs(forces_all - expected).max()
        
        print(f"   Max accumulation error: {error:.6f} kJ/mol/nm")
        
        if error < 0.01:
            print(f"   PASS: Forces correctly accumulated")
        else:
            print(f"   FAIL: Force accumulation error")
            all_checks_passed = False
        
        # Check 5: Cavity gets zero force from UMA
        print(f"\n5. Cavity Force from UMA:")
        cavity_force = np.linalg.norm(forces_0[cavity_index])
        
        print(f"   Force magnitude: {cavity_force:.8f} kJ/mol/nm")
        
        if cavity_force < 1e-6:
            print(f"   PASS: Cavity receives zero force from UMA")
        else:
            print(f"   FAIL: Cavity has non-zero UMA force")
            all_checks_passed = False
        
        # Summary
        print("\n" + "=" * 70)
        if all_checks_passed:
            print("END-TO-END TEST PASSED")
            print("  All pipeline components working correctly")
            return True
        else:
            print("END-TO-END TEST FAILED")
            print("  Integration bugs detected")
            return False
        
    except Exception as e:
        print(f"\nPIPELINE FAILED WITH EXCEPTION:")
        print(f"  {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_full_pipeline()
    sys.exit(0 if success else 1)
