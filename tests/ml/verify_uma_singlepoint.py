#!/usr/bin/env python3
"""
Single-point calculation verification for UMA integration.

This script:
1. Creates a water molecule
2. Performs single-point calculation with UMA (OpenMM-ML)
3. Performs reference calculation with FAIRChem ASE calculator
4. Compares energies and forces
5. Reports validation results
"""

import sys
import openmm
from openmm import app, unit
import numpy as np

def create_water_molecule():
    """Create a water molecule topology and positions."""
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("WAT", chain)
    O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
    H1 = topology.addAtom("H1", app.Element.getBySymbol("H"), residue)
    H2 = topology.addAtom("H2", app.Element.getBySymbol("H"), residue)
    topology.addBond(O, H1)
    topology.addBond(O, H2)
    
    # Equilibrium water geometry (in nm)
    positions = np.array([
        [0.0, 0.0, 0.0],      # O
        [0.0957, 0.0, 0.0],   # H1
        [-0.024, 0.093, 0.0]  # H2
    ])
    
    return topology, positions

def openmm_calculation(topology, positions, model_name='uma-s-1p1'):
    """Perform single-point calculation with OpenMM-ML UMA."""
    print("\n" + "="*70)
    print(f"OpenMM-ML Calculation ({model_name})")
    print("="*70)
    
    try:
        from openmmml import MLPotential
        
        # Create UMA potential
        print(f"Loading model: {model_name}")
        potential = MLPotential(model_name)
        
        # Create system
        print("Creating OpenMM system...")
        system = potential.createSystem(
            topology,
            task_name='omol',
            charge=0,
            spin=1
        )
        
        # Create context
        platform = openmm.Platform.getPlatform(0)  # Use first available
        print(f"Using platform: {platform.getName()}")
        
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions * unit.nanometers)
        
        # Get energy and forces
        state = context.getState(getEnergy=True, getForces=True)
        energy_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        forces = state.getForces(asNumpy=True)
        
        print(f"\nResults:")
        print(f"  Energy: {energy_kj:.6f} kJ/mol")
        print(f"  Forces (O):  [{forces[0][0]:.6f}, {forces[0][1]:.6f}, {forces[0][2]:.6f}] kJ/mol/nm")
        print(f"  Forces (H1): [{forces[1][0]:.6f}, {forces[1][1]:.6f}, {forces[1][2]:.6f}] kJ/mol/nm")
        print(f"  Forces (H2): [{forces[2][0]:.6f}, {forces[2][1]:.6f}, {forces[2][2]:.6f}] kJ/mol/nm")
        print(f"  Force magnitude: {np.linalg.norm(forces):.6f} kJ/mol/nm")
        
        return {
            'energy_kj': energy_kj,
            'forces': forces,
            'success': True
        }
        
    except Exception as e:
        print(f"\n OpenMM calculation failed: {e}")
        import traceback
        traceback.print_exc()
        return {'success': False, 'error': str(e)}

def fairchem_calculation(positions, model_name='uma-s-1p1'):
    """Perform reference calculation with FAIRChem ASE calculator."""
    print("\n" + "="*70)
    print("FAIRChem ASE Calculator (Reference)")
    print("="*70)
    
    try:
        from fairchem.core.calculate import pretrained_mlip
        from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
        from ase import Atoms
        
        # Convert positions to Angstrom
        positions_ang = positions * 10.0
        
        # Load predict unit
        print(f"Loading FAIRChem model: {model_name}")
        predict_unit = pretrained_mlip.get_predict_unit(
            model_name,
            inference_settings='default',
            device='cpu'
        )
        
        # Create ASE calculator
        calc = FAIRChemCalculator(
            predict_unit=predict_unit,
            task_name='omol'
        )
        
        # Create ASE Atoms
        atoms = Atoms(
            symbols=['O', 'H', 'H'],
            positions=positions_ang,
            pbc=False
        )
        atoms.info['charge'] = 0
        atoms.info['spin'] = 1
        atoms.calc = calc
        
        # Calculate
        print("Computing energy and forces...")
        energy_ev = atoms.get_potential_energy()
        forces_ev_ang = atoms.get_forces()
        
        # Convert to kJ/mol
        energy_kj = energy_ev * 96.4853
        
        print(f"\nResults:")
        print(f"  Energy: {energy_kj:.6f} kJ/mol ({energy_ev:.6f} eV)")
        print(f"  Forces (O):  [{forces_ev_ang[0][0]:.6f}, {forces_ev_ang[0][1]:.6f}, {forces_ev_ang[0][2]:.6f}] eV/Å")
        print(f"  Forces (H1): [{forces_ev_ang[1][0]:.6f}, {forces_ev_ang[1][1]:.6f}, {forces_ev_ang[1][2]:.6f}] eV/Å")
        print(f"  Forces (H2): [{forces_ev_ang[2][0]:.6f}, {forces_ev_ang[2][1]:.6f}, {forces_ev_ang[2][2]:.6f}] eV/Å")
        
        return {
            'energy_kj': energy_kj,
            'energy_ev': energy_ev,
            'forces_ev_ang': forces_ev_ang,
            'success': True
        }
        
    except Exception as e:
        print(f"\n FAIRChem calculation failed: {e}")
        import traceback
        traceback.print_exc()
        return {'success': False, 'error': str(e)}

def compare_results(openmm_result, fairchem_result):
    """Compare OpenMM and FAIRChem results."""
    print("\n" + "="*70)
    print("Comparison & Validation")
    print("="*70)
    
    if not openmm_result['success']:
        print(" OpenMM calculation failed, cannot compare")
        return False
    
    if not fairchem_result['success']:
        print(" FAIRChem calculation failed, cannot compare")
        return False
    
    # Energy comparison
    energy_diff = abs(openmm_result['energy_kj'] - fairchem_result['energy_kj'])
    energy_rel_diff = energy_diff / abs(fairchem_result['energy_kj']) * 100
    
    print(f"\nEnergy Comparison:")
    print(f"  OpenMM:   {openmm_result['energy_kj']:.6f} kJ/mol")
    print(f"  FAIRChem: {fairchem_result['energy_kj']:.6f} kJ/mol")
    print(f"  Absolute difference: {energy_diff:.6f} kJ/mol")
    print(f"  Relative difference: {energy_rel_diff:.4f} %")
    
    # Force comparison
    # Convert OpenMM forces (kJ/mol/nm) to FAIRChem forces (eV/Å)
    # 1 kJ/mol/nm = 0.0103642 eV/Å
    forces_openmm_ev_ang = openmm_result['forces'] * 0.0103642
    forces_fairchem = fairchem_result['forces_ev_ang']
    
    force_diff = np.linalg.norm(forces_openmm_ev_ang - forces_fairchem)
    force_magnitude = np.linalg.norm(forces_fairchem)
    force_rel_diff = force_diff / force_magnitude * 100
    
    print(f"\nForce Comparison:")
    print(f"  Force difference magnitude: {force_diff:.6f} eV/Å")
    print(f"  Reference force magnitude: {force_magnitude:.6f} eV/Å")
    print(f"  Relative difference: {force_rel_diff:.4f} %")
    
    # Validation
    energy_threshold = 0.1  # 0.1% relative difference
    force_threshold = 5.0   # 5% relative difference
    
    energy_pass = energy_rel_diff < energy_threshold
    force_pass = force_rel_diff < force_threshold
    
    print(f"\nValidation Results:")
    print(f"  Energy: {'PASS' if energy_pass else ' FAIL'} (threshold: {energy_threshold}%)")
    print(f"  Forces: {'PASS' if force_pass else ' FAIL'} (threshold: {force_threshold}%)")
    
    if energy_pass and force_pass:
        print(f"\n{'='*70}")
        print("ALL VALIDATIONS PASSED")
        print(f"{'='*70}")
        return True
    else:
        print(f"\n{'='*70}")
        print(" VALIDATION FAILED")
        print(f"{'='*70}")
        if not energy_pass:
            print(f"Energy difference ({energy_rel_diff:.4f}%) exceeds threshold ({energy_threshold}%)")
        if not force_pass:
            print(f"Force difference ({force_rel_diff:.4f}%) exceeds threshold ({force_threshold}%)")
        return False

def main():
    print("="*70)
    print("UMA Single-Point Calculation Verification")
    print("="*70)
    
    # Create water molecule
    print("\nCreating water molecule...")
    topology, positions = create_water_molecule()
    print("  3 atoms: O, H, H")
    print(f"  Positions (nm):")
    print(f"    O:  [{positions[0][0]:.4f}, {positions[0][1]:.4f}, {positions[0][2]:.4f}]")
    print(f"    H1: [{positions[1][0]:.4f}, {positions[1][1]:.4f}, {positions[1][2]:.4f}]")
    print(f"    H2: [{positions[2][0]:.4f}, {positions[2][1]:.4f}, {positions[2][2]:.4f}]")
    
    # Run calculations
    openmm_result = openmm_calculation(topology, positions)
    fairchem_result = fairchem_calculation(positions)
    
    # Compare results
    success = compare_results(openmm_result, fairchem_result)
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
