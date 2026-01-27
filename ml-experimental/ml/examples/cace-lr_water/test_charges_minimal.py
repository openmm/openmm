#!/usr/bin/env python3
"""Test eSEN charge prediction - direct fairchem API."""

import sys
sys.path.insert(0, '/media/extradrive/Trajectories/openmm/fairchem-stable/src')

import numpy as np
import torch
from ase.build import molecule
from fairchem.core import FAIRChemCalculator
from fairchem.core.calculate import pretrained_mlip

def test_model_charges(model_name):
    """Test charge extraction from a specific model."""
    
    print(f"\n{'='*70}")
    print(f"Testing: {model_name}")
    print('='*70)
    
    # Create water molecule
    water = molecule('H2O')
    water.info['charge'] = 0
    water.info['spin'] = 1
    
    print(f"\nWater molecule:")
    for i, (sym, pos) in enumerate(zip(water.get_chemical_symbols(), water.get_positions())):
        print(f"  {i}: {sym:2s} [{pos[0]:7.3f}, {pos[1]:7.3f}, {pos[2]:7.3f}]")
    
    try:
        # Load model
        print(f"\nLoading model...")
        predictor = pretrained_mlip.get_predict_unit(model_name, device='cpu')  # Use CPU to avoid OOM
        print(f"  ✓ Model loaded")
        
        # Check task support
        if "omol" not in predictor.dataset_to_tasks:
            print(f"  ✗ No 'omol' task support")
            return False
        
        # Create calculator
        calc = FAIRChemCalculator(predictor, task_name="omol")
        water.calc = calc
        
        # Get energy
        print(f"\nComputing energy...")
        energy = water.get_potential_energy()
        print(f"  ✓ Energy: {energy:.4f} eV")
        
        # Check model structure for charges
        print(f"\nChecking model structure...")
        model = predictor.model
        print(f"  Model type: {type(model).__name__}")
        
        if hasattr(model, 'heads'):
            print(f"  ✓ Has {len(model.heads)} head(s)")
            for i, head in enumerate(model.heads):
                print(f"    Head {i}: {type(head).__name__}")
                
                # Check for charge-related methods
                charge_methods = [m for m in dir(head) if 'charge' in m.lower() and not m.startswith('_')]
                if charge_methods:
                    print(f"      Charge methods: {charge_methods}")
                    
                if hasattr(head, 'get_charges'):
                    print(f"      ✓✓✓ HAS get_charges()! ✓✓✓")
                    return True
                elif hasattr(head, 'get_lr_energies'):
                    print(f"      ✓ Has get_lr_energies() - may support charges")
                    return True
                    
        print(f"  ✗ No charge prediction capability found")
        return False
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("=" * 70)
    print("ESEN/UMA Charge Prediction Test")
    print("=" * 70)
    
    models_to_test = [
        "uma-s-1p1",
        "esen-md-direct-all-omol",
        "esen-sm-conserving-all-omol"
    ]
    
    for model_name in models_to_test:
        result = test_model_charges(model_name)
        if result:
            print(f"\n✓✓✓ FOUND WORKING MODEL: {model_name} ✓✓✓")
            break
    else:
        print(f"\n✗ No models with charge prediction found")
