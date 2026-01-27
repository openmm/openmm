#!/usr/bin/env python3
"""Simple test for eSEN charge prediction on water - uses existing UMA wrapper."""

import sys
import numpy as np

try:
    from openmmml import MLPotential
    from ase.build import molecule
    import torch
    print("✓ Imports successful")
except ImportError as e:
    print(f"✗ Import error: {e}")
    sys.exit(1)

print("\n" + "=" * 70)
print("eSEN Charge Prediction Test - Using UMA Wrapper")
print("=" * 70)

# Create water molecule
water = molecule('H2O')
print(f"\nWater molecule created:")
print(f"  Symbols: {water.get_chemical_symbols()}")
print(f"  Positions (Å):")
for i, (sym, pos) in enumerate(zip(water.get_chemical_symbols(), water.get_positions())):
    print(f"    {i}: {sym:2s} [{pos[0]:7.3f}, {pos[1]:7.3f}, {pos[2]:7.3f}]")

# Try UMA models through openmmml
print("\n" + "=" * 70)
print("Testing UMA models via openmmml...")
print("=" * 70)

model_names = ['uma-s-1p1-pythonforce']

for model_name in model_names:
    print(f"\nTrying {model_name}...")
    try:
        # Load potential
        potential = MLPotential(model_name)
        print(f"  ✓ Model loaded")
        
        # Check if we can access the underlying predictor
        if hasattr(potential, 'predictor'):
            predictor = potential.predictor
            print(f"  ✓ Has predictor attribute")
            
            # Check model structure
            if hasattr(predictor, 'model'):
                model = predictor.model
                print(f"  ✓ Has model attribute")
                print(f"    Model type: {type(model)}")
                
                # Check for heads
                if hasattr(model, 'heads'):
                    print(f"  ✓ Model has heads: {len(model.heads)} head(s)")
                    head = model.heads[0]
                    print(f"    Head type: {type(head).__name__}")
                    
                    # Check for charge methods
                    if hasattr(head, 'get_charges'):
                        print(f"  ✓✓✓ HEAD HAS get_charges() METHOD! ✓✓✓")
                        print(f"\n  This model supports charge prediction!")
                        print(f"  Model: {model_name}")
                    else:
                        print(f"  ✗ Head does not have get_charges()")
                        print(f"    Available methods: {[m for m in dir(head) if not m.startswith('_') and 'charge' in m.lower()]}")
                else:
                    print(f"  ✗ Model has no 'heads' attribute")
                    
        else:
            print(f"  ✗ No predictor attribute")
            
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()

print("\n" + "=" * 70)
print("Test complete - check results above")
print("=" * 70)
