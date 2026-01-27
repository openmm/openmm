#!/usr/bin/env python3
"""Investigate model structure to find charges."""

import sys
sys.path.insert(0, '/media/extradrive/Trajectories/openmm/fairchem-stable/src')

from fairchem.core.calculate import pretrained_mlip

print("=" * 70)
print("Deep Model Structure Investigation")
print("=" * 70)

model_name = "esen-md-direct-all-omol"

print(f"\nLoading {model_name}...")
predictor = pretrained_mlip.get_predict_unit(model_name, device='cpu')
print(f"✓ Loaded")

model = predictor.model

def inspect_object(obj, name="object", indent=0):
    """Recursively inspect an object."""
    prefix = "  " * indent
    print(f"{prefix}{name}: {type(obj).__name__}")
    
    # Check for common model attributes
    for attr_name in ['module', 'model', 'modules', 'heads', 'head', 'backbone', 
                      'get_charges', 'get_lr_energies', '_modules']:
        if hasattr(obj, attr_name):
            attr = getattr(obj, attr_name)
            if callable(attr) and not isinstance(attr, type):
                print(f"{prefix}  ✓ Method: {attr_name}()")
            elif attr_name == '_modules':
                # PyTorch module dict
                print(f"{prefix}  ✓ _modules: {list(attr.keys())}")
                for key, submodule in list(attr.items())[:3]:  # First 3
                    inspect_object(submodule, f"_modules['{key}']", indent+2)
            else:
                try:
                    if hasattr(attr, '__len__') and not isinstance(attr, str):
                        print(f"{prefix}  ✓ Attribute: {attr_name} (len={len(attr)})")
                        if len(attr) > 0 and indent < 3:
                            inspect_object(attr[0] if hasattr(attr, '__getitem__') else list(attr)[0], 
                                         f"{attr_name}[0]", indent+2)
                    else:
                        inspect_object(attr, attr_name, indent+1)
                except:
                    print(f"{prefix}  ✓ Attribute: {attr_name}")

inspect_object(model, "predictor.model", 0)

# Try to access inner model
print("\n" + "=" * 70)
print("Attempting to access inner model...")
print("=" * 70)

if hasattr(model, 'module'):
    inner = model.module
    print(f"\n✓ Found model.module: {type(inner).__name__}")
    inspect_object(inner, "model.module", 0)
    
    if hasattr(inner, 'heads'):
        print(f"\n✓✓ Found model.module.heads!")
        for i, head in enumerate(inner.heads):
            print(f"\n  Head {i}: {type(head).__name__}")
            if hasattr(head, 'get_charges'):
                print(f"    ✓✓✓ HAS get_charges()!")
            if hasattr(head, 'get_lr_energies'):
                print(f"    ✓✓✓ HAS get_lr_energies()!")
