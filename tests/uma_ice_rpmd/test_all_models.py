#!/usr/bin/env python3
"""
Test script to verify all UMA/eSEN models with pythonforce-batch variants are registered.
"""

from openmmml import MLPotential

# List of all models that should be registered
models = [
    # UMA models
    'uma-s-1-pythonforce-batch',
    'uma-s-1p1-pythonforce-batch',
    'uma-m-1p1-pythonforce-batch',
    
    # OMol25 eSEN models
    'esen-md-direct-all-omol-pythonforce-batch',
    'esen-sm-conserving-all-omol-pythonforce-batch',
    'esen-sm-direct-all-omol-pythonforce-batch',
    
    # OC25 eSEN models
    'esen-sm-conserving-all-oc25-pythonforce-batch',
    'esen-md-direct-all-oc25-pythonforce-batch',
]

print("=" * 80)
print("Testing UMA/eSEN Model Registration (pythonforce-batch variants)")
print("=" * 80)
print()

passed = 0
failed = 0
failed_models = []

for model_name in models:
    try:
        # Try to create the potential (this will check if it's registered)
        # We don't actually load the model yet, just test registration
        print(f"Testing: {model_name:<50}", end=" ")
        potential = MLPotential(model_name)
        print("✓ REGISTERED")
        passed += 1
    except Exception as e:
        print(f"✗ FAILED: {e}")
        failed += 1
        failed_models.append(model_name)

print()
print("=" * 80)
print(f"Results: {passed}/{len(models)} models registered successfully")
if failed > 0:
    print(f"Failed models: {', '.join(failed_models)}")
print("=" * 80)
