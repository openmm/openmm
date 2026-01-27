#!/bin/bash
set -e  # Exit on error

echo "=================================="
echo "Installing OpenMM-ML with UMA"
echo "=================================="

# Navigate to openmm-ml directory
cd /media/extradrive/Trajectories/openmm/fairchem/openmm-ml

# Uninstall existing version (if any)
echo "Removing existing installation..."
pip uninstall -y openmmml || true

# Install in editable mode
echo "Installing OpenMM-ML in editable mode..."
pip install -e .

# Verify installation
echo -e "\nVerifying installation..."
python -c "
from openmmml import MLPotential
import sys

# Check if UMA models are registered
try:
    # Try to instantiate a UMA potential (won't download yet, just check registration)
    print('Testing UMA model registration...')
    potential = MLPotential('uma-s-1p1')
    print('✓ uma-s-1p1 registered successfully')
    
    # List all available models
    from openmmml.mlpotential import MLPotential
    available = list(MLPotential._implFactories.keys())
    uma_models = [m for m in available if m.startswith('uma') or m.startswith('esen')]
    
    print(f'\n✓ Found {len(uma_models)} UMA models:')
    for model in sorted(uma_models):
        print(f'  - {model}')
    
    print('\n✓ Installation successful!')
    sys.exit(0)
    
except KeyError as e:
    print(f'✗ UMA models not found in registry: {e}')
    print('Available models:', list(MLPotential._implFactories.keys()))
    sys.exit(1)
except Exception as e:
    print(f'✗ Error: {e}')
    import traceback
    traceback.print_exc()
    sys.exit(1)
"

echo -e "\n=================================="
echo "Installation complete!"
echo "=================================="
