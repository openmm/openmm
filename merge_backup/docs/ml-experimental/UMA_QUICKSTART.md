# UMA Quick Start Guide

This guide will help you get started with UMA models in OpenMM-ML quickly.

## Installation

```bash
# 1. Install OpenMM and OpenMM-Torch
conda install -c conda-forge openmm openmmtorch

# 2. Install FAIRChem
pip install fairchem-core

# 3. Install OpenMM-ML (from this repository)
cd fairchem/openmm-ml
pip install -e .
```

## 5-Minute Example: Basic UMA Simulation

```python
import openmm
from openmm import app, unit
from openmmml import MLPotential
import numpy as np

# Create a water molecule
topology = app.Topology()
chain = topology.addChain()
residue = topology.addResidue("WAT", chain)
O = topology.addAtom("O", app.Element.getBySymbol("O"), residue)
H1 = topology.addAtom("H1", app.Element.getBySymbol("H"), residue)
H2 = topology.addAtom("H2", app.Element.getBySymbol("H"), residue)
topology.addBond(O, H1)
topology.addBond(O, H2)

positions = np.array([
    [0.0, 0.0, 0.0],
    [0.0957, 0.0, 0.0],
    [-0.024, 0.093, 0.0]
]) * unit.nanometers

# Create UMA system
potential = MLPotential('uma-s-1p1')
system = potential.createSystem(topology, task_name='omol')

# Run simulation on GPU
platform = openmm.Platform.getPlatformByName('CUDA')
integrator = openmm.LangevinMiddleIntegrator(
    300*unit.kelvin, 1.0/unit.picosecond, 1.0*unit.femtoseconds
)
context = openmm.Context(system, integrator, platform)
context.setPositions(positions)
context.setVelocitiesToTemperature(300*unit.kelvin)

# Run 1000 steps
integrator.step(1000)

# Get results
state = context.getState(getEnergy=True, getPositions=True)
print(f"Final energy: {state.getPotentialEnergy()}")
```

## RPMD Example

```python
# After creating system...
integrator = openmm.RPMDIntegrator(
    8,  # number of beads
    300*unit.kelvin,
    1.0/unit.picosecond,
    0.5*unit.femtoseconds
)

# Use PILE_G thermostat
integrator.setThermostatType(openmm.RPMDIntegrator.PileG)

# Initialize beads
for bead in range(8):
    bead_pos = [pos + np.random.randn(3)*0.001*unit.nanometers 
                for pos in positions]
    integrator.setPositions(bead, bead_pos)

# Run RPMD
platform = openmm.Platform.getPlatformByName('CUDA')
context = openmm.Context(system, integrator, platform)
integrator.step(1000)
```

## Available Models

- `uma-s-1p1` - Small, fast (recommended for testing)
- `uma-m-1p1` - Medium, more accurate
- `esen-md-direct-all-omol` - Optimized for molecules

## Task Names

- `omol` - Organic molecules (most common)
- `omat` - Inorganic materials
- `oc20` - Catalysis systems

## Common Parameters

```python
system = potential.createSystem(
    topology,
    task_name='omol',           # Required
    charge=0,                   # Molecular charge
    spin=1,                     # Spin multiplicity
    inference_settings='turbo'  # 'default' or 'turbo'
)
```

## Troubleshooting

**Model download fails:**
```python
# Check your cache directory
from fairchem.core._config import CACHE_DIR
print(CACHE_DIR)  # Should show ~/.cache/fairchem
```

**GPU not used:**
```bash
# Verify OpenMM-Torch installation
conda list | grep openmmtorch

# Check PyTorch CUDA
python -c "import torch; print(torch.cuda.is_available())"
```

**Memory error with RPMD:**
- Reduce number of beads (try 4 or 8 instead of 32)
- Use smaller molecule
- Enable gradient checkpointing

## Next Steps

1. Try the examples in `fairchem/openmm-ml/examples/uma/`
2. Read the full documentation in `fairchem/openmm-ml/doc/uma.rst`
3. Run the validation script to verify your setup
4. Benchmark GPU performance with `benchmark_gpu.py`

## Getting Help

- OpenMM-ML issues: https://github.com/openmm/openmm-ml/issues
- FAIRChem documentation: https://fair-chem.github.io/
- OpenMM documentation: https://openmm.org/
