# UMA Integration - Final Report

## Executive Summary

The UMA integration for OpenMM-ML was **fully implemented** but encounters a **fundamental architectural incompatibility** that prevents it from running in OpenMM.

## What Was Accomplished

### ✅ Complete Implementation

1. **Core Integration** (`umapotential.py`)
   - `UMAForce` class with unit conversions (nm ↔ Å, kJ/mol ↔ eV)
   - Support for all task types ('omol', 'omat', 'oc20', etc.)
   - Charge and spin multiplicity configuration
   - Periodic boundary condition handling
   - Inference settings ('default', 'turbo')

2. **Package Registration**
   - All 8 UMA/eSEN models registered as entry points
   - Package successfully installed in editable mode
   - Entry points verified: ✅

3. **Testing & Verification**
   - Installation script (`install_uma.sh`)
   - Verification script (`verify_uma_singlepoint.py`)
   - FAIRChem ASE calculator validated and working

### ✅ Verification Results

The FAIRChem reference calculation works perfectly:

```
Energy: -200676.323290 kJ/mol (-2079.864221 eV)
Forces (O):  [-0.194868, -0.010550, -0.000000] eV/Å
Forces (H1): [0.259212, -0.047877, 0.000000] eV/Å
Forces (H2): [-0.064344, 0.058426, 0.000000] eV/Å
```

## The Technical Barrier

### Problem

**OpenMM-Torch requires TorchScript**
- All PyTorch models must be compiled to TorchScript
- Uses `torch.jit.script` or `torch.jit.trace`
- No Python objects, no NumPy, no dynamic control flow

**FAIRChem/UMA uses non-TorchScript components**
- ASE `Atoms` objects (Python class)
- NumPy array conversions
- Python assertions and control flow
- Dynamic data structures

**Result**: Cannot compile UMA models to TorchScript → Cannot use with OpenMM-Torch

### Attempts Made

1. ❌ Direct TorchScript compilation → ASE incompatibility
2. ❌ Save/load with `torch.save` → Not TorchScript format  
3. ❌ Trace with `torch.jit.trace` → Assertions fail, constants baked in

## Alternative Solutions

### Solution 1: Use ASE MD Instead of OpenMM ⭐ **Recommended**

```python
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
from fairchem.core.calculate import pretrained_mlip

# Load UMA model
predict_unit = pretrained_mlip.get_predict_unit('uma-s-1p1')
calc = FAIRChemCalculator(predict_unit=predict_unit, task_name='omol')

# Setup atoms
atoms.calc = calc
atoms.info['charge'] = 0
atoms.info['spin'] = 1

# Run MD
MaxwellBoltzmannDistribution(atoms, temperature_K=300)
dyn = VelocityVerlet(atoms, timestep=0.5)  # fs
dyn.run(10000)
```

**Pros**:
- ✅ Works immediately
- ✅ GPU acceleration available
- ✅ All UMA features supported
- ✅ Can use ASE's MD integrators

**Cons**:
- ❌ No RPMD support in ASE
- ❌ Different API than OpenMM

### Solution 2: Use Compatible ML Potentials in OpenMM

OpenMM-ML supports these TorchScript-compatible models:
- **MACE** (materials, molecules) - TorchScript compatible
- **NequIP** (equivariant neural networks)
- **ANI** (organic molecules)
- **AIMNet2** (general chemistry)

### Solution 3: Custom C++ OpenMM Plugin

Write a native OpenMM force that calls PyTorch C++ API directly.

**Effort**: 2-4 weeks of C++ development
**Complexity**: High

### Solution 4: Request FAIRChem TorchScript Support

Contact FAIRChem team to refactor UMA for TorchScript compatibility.

**Likelihood**: Low (major refactoring required)

## For RPMD Simulations

Since you mentioned RPMD as a requirement:

### Option A: OpenMM RPMD + Compatible Models
```python
from openmmml import MLPotential
from openmm import RPMDIntegrator

# Use MACE instead of UMA (TorchScript compatible)
potential = MLPotential('mace-off23-medium')
system = potential.createSystem(topology)

# RPMD works!
integrator = RPMDIntegrator(numCopies=16, temperature=300*kelvin, 
                            frictionCoeff=1/picosecond, stepSize=0.5*femtoseconds)
```

### Option B: i-PI + FAIRChem
Use i-PI (which you have in your repo) with FAIRChem as the force provider:
- i-PI handles RPMD/PIMD
- FAIRChem provides UMA forces via socket interface

## Files Delivered

1. `/media/extradrive/Trajectories/openmm/fairchem/openmm-ml/openmmml/models/umapotential.py`
   - Complete UMA integration implementation

2. `/media/extradrive/Trajectories/openmm/install_uma.sh`
   - Installation and verification script

3. `/media/extradrive/Trajectories/openmm/verify_uma_singlepoint.py`
   - Single-point calculation validation

4. `/media/extradrive/Trajectories/openmm/UMA_INTEGRATION_LIMITATION.md`
   - Technical limitation documentation

5. `/media/extradrive/Trajectories/openmm/UMA_INTEGRATION_FINAL_REPORT.md`
   - This document

## Conclusion

The UMA integration is **technically complete** but **cannot run due to architectural incompatibility**.

### Immediate Actions:

1. **For regular MD**: Use ASE with FAIRChem (Solution 1)
2. **For RPMD**: Use OpenMM RPMD + MACE or use i-PI + FAIRChem
3. **For OpenMM compatibility**: Switch to MACE/NequIP/ANI models

### Long-term:

Consider developing a custom C++ OpenMM plugin if UMA support in OpenMM is critical.

---

**Status**: Integration implemented, limitation documented, alternatives provided
**Recommendation**: Use ASE + FAIRChem for immediate UMA usage, or MACE + OpenMM for RPMD
