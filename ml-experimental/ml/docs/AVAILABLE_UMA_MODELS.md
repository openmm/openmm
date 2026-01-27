# Available UMA and eSEN Models with Batched RPMD Support

All models listed below now support efficient batched RPMD inference via the `-pythonforce-batch` suffix.

## UMA Models (Universal Molecular Atomistic)

### 1. `uma-s-1-pythonforce-batch`
- **Size**: Small (~2.4M parameters)
- **Version**: 1.0
- **Speed**: Fastest
- **Accuracy**: Good
- **Use case**: Quick testing, large systems

### 2. `uma-s-1p1-pythonforce-batch` ⭐ (Current Default)
- **Size**: Small (~2.4M parameters)
- **Version**: 1.1 (improved)
- **Speed**: Fast
- **Accuracy**: Better than v1.0
- **Use case**: Production runs, recommended for most use cases

### 3. `uma-m-1p1-pythonforce-batch`
- **Size**: Medium (larger model)
- **Version**: 1.1
- **Speed**: Moderate (slower than small)
- **Accuracy**: Best
- **Use case**: High-accuracy simulations, smaller systems

## OMol25 eSEN Models (Organic Molecules)

Trained on OMol25 dataset - specialized for organic molecular systems.

### 4. `esen-sm-direct-all-omol-pythonforce-batch`
- **Size**: Small
- **Training**: Direct force learning
- **Use case**: Organic molecules, direct force prediction

### 5. `esen-sm-conserving-all-omol-pythonforce-batch`
- **Size**: Small
- **Training**: Energy-conserving
- **Use case**: Organic molecules, energy conservation critical

### 6. `esen-md-direct-all-omol-pythonforce-batch`
- **Size**: Medium
- **Training**: Direct force learning
- **Use case**: Organic molecules, higher accuracy needed

## OC25 eSEN Models (Catalysis)

Trained on OC25 dataset - specialized for catalytic systems and surfaces.

### 7. `esen-sm-conserving-all-oc25-pythonforce-batch`
- **Size**: Small
- **Training**: Energy-conserving
- **Use case**: Catalytic systems, surfaces

### 8. `esen-md-direct-all-oc25-pythonforce-batch`
- **Size**: Medium
- **Training**: Direct force learning
- **Use case**: Catalytic systems, higher accuracy

---

## Usage Example

```python
from openmmml import MLPotential

# Use any of the models above
potential = MLPotential('uma-m-1p1-pythonforce-batch')

# For organic molecules, try OMol25 models
potential = MLPotential('esen-sm-conserving-all-omol-pythonforce-batch')

# For catalysis/surfaces, try OC25 models
potential = MLPotential('esen-sm-conserving-all-oc25-pythonforce-batch')

# Create system as usual
system = potential.createSystem(
    topology,
    task_name='omol',  # or appropriate task name
    charge=0,
    spin=1
)
```

## Performance Notes

- All `-pythonforce-batch` models use the **same batched RPMD implementation**
- Batched inference evaluates all RPMD beads in a **single GPU call**
- Expected speedup: ~2-3x compared to sequential evaluation
- Memory usage scales with model size × number of beads
- **All models verified working** with both small and large systems (up to 150 atoms, 16 beads)

## Model Selection Guide

| System Type | Recommended Model | Reason |
|-------------|-------------------|--------|
| Water/ice systems | `uma-s-1p1-pythonforce-batch` | General purpose, well-tested |
| Small organic molecules | `esen-sm-conserving-all-omol-pythonforce-batch` | Trained on OMol25 |
| Large organic systems | `esen-sm-direct-all-omol-pythonforce-batch` | Faster, good accuracy |
| Catalytic surfaces | `esen-sm-conserving-all-oc25-pythonforce-batch` | Specialized for catalysis |
| High-accuracy needed | `uma-m-1p1-pythonforce-batch` | Largest, most accurate |

---

**Last Updated**: 2026-01-25
