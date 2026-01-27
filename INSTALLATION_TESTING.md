# UMA Integration - Installation & Testing Guide

## Installation

### Step 1: Install OpenMM and OpenMM-Torch

```bash
# Install OpenMM
conda install -c conda-forge openmm

# Install OpenMM-Torch (required for GPU support)
conda install -c conda-forge openmmtorch
```

### Step 2: Install FAIRChem

```bash
pip install fairchem-core
```

### Step 3: Install OpenMM-ML with UMA Support

```bash
# Navigate to the openmm-ml directory
cd fairchem/openmm-ml

# Install in editable mode
pip install -e .
```

### Step 4: Verify Installation

```python
# Test basic import
from openmmml import MLPotential

# Try loading a UMA model (will download on first use)
potential = MLPotential('uma-s-1p1')
print("UMA integration successful!")
```

## Testing

### Run Unit Tests

```bash
# Navigate to test directory
cd fairchem/openmm-ml/test

# Run UMA tests with pytest
pytest TestUMAPotential.py -v

# Run specific test
pytest TestUMAPotential.py::TestUMA::testCreatePureMLSystem -v
```

### Run RPMD Integration Test

```bash
# From repository root
python test_rpmd_uma.py
```

Expected output:
- Basic RPMD test passes
- PILE_G thermostat test passes
- All energies are finite
- No errors during simulation

### Run Validation

```bash
cd fairchem/openmm-ml/examples/uma
python validate_uma.py
```

This compares OpenMM-ML energies with FAIRChem ASE calculator reference.

### Run GPU Benchmark

```bash
cd fairchem/openmm-ml/examples/uma
python benchmark_gpu.py
```

This will test performance on available platforms (CPU, CUDA, OpenCL).

## Quick Examples

### Example 1: Basic MD Simulation

```bash
cd fairchem/openmm-ml/examples/uma
python run_uma_basic.py
```

### Example 2: RPMD Simulation

```bash
cd fairchem/openmm-ml/examples/uma
python run_uma_rpmd.py
```

## Troubleshooting

### Issue: Model download fails

**Solution:**
```python
# Check cache directory
from fairchem.core._config import CACHE_DIR
print(CACHE_DIR)

# Ensure you have internet connection
# Models are ~100MB each and download from HuggingFace
```

### Issue: GPU not detected

**Solution:**
```bash
# Check OpenMM-Torch installation
conda list | grep openmmtorch

# Check PyTorch CUDA
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"

# Reinstall OpenMM-Torch if needed
conda install -c conda-forge openmmtorch --force-reinstall
```

### Issue: ImportError for FAIRChem

**Solution:**
```bash
# Install FAIRChem with all dependencies
pip install fairchem-core

# Or install from source
git clone https://github.com/FAIR-Chem/fairchem.git
cd fairchem
pip install -e packages/fairchem-core
```

### Issue: Memory error with RPMD

**Solution:**
- Reduce number of beads (use 4-8 instead of 32)
- Use smaller system
- Use CPU platform if GPU memory is insufficient
- Enable gradient checkpointing (future feature)

### Issue: Test failures

**Solution:**
```bash
# Ensure all dependencies are up to date
pip install --upgrade openmm fairchem-core torch

# Check if test data exists
ls fairchem/openmm-ml/test/data/

# Run tests with verbose output
pytest TestUMAPotential.py -v -s
```

## Verifying GPU Acceleration

```python
import openmm
from openmmml import MLPotential

# Create system
potential = MLPotential('uma-s-1p1')
system = potential.createSystem(topology, task_name='omol')

# Try CUDA platform
try:
    platform = openmm.Platform.getPlatformByName('CUDA')
    context = openmm.Context(system, integrator, platform)
    print("✓ GPU acceleration working!")
except:
    print("✗ CUDA not available, using CPU")
```

## Performance Expectations

Typical performance on a modern GPU (NVIDIA RTX 3090):
- Small molecule (10-20 atoms): 1000-5000 steps/s
- Medium molecule (50-100 atoms): 500-1000 steps/s
- RPMD (8 beads, 20 atoms): 100-500 steps/s

CPU performance is typically 10-50x slower.

## Known Limitations

1. **First run slow**: Model download and compilation happens on first use
2. **Memory usage**: RPMD with many beads requires significant GPU memory
3. **Turbo mode**: Only works with fixed atomic composition
4. **Task specification**: Must match model capabilities

## Getting Help

- **OpenMM-ML Issues**: https://github.com/openmm/openmm-ml/issues
- **FAIRChem Documentation**: https://fair-chem.github.io/
- **OpenMM Forum**: https://github.com/openmm/openmm/discussions

## Success Checklist

- [ ] OpenMM installed (version >= 8.4)
- [ ] OpenMM-Torch installed
- [ ] FAIRChem-core installed
- [ ] OpenMM-ML installed in editable mode
- [ ] Can import MLPotential
- [ ] Can load UMA model
- [ ] Unit tests pass
- [ ] RPMD test passes
- [ ] GPU detected (if available)
- [ ] Example scripts run successfully

If all items are checked, your UMA integration is ready to use!
