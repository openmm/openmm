# OMol25 Dataset + LES Integration Guide for UMA

This guide explains how to access the OMol25 dataset and supplement UMA models with LES (Latent Ewald Summation) for enhanced long-range electrostatic interactions.

## 📊 Overview

**OMol25 Dataset:**
- 100M+ DFT single-point calculations and relaxations
- Organic/inorganic molecules, transition metal complexes, electrolytes
- Level of theory: wB97M-V/def2-TZVPD (ORCA6) with non-local dispersion
- Format: ASE DB compatible LMDB files (*.aselmdb)
- License: CC BY 4.0

**LES (Latent Ewald Summation):**
- Universal augmentation framework for long-range electrostatics
- Already integrated with MACE, NequIP, Allegro, CACE, CHGNet
- Enables Born Effective Charge (BEC) inference
- Handles electric field responses and polarization

**Goal:** Combine UMA's molecular accuracy with LES's long-range electrostatics for improved training.

---

## 🔑 Step 1: Accessing OMol25 Dataset

### HuggingFace Access

The dataset requires accepting a gated access agreement.

```bash
# 1. Go to HuggingFace and request access
# https://huggingface.co/facebook/OMol25

# 2. After approval, install huggingface-cli
pip install huggingface_hub[cli]

# 3. Login with your token
huggingface-cli login

# 4. Download the dataset splits
# Option A: Download specific splits (recommended to start)
huggingface-cli download facebook/OMol25 \
  --repo-type dataset \
  --include "train_4M/*" \
  --local-dir ./omol25_data

# Option B: Full dataset (very large - ~petabytes)
# Start with train_4M (4M subset) first
```

### Dataset Structure

```
omol25_data/
├── train_4M/          # 4M randomly sampled subset
│   └── *.aselmdb      # ASE-compatible LMDB files
├── train_all/         # Full training set (very large)
├── val/               # Validation set
└── test/              # Test set
```

### Loading OMol25 in Python

```python
from fairchem.core.datasets import AseDBDataset

# Load dataset
dataset = AseDBDataset({
    "src": "omol25_data/train_4M/"
})

# Access individual structures
for idx in range(len(dataset)):
    atoms = dataset.get_atoms(idx)
    
    # Check metadata
    print(f"Charge: {atoms.info.get('charge', 0)}")
    print(f"Spin: {atoms.info.get('spin', 1)}")
    print(f"Energy: {atoms.get_potential_energy()} eV")
    print(f"Forces: {atoms.get_forces()} eV/Å")
```

### Dataset Metadata

Each structure contains:
- **Atomic positions** (Å)
- **Energies** (eV)
- **Forces** (eV/Å)
- **Charge** (stored in `atoms.info['charge']`)
- **Spin multiplicity** (stored in `atoms.info['spin']`)
- **Source path** (for tracking)

---

## 🧬 Step 2: Understanding LES Integration

### What LES Provides

From the `LES-BEC` directory, LES adds:

1. **Long-range electrostatics** via Ewald summation
2. **Born Effective Charges (BEC)** inference
3. **Electric field response**
4. **Polarization currents**

### LES Implementation (CACE Example)

The LES method is already implemented in CACE with these modules:

```python
# From cace/modules/
from cace.modules import EwaldPotential      # Long-range electrostatics
from cace.modules import Polarization        # Polarization effects
from cace.modules import AngularComponent    # Angular decomposition
```

### Example: CACE-LR with LES

```python
# Training example from cace/water/
from cace.modules import (
    EwaldPotential,
    Polarization,
    BesselRBF,
    PolynomialCutoff
)

# Ewald potential module
ewald = EwaldPotential(
    dl=2,                    # Decomposition level
    sigma=1.0,               # Ewald parameter
    exponent=1,              # Power for 1/r potential
    feature_key='q',         # Charge feature
    aggregation_mode='sum',
    compute_field=True       # Enable electric field
)

# Polarization module
polarization = Polarization(
    # Configuration for induced dipoles
)
```

---

## 🔄 Step 3: Integrating LES with UMA

### Current UMA Architecture

UMA uses:
- **Mixture-of-Linear-Experts (MoLE)** routing
- **Task-specific embeddings** (omol, omat, oc20, odac, omc)
- **Equivariant graph neural network**
- **Single output head** for all tasks

### Proposed Integration Strategy

There are **three approaches** to supplement UMA with LES:

#### Option A: Pre-training Data Augmentation (Recommended)

Augment OMol25 training data with LES-computed long-range terms:

```python
from fairchem.core.scripts.create_uma_finetune_dataset import create_dataset
import torch

def augment_omol25_with_les(omol25_path, output_path):
    """
    Add LES long-range electrostatic features to OMol25 data
    """
    from cace.modules import EwaldPotential
    
    # Load OMol25
    dataset = AseDBDataset({"src": omol25_path})
    
    # Initialize Ewald module
    ewald = EwaldPotential(dl=2, sigma=1.0, exponent=1, 
                          feature_key='q', compute_field=True)
    
    augmented_data = []
    
    for idx in range(len(dataset)):
        atoms = dataset.get_atoms(idx)
        
        # Compute LES correction
        positions = torch.tensor(atoms.get_positions())
        charges = estimate_charges(atoms)  # From UMA charge prediction
        cell = torch.tensor(atoms.get_cell())
        
        # Get Ewald energy + field
        ewald_energy, ewald_field = ewald.compute_potential_optimized(
            positions, charges, cell, compute_field=True
        )
        
        # Add to atoms.info for training
        atoms.info['ewald_energy'] = ewald_energy.item()
        atoms.arrays['ewald_field'] = ewald_field.numpy()
        
        augmented_data.append(atoms)
    
    # Write augmented dataset
    write_lmdb(augmented_data, output_path)
```

#### Option B: Fine-tuning UMA on LES-augmented Dataset

Fine-tune existing UMA checkpoint on LES-enhanced data:

```bash
# 1. Create LES-augmented dataset
python augment_omol25_with_les.py \
  --omol25-path ./omol25_data/train_4M/ \
  --output-path ./omol25_les_augmented/

# 2. Create fine-tuning dataset
python src/fairchem/core/scripts/create_uma_finetune_dataset.py \
  --train-dir ./omol25_les_augmented/train \
  --val-dir ./omol25_les_augmented/val \
  --output-dir ./omol25_les_finetune \
  --uma-task omol \
  --regression-task ef

# 3. Fine-tune UMA
fairchem -c ./omol25_les_finetune/uma_sm_finetune_template.yaml \
  base_model_name=uma-s-1p1 \
  epochs=10 \
  batch_size=16 \
  lr=4e-4
```

#### Option C: Hybrid Architecture (Advanced)

Add LES as an auxiliary loss during UMA training:

```python
# Custom training loop modification
class UMAWithLESLoss:
    def __init__(self, uma_model, ewald_module):
        self.uma_model = uma_model
        self.ewald = ewald_module
    
    def compute_loss(self, batch):
        # Standard UMA loss
        uma_pred = self.uma_model(batch)
        energy_loss = F.mse_loss(uma_pred['energy'], batch['energy'])
        force_loss = F.mse_loss(uma_pred['forces'], batch['forces'])
        
        # LES electrostatic correction
        charges = uma_pred.get('charges', estimate_charges(batch))
        ewald_energy, _ = self.ewald.compute_potential_optimized(
            batch['positions'], charges, batch['cell']
        )
        
        # Auxiliary LES loss
        les_loss = F.mse_loss(
            uma_pred['energy'] + ewald_energy,
            batch['energy']
        )
        
        # Combined loss
        total_loss = energy_loss + force_loss + 0.1 * les_loss
        return total_loss
```

---

## 📝 Step 4: Training Configuration

### ASE LMDB Dataset Configuration

```yaml
# omol25_les_config.yaml
datasets:
  omol_train:
    format: ase_db
    splits:
      train:
        src: ./omol25_les_augmented/train/
        a2g_args:
          r_energy: true
          r_forces: true
          r_stress: false
      val:
        src: ./omol25_les_augmented/val/
        a2g_args:
          r_energy: true
          r_forces: true
```

### UMA Fine-tuning Parameters

```yaml
base_model_name: uma-s-1p1
max_neighbors: 300

# Training hyperparameters
epochs: 10
batch_size: 16
lr: 4e-4
weight_decay: 1e-3

# Task configuration
task_name: omol
charge: 0
spin: 1

# LES-specific augmentation
augmentations:
  ewald:
    enabled: true
    sigma: 1.0
    decomposition_level: 2
```

---

## 🔬 Step 5: Validation & Testing

### Test LES Integration

```python
import openmm
from openmm import app, unit
from openmmml import MLPotential
import numpy as np

# Load fine-tuned UMA+LES model
potential = MLPotential('path/to/finetuned_uma_les.pt')

# Create test system (water molecule in cavity)
topology = create_water_topology()
positions = load_positions()

system = potential.createSystem(
    topology, 
    task_name='omol',
    charge=0,
    spin=1
)

# Run MD with electric field
platform = openmm.Platform.getPlatformByName('CUDA')
integrator = openmm.LangevinMiddleIntegrator(
    300*unit.kelvin, 1/unit.picosecond, 1*unit.femtoseconds
)
context = openmm.Context(system, integrator, platform)
context.setPositions(positions)

# Apply external electric field
# (This requires custom force implementation)
electric_field = openmm.CustomExternalForce("q*Ex*x + q*Ey*y + q*Ez*z")
electric_field.addPerParticleParameter("q")
electric_field.addGlobalParameter("Ex", 0.0)
electric_field.addGlobalParameter("Ey", 0.0)
electric_field.addGlobalParameter("Ez", 0.01)  # V/Å

# Run simulation
integrator.step(10000)
```

### Benchmark Against Pure UMA

```python
from fairchem.core import FAIRChemCalculator
from ase import Atoms

atoms = Atoms('H2O', positions=[[0,0,0], [0.96,0,0], [-0.24,0.93,0]])

# UMA baseline
calc_uma = FAIRChemCalculator('uma-s-1p1', task_name='omol')
atoms.calc = calc_uma
energy_uma = atoms.get_potential_energy()
forces_uma = atoms.get_forces()

# UMA + LES
calc_les = FAIRChemCalculator('uma_les_finetuned.pt', task_name='omol')
atoms.calc = calc_les
energy_les = atoms.get_potential_energy()
forces_les = atoms.get_forces()

print(f"Energy difference: {energy_les - energy_uma:.4f} eV")
```

---

## 📦 Step 6: Data Pipeline Summary

```
┌─────────────────┐
│  OMol25 Dataset │  (HuggingFace: facebook/OMol25)
│   100M+ DFT     │
│   wB97M-V       │
└────────┬────────┘
         │
         │ Download & Extract
         │
         ▼
┌─────────────────┐
│  ASE LMDB Files │  train_4M/*.aselmdb
│  - positions    │
│  - energies     │
│  - forces       │
│  - charge/spin  │
└────────┬────────┘
         │
         │ Augment with LES
         │
         ▼
┌─────────────────┐
│   LES Module    │  (from cace.modules)
│  - EwaldPot     │  Compute long-range terms
│  - Polarization │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Augmented LMDB  │  + ewald_energy
│                 │  + ewald_field
│                 │  + polarization
└────────┬────────┘
         │
         │ Fine-tune
         │
         ▼
┌─────────────────┐
│  UMA + LES      │  Enhanced model
│  Fine-tuned     │  - Better long-range
│  Checkpoint     │  - BEC inference
└────────┬────────┘
         │
         │ Deploy
         │
         ▼
┌─────────────────┐
│  OpenMM-ML      │  Production simulations
│  Simulation     │  with electric field
└─────────────────┘
```

---

## 🛠️ Implementation Checklist

### Phase 1: Data Access (Week 1)
- [ ] Request OMol25 access on HuggingFace
- [ ] Download train_4M subset (~hundreds of GB)
- [ ] Verify LMDB files load correctly
- [ ] Explore dataset statistics (charge/spin distributions)

### Phase 2: LES Integration (Week 2-3)
- [ ] Install CACE with LES modules
- [ ] Test EwaldPotential on sample molecules
- [ ] Develop charge estimation from UMA
- [ ] Implement data augmentation script
- [ ] Augment 10k samples as pilot

### Phase 3: Fine-tuning (Week 3-4)
- [ ] Create fine-tuning dataset with augmented data
- [ ] Configure UMA fine-tuning YAML
- [ ] Run fine-tuning on pilot dataset
- [ ] Monitor convergence and validation metrics
- [ ] Scale to full train_4M

### Phase 4: Validation (Week 4-5)
- [ ] Test on molecules with known BEC
- [ ] Compare to pure UMA baseline
- [ ] Validate electric field response
- [ ] Benchmark OpenMM-ML integration
- [ ] Run RPMD simulations

### Phase 5: Production (Week 6+)
- [ ] Deploy to cavity dimer system
- [ ] Compute IR spectra with electric field
- [ ] Compare to current CACE-LR results
- [ ] Document performance improvements
- [ ] Publish enhanced model checkpoint

---

## 📚 Key References

### OMol25
- Paper: [The Open Molecules 2025 (OMol25) Dataset](https://arxiv.org/abs/2505.08762)
- HuggingFace: https://huggingface.co/facebook/OMol25
- Documentation: https://fair-chem.github.io/molecules/datasets/omol25.html

### UMA
- Paper: [UMA: A Family of Universal Models for Atoms](https://arxiv.org/abs/2506.23971)
- Models: https://huggingface.co/facebook/UMA
- Documentation: https://fair-chem.github.io/core/uma.html

### LES
- Paper: [Latent Ewald summation for machine learning](https://www.nature.com/articles/s41524-025-01577-7)
- Universal Framework: [arXiv:2507.14302](https://arxiv.org/abs/2507.14302)
- BEC Paper: [Machine learning interatomic potential can infer electrical response](https://arxiv.org/abs/2504.05169)
- Code: https://github.com/ChengUCB/les

### FAIRChem Training
- Fine-tuning Guide: https://fair-chem.github.io/core/common_tasks/fine_tuning.html
- Training Docs: https://fair-chem.github.io/core/common_tasks/training.html

---

## ⚡ Quick Start Commands

```bash
# 1. Access OMol25
huggingface-cli login
huggingface-cli download facebook/OMol25 --include "train_4M/*" --local-dir ./omol25_data

# 2. Install dependencies
pip install fairchem-core
cd cace && pip install -e .

# 3. Augment dataset (create this script)
python scripts/augment_omol25_with_les.py \
  --input ./omol25_data/train_4M \
  --output ./omol25_les

# 4. Create fine-tuning config
python src/fairchem/core/scripts/create_uma_finetune_dataset.py \
  --train-dir ./omol25_les/train \
  --val-dir ./omol25_les/val \
  --output-dir ./finetune_config \
  --uma-task omol --regression-task ef

# 5. Fine-tune UMA
fairchem -c ./finetune_config/uma_sm_finetune_template.yaml \
  base_model_name=uma-s-1p1 epochs=10 batch_size=16

# 6. Test in OpenMM
python test_uma_les_openmm.py --checkpoint ./finetune_runs/final/inference_ckpt.pt
```

---

## 🎯 Expected Outcomes

1. **Enhanced Long-Range Accuracy**: Better electrostatics beyond UMA's cutoff
2. **BEC Prediction**: Ability to infer Born Effective Charges
3. **Electric Field Response**: Accurate polarization under external fields
4. **IR Spectra Quality**: Improved cavity coupling predictions
5. **Performance**: Minimal overhead compared to pure UMA (<5% slower)

---

## 💡 Next Steps

After this integration, you can:
1. Apply to cavity-molecule systems (your dimer tests)
2. Compute field-dependent IR spectra
3. Study polariton formation with ML potential
4. Scale RPMD to larger systems with UMA+LES
5. Publish methodology and benchmarks

---

## 🤝 Support & Contact

- FAIRChem Issues: https://github.com/facebookresearch/fairchem/issues
- LES/CACE: https://github.com/BingqingCheng/cace/issues
- OpenMM-ML: https://github.com/openmm/openmm-ml/issues

For questions specific to this integration:
- Check existing issues in the repos above
- Consider posting on OpenMM forums
- Review the paper discussions

---

## 🎯 Implementation Status

**IMPLEMENTED!** All configuration files and scripts have been created.

### Created Files

#### Configuration Files (in `fairchem/configs/uma/lr/`)
- ✅ `dataset/omol25.yaml` - OMol25 dataset configuration
- ✅ `tasks/omol25.yaml` - Energy and force tasks
- ✅ `element_refs/omol25_element_reference.yaml` - Element reference energies
- ✅ `backbone/uma_sm_les.yaml` - UMA+LES backbone configuration
- ✅ `uma_sm_omol25_les_train.yaml` - Main training configuration
- ✅ `README.md` - Complete training documentation

#### Scripts (in `scripts/`)
- ✅ `download_omol25.sh` - Download OMol25 from HuggingFace
- ✅ `setup_uma_les_training.sh` - Complete setup and dependency checker
- ✅ `validate_uma_les.py` - Model validation on test set
- ✅ `test_uma_les_openmm.py` - OpenMM-ML integration test
- ✅ `augment_omol25_with_les.py` - Data augmentation (optional)

### Quick Start

```bash
# 1. Run setup script
cd /media/extradrive/Trajectories/openmm
./scripts/setup_uma_les_training.sh

# 2. Download data (if not already done)
./scripts/download_omol25.sh

# 3. Test configuration
cd fairchem
fairchem -c configs/uma/lr/uma_sm_omol25_les_train.yaml \
  job.debug=True epochs=1 batch_size=2

# 4. Run training
fairchem -c configs/uma/lr/uma_sm_omol25_les_train.yaml
```

**Status**: Ready for training
**Last Updated**: 2026-01-23
**Version**: 2.0 - Implementation Complete
