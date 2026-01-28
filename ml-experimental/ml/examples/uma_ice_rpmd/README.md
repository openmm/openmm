# UMA Ice RPMD Simulation

Test if ice remains frozen at 243K using UMA potential with RPMD.

## Overview

This simulation tests the stability of ice at sub-freezing temperatures using:
- **UMA potential**: Universal Molecular Atomistic ML potential (uma-s-1p1-pythonforce)
- **RPMD**: Ring Polymer Molecular Dynamics with 8 beads
- **Temperature**: 243 K (30 degrees below freezing)
- **System**: 32 water molecules in ice Ih structure

## Key Physics

- **Ice Ih structure**: Hexagonal ice, the common form of ice
- **RPMD**: Captures quantum nuclear effects important for hydrogen bonding in water
- **Temperature**: 243 K is well below the melting point (273 K), so ice should remain frozen

## Expected Results

### If Ice Remains Frozen
- RMSD < 0.1 nm (atoms stay close to initial positions)
- Density ~0.92 g/cm³ (ice density)
- Sharp peaks in O-O radial distribution function
- Low mean-squared displacement (MSD)

### If Ice Melts
- RMSD > 0.1 nm (significant atomic rearrangement)
- Density ~1.0 g/cm³ (liquid water density)
- Broader peaks in RDF
- Higher MSD

## Usage

### Quick Test (10 ps equilibration + 100 ps production)
```bash
python test_uma_ice_rpmd.py --molecules 32 --beads 8 --temperature 243 \
                             --equil 10 --prod 100
```

### Longer Production Run
```bash
python test_uma_ice_rpmd.py --molecules 32 --beads 8 --temperature 243 \
                             --equil 50 --prod 500
```

### Custom Parameters
```bash
python test_uma_ice_rpmd.py \
    --molecules 32 \        # Number of water molecules
    --beads 8 \             # RPMD beads
    --temperature 243 \     # Temperature (K)
    --dt 1.0 \              # Timestep (fs)
    --equil 10.0 \          # Equilibration (ps)
    --prod 100.0 \          # Production (ps)
    --output ice_test       # Output prefix
```

## Output Files

### Data File (NPZ)
`uma_ice_rpmd_T243K_b8.npz` contains:
- `times`: Time points (ps)
- `rmsds`: RMSD from initial structure (nm)
- `densities`: System density (g/cm³)
- `potential_energies`: Mean PE per bead (kJ/mol)
- `kinetic_energies`: Mean KE per bead (kJ/mol)
- `temperatures`: Temperature per bead (K)
- `msds`: Mean-squared displacement (nm²)
- `r_initial`, `g_initial`: Initial RDF
- `r_final`, `g_final`: Final RDF
- `initial_positions`, `final_positions`: Positions
- `is_frozen`: Boolean flag for ice status

### Analysis Plot
`uma_ice_rpmd_T243K_b8_analysis.png` shows:
1. RMSD vs time
2. Density vs time
3. Potential energy vs time
4. Temperature vs time
5. MSD vs time
6. RDF comparison (initial vs final)

## Analysis

### Loading Data
```python
import numpy as np
import matplotlib.pyplot as plt

data = np.load('uma_ice_rpmd_T243K_b8.npz')
times = data['times']
rmsds = data['rmsds']
is_frozen = data['is_frozen']

print(f"Ice frozen: {is_frozen}")
print(f"Final RMSD: {rmsds[-1]:.4f} nm")
```

### Indicators of Frozen Ice
- RMSD stays below 0.1 nm throughout simulation
- Density remains around 0.92 g/cm³
- Sharp first peak in RDF at ~0.275 nm (nearest neighbor O-O distance)
- Low MSD indicates restricted motion

## Requirements

- OpenMM with RPMD plugin
- openmm-ml with UMA support
- ASE (Atomic Simulation Environment)
- NumPy
- Matplotlib

## Notes

- Ice Ih is generated using ASE with appropriate lattice parameters
- RPMD with 8 beads captures quantum effects at 243 K
- 1 fs timestep is appropriate for flexible water models
- UMA potential provides accurate ML-based forces and energies
- Simulation uses PILE thermostat (default for RPMD)

## Performance

- Typical performance: 50-200 steps/s on GPU (CUDA)
- 10 ps + 100 ps simulation: ~5-20 minutes on GPU
- CPU is slower but functional for small systems

## Troubleshooting

### ModuleNotFoundError: No module named 'openmm'

OpenMM is not installed in the **current Python/conda environment**. You are likely in `(base)` or a env without OpenMM.

**Fix:** Install OpenMM and related packages, then run from that environment.

**Option A – Conda (recommended):**
```bash
# Create a dedicated env with OpenMM
conda create -n openmm-ml python=3.10 -y
conda activate openmm-ml
conda install -c conda-forge openmm openmmtorch
pip install openmmml  # or: pip install -e /path/to/openmm-ml

# Run from this env
cd /media/extradrive/Trajectories/openmm/ml-experimental/ml/examples/uma_ice_rpmd
python test_uma_ice_rpmd.py -h
```

**Option B – Install into current env:**
```bash
conda install -c conda-forge openmm openmmtorch
# Then ensure openmmml is available (see below).
```

**Check which env has OpenMM:**
```bash
conda activate your_openmm_env   # e.g. the one you use for OpenMM work
python -c "import openmm; print(openmm.__file__)"
```

Always run `python test_uma_ice_rpmd.py` from the **same** environment where OpenMM is installed.

### Import Error: openmmml
```bash
cd /media/extradrive/Trajectories/openmm/fairchem/openmm-ml
pip install -e .
```

### Import Error: ase
```bash
pip install ase
```

### CUDA Out of Memory
Reduce system size or use CPU platform:
- Reduce `--molecules` to 16 or 8
- Script will automatically fall back to CPU if CUDA fails

## References

- UMA: Universal Molecular Atomistic potentials (FAIRChem)
- RPMD: Ring Polymer Molecular Dynamics
- Ice Ih: Hexagonal ice structure
