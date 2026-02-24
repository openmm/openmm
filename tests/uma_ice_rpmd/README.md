# UMA Ice RPMD Simulation

Test if ice remains frozen at 243K using UMA potential with RPMD.

## Overview

This simulation tests the stability of ice at sub-freezing temperatures using:
- **UMA potential**: Smallest model (uma-s-1-pythonforce-batch) for fast runs
- **RPMD**: Ring Polymer Molecular Dynamics with 8 beads
- **Temperature**: 243 K (30 degrees below freezing)
- **Ice Ih structure**: Proper tetrahedral hydrogen-bonding network (GenIce2 or embedded CIF fallback)

## Ice Structure

The initial structure uses proper ice Ih:
- **GenIce2** (optional): Proton-disordered ice via `pip install genice2`
- **Embedded CIF fallback**: Ice Ih from Avogadro/COD (12 molecules) replicated to target size

## Key Physics

- **Ice Ih**: Hexagonal ice, O-O ~2.75 Å, density ~0.92 g/cm³
- **RPMD**: Captures quantum nuclear effects for hydrogen bonding
- **Timestep**: 0.5-1.0 fs recommended for flexible water + RPMD

## Expected Results

### If Ice Remains Frozen
- **NPT**: Density > 0.85 g/cm³ (ice ~0.92)
- **NVT**: RMSD < 0.15 nm from minimized structure
- Sharp O-O RDF first peak at ~0.275 nm

### If Ice Melts
- **NPT**: Density ~1.0 g/cm³ (liquid)
- **NVT**: RMSD > 0.15 nm
- Broader RDF

## Usage

### Quick Test
```bash
python test_uma_ice_rpmd.py --molecules 32 --beads 8 --temperature 243 \
    --dt 1.0 --equil 5 --prod 50 --pressure 1 --model uma-s-1-pythonforce-batch
```

### NVT (no barostat)
```bash
python test_uma_ice_rpmd.py --molecules 32 --beads 8 --temperature 243 \
    --equil 10 --prod 100 --pressure 0
```

### Custom Parameters
```bash
python test_uma_ice_rpmd.py \
    --molecules 32 \           # Number of water molecules
    --beads 8 \                # RPMD beads
    --temperature 243 \        # Temperature (K)
    --dt 1.0 \                 # Timestep in fs (0.5-1.0 for stability)
    --equil 5.0 \              # Equilibration (ps)
    --prod 50.0 \              # Production (ps)
    --model uma-s-1-pythonforce-batch \   # Smallest UMA model
    --output ice_test          # Output directory
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

### CUDA Error (CUDA_ERROR_ILLEGAL_ADDRESS)
When OpenMM and PyTorch both use CUDA on the same GPU, a context conflict can cause
`CUDA_ERROR_ILLEGAL_ADDRESS`. To avoid this:
- **ML model uses CPU by default** when OpenMM uses CUDA (OpenMM on GPU, UMA on CPU).
- Use `--ml-device cuda` to try GPU for ML (may fail on some systems).
- Use `--platform cpu` to run everything on CPU if GPU issues persist.

## References

- UMA: Universal Molecular Atomistic potentials (FAIRChem)
- RPMD: Ring Polymer Molecular Dynamics
- Ice Ih: Hexagonal ice structure
