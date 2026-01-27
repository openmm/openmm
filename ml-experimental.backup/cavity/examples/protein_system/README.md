# Protein Cavity Test (3UTL)

This test mirrors the water-system workflow but uses a protein in explicit
solvent and couples the cavity mode to the OH stretch frequency.

## Dependencies
- `openmm`
- `pdbfixer`
- `openmmforcefields` (for AMBER ff14SB + TIP4P-Ew XMLs)

## Running
Default (downloads 3UTL, cleans with PDBFixer, solvates, runs 100+900 ps):

```
python /media/extradrive/Trajectories/openmm/tests/protein_system/run_simulation.py
```

Quick test mode (short run):

```
python /media/extradrive/Trajectories/openmm/tests/protein_system/run_simulation.py --test
```

Use a local PDB instead of downloading:

```
python /media/extradrive/Trajectories/openmm/tests/protein_system/run_simulation.py --pdb-path /path/to/3UTL.pdb
```

## Output
The script writes `protein_cavity_lambdaXXXX.npz` with the same schema as the
water-case:

- `time_ps`: time array (ps) for the production window
- `dipole_nm`: dipole trajectory in e·nm
- `cavity_nm`: cavity coordinate trajectory (nm)
- `metadata`: run parameters (lambda, cavity freq, force field, padding, etc.)
