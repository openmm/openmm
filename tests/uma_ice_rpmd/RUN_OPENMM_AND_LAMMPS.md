# Run OpenMM and LAMMPS ice + UMA (same setup)

## What was fixed on your machine

1. **OpenMM UMA batch dtype** — FAIRChem requires `pos`, `cell`, and `cell_offsets` to share one dtype. The installed `openmmml` batch path was updated (copy from repo):
   ```bash
   cp wrappers/python/openmmml/models/umapotential_pythonforce_batch.py \
      "$(python -c "import openmmml,os; print(os.path.join(os.path.dirname(openmmml.__file__),'models')))")/"
   ```

2. **GPU OOM** — If the GPU is full (other processes), use **CPU** for OpenMM and UMA:
   ```bash
   --platform cpu --ml-device cpu
   ```
   Or free GPU memory, then use CUDA (much faster).

3. **LAMMPS Python** — `pip install lammps` provides `from lammps import lammps`. The wheel links **MPI**:
   ```bash
   conda install -y -c conda-forge openmpi
   export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
   ```
   Then `python -c "from lammps import lammps; lammps()"` should start without `libmpi.so.12` errors.

---

## OpenMM (trajectory: PDB)

From `tests/uma_ice_rpmd/`:

```bash
# GPU (default) — needs free VRAM
python test_uma_ice_rpmd.py --input ice.cif --molecules 16 --beads 1 \
  --model uma-s-1p1-pythonforce-batch --temperature 243 --dt 1.0 \
  --equil 1 --prod 10 --pressure 0 --output run_openmm --minimal

# CPU — slow but works when GPU is busy
python test_uma_ice_rpmd.py --input ice.cif --molecules 8 --beads 1 \
  --model uma-s-1p1-pythonforce-batch --temperature 243 --dt 1.0 \
  --equil 0.2 --prod 1 --pressure 0 --output run_openmm --minimal \
  --platform cpu --ml-device cpu
```

**Output:** `run_openmm/ice_rpmd_uma_T243_b1.pdb` (and NPZ if not `--minimal`).

---

## LAMMPS (trajectory: `.lammpstrj`)

```bash
cd lammps
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"   # after openmpi

# Quick test (8 molecules, ~700 steps)
python run_lammps_uma_ice.py --molecules 8 --device cpu \
  --infile in.ice_uma_quick.lmp
# → dump.ice_uma_quick.lammpstrj (quick input) or edit dump name in infile

# Longer / full box
python build_lammps_ice_data.py -i ../ice.cif -o data.ice_uma
python run_lammps_uma_ice.py --device cuda   # if GPU free
```

**Output:** `lammps/dump.ice_uma*.lammpstrj` — open in Ovito.

---

## One-liner checklist

| Step | Command |
|------|--------|
| openmmml fix | `cp .../umapotential_pythonforce_batch.py` → site-packages `models/` |
| LAMMPS + MPI | `pip install lammps` + `conda install openmpi` + `LD_LIBRARY_PATH` |
| FAIRChem | `pip install fairchem-core fairchem-lammps` |

If LAMMPS still fails on MPI, use only the **OpenMM** line above; same physics (UMA, same `ice.cif`).
