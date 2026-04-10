# `libgfortran.so.3` without `apt install libgfortran3`

On **Ubuntu 24.04+** the Debian package **`libgfortran3` is gone**, so `sudo apt install libgfortran3` fails with *Unable to locate package*.

Your **omnia `mpich 3.2`** was built against **old libgfortran (.so.3)**. The pip **LAMMPS** wheel loads that MPI → missing `.so.3`.

## Fix A — Conda legacy Fortran (no sudo)

```bash
conda install -y -c conda-forge "libgfortran=3.0.0"
```

That package ships the old runtime. Then:

```bash
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$HOME/miniconda3/lib/python3.12/site-packages/lammps.libs:$LD_LIBRARY_PATH"
cd /media/extradrive/Trajectories/openmm/tests/uma_ice_rpmd/lammps
python -c "from lammps import lammps; lammps(cmdargs=['-log','none']); print('OK')"
python run_lammps_uma_ice.py --molecules 8 --device cuda --infile in.ice_uma_quick.lmp
```

## Fix B — Drop old MPICH, use conda-forge MPICH (often cleaner)

Removes the `.so.3` requirement if nothing else needs omnia mpich:

```bash
conda remove -y mpich --force
conda install -y -c conda-forge mpich
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$HOME/miniconda3/lib/python3.12/site-packages/lammps.libs:$LD_LIBRARY_PATH"
```

Then the same `python -c "from lammps import lammps"` test.

**Warning:** other tools in `base` that depended on omnia mpich could break; use a dedicated env if unsure.

## If LAMMPS still fails

Use **OpenMM** only (same UMA ice physics):

```bash
cd /media/extradrive/Trajectories/openmm/tests/uma_ice_rpmd
python test_uma_ice_rpmd.py --input ice.cif --molecules 8 --beads 1 \
  --model uma-s-1p1-pythonforce-batch --temperature 243 --dt 1.0 \
  --equil 0.1 --prod 0.5 --pressure 0 --output run_now --minimal \
  --platform cpu --ml-device cpu
```
