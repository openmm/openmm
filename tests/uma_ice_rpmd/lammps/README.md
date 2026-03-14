# LAMMPS + UMA ice (same system as OpenMM `test_uma_ice_rpmd.py`)

This folder runs **the same ice Ih cell and UMA potential** as the OpenMM RPMD test, but in **LAMMPS**
via FAIRChem’s **`fix external`** driver (`fairchem-lammps`).

**Ovito-friendly output**

| File | Use |
|------|-----|
| `dump.ice_uma_quick_ovito.xyz` | **XYZ** — first column is **O** or **H** (via `dump_modify element O H`). |
| `dump.*.extxyz` | **Extended XYZ** (auto after run) — **`Lattice=`** + **`species`** so Ovito gets PBC + correct elements. |

If you only have `.lammpstrj`:  
`python lammpstrj_to_extxyz.py dump.ice_uma_quick.lammpstrj data.ice_uma`

## Energy minimization

- **`minimize`** runs **after** the UMA `fix external` callback is registered (runner splits the input at the first `minimize` line). Same idea as OpenMM’s prerun minimize.
- Quick deck: `minimize 0.0 1.0e-6 500 10000` (force tol ~1e-6 eV/Å). Full deck uses more iterations.

## Thermostat / temperature

- **`fix nvt` + ML external forces** often shows **huge early T spikes** (e.g. 243 K → 1000+ K): stiff UMA forces + tight Nose–Hoover coupling, not a wrong box.
- Inputs now use **`fix nve` + `fix langevin 243 243 1.0`** (damping ~1 ps), like many MLIP workflows; **T should stay near 243 K** after a short transient.
- **Box** is correct: same reduced cell as OpenMM (`build_lammps_ice_data.py`). **`--molecules 8`** is a **tiny** cell (~6×6×7 Å) — larger T noise than full **128-molecule** ice; use full `data.ice_uma` for production.

## Difference vs OpenMM run

| OpenMM | This LAMMPS setup |
|--------|-------------------|
| RPMD (many beads) | **Classical** MD (1 bead) — same UMA forces. To match cost in OpenMM use **`--beads 1`** (see main `../README.md` “Speed vs LAMMPS”). |
| PILE | **Langevin + NVE** (~NVT sampling) |
| PDB trajectory | `dump custom` → `.lammpstrj` |

RPMD in LAMMPS needs extra packages (`fix pimd`, etc.); this deck is the **direct classical analogue**
of a **1-bead** run so forces and geometry match.

## Prerequisites

1. **LAMMPS with Python interface** — `pip install lammps` is enough. The wheel links **MPICH** (`libmpi.so.12`). Conda **OpenMPI** only provides `libmpi.so.40`, so you still get *cannot open libmpi.so.12* unless you add MPICH:
   ```bash
   conda install -y -c conda-forge mpich
   export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
   python -c "from lammps import lammps; lammps(cmdargs=['-log','none']); print('OK')"
   ```
   Or use the helper script (edit `fenics` path if needed): `./run_uma_ice.sh`

2. **FAIRChem + fairchem-lammps**
   ```bash
   pip install fairchem-core fairchem-lammps
   ```

3. **GPU** (recommended): UMA on CPU is very slow for 128 molecules.

## Run

```bash
cd tests/uma_ice_rpmd/lammps

# Data file from same GenIce CIF as OpenMM (default ../ice.cif)
python build_lammps_ice_data.py -i ../ice.cif -o data.ice_uma

# Full 128 molecules (~55 ps @ 1 fs as in in.ice_uma.lmp)
python run_lammps_uma_ice.py

# Smaller / faster
python run_lammps_uma_ice.py --molecules 16
```

Outputs (under this directory):

- **`dump.ice_uma.lammpstrj`** — multi-frame trajectory (Ovito: File → Open)
- Edit **`in.ice_uma.lmp`** to change `run`, `dump` frequency, temperature (243 K), timestep (0.001 ps = 1 fs)

## Input script rules (fairchem-lammps)

- **Do not** use `pair_style`, `bond_style`, etc. ML forces replace all of that.
- **`units metal`**: Å, eV, ps, eV/Å forces internally in the callback.
- Atom **masses** must be physical (O/H); the callback maps mass → atomic number.

## License note

`fairchem-lammps` is **GPL-2.0** (LAMMPS-compatible). The rest of this repo may use other licenses; only this workflow pulls in GPL code when you install/run the driver.
