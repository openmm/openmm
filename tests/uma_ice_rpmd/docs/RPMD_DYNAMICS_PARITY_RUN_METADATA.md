# RPMD dynamics parity — run metadata (fill when saving figures)

Copy this block into your lab notebook or rename to `pipeline_out/RUN_<date>_metadata.txt` when you generate comparison plots. **Dynamics** (Q6, q_tet, centroid T) can diverge even when **static forces** agree if any of these differ.

## Template (copy and fill)

```
Date:
Figure output path:
Git commit:

--- OpenMM leg ---
Command (exact):
  (e.g. bash run_rpmd_32x32.sh  OR  python run_openmm_rpmd_reference.py ...)

Environment:
  CONDA_ENV=
  OPENMM_PLUGIN_DIR=
  OPENMMML_UMA_RPMD_CHUNK=   # e.g. 4
  PYTORCH_CUDA_ALLOC_CONF=
  RPMD_STEPS=
  RPMD_DT_FS=
  RPMD_EQUIL_PS=             # 0 vs 2 strongly affects early centroid T_K
  RPMD_REPORT_EVERY_STEPS=

OpenMM / test_uma_ice_rpmd:
  --beads
  --dt (fs)
  --equil (ps)
  --temperature (K)
  --seed
  --rpmd-thermostat
  --rpmd-friction (1/ps)
  --rpmd-centroid-friction (1/ps)
  --platform
  --model

Output CSV:
  pipeline_out/ice_order_openmm_rpmd.csv  (or path)

--- i-PI + client leg ---
Command:
  (e.g. python run_ipi_lammps_uma_rpmd.py ...)

RPMD_IPI_CLIENT=   # lammps | python

ipi/input.xml snapshot (or attach):
  total_steps:
  timestep (fs):
  thermostat:
  tau (fs):
  nbeads:
  seed:

Output CSV:
  pipeline_out/ice_order_ipi_rpmd.csv

--- Post-processing ---
  ipi_order_from_traj: --thermo default uses temperature(nm=0) for T_K when present
  plot: python plot_rpmd_comparison.py ...
```

## Example: prior “fast” 32×32 pipeline (asymmetric equil)

Settings that produced **OpenMM with no extra NVT equilibration** while i-PI ran a single continuous segment (see [`run_rpmd_32x32.sh`](../run_rpmd_32x32.sh) history):

| Variable | Typical value |
|----------|----------------|
| `RPMD_STEPS` | 10000 |
| `RPMD_DT_FS` | 0.1 |
| `RPMD_EQUIL_PS` | **0** (OpenMM production CSV time starts after minimization + velocity init only) |
| `OPENMMML_UMA_RPMD_CHUNK` | 4 |
| i-PI `tau` | 1000 fs |
| OpenMM centroid / internal friction | 0.5 / 1.0 1/ps |

For **closer thermostat-time alignment** with [`run_openmm_rpmd_reference.py`](../run_openmm_rpmd_reference.py) alone, use **`RPMD_EQUIL_PS=2`** (this is now the default in [`run_rpmd_32x32.sh`](../run_rpmd_32x32.sh) and [`run_openmm_uma_rpmd_only.sh`](../run_openmm_uma_rpmd_only.sh)).

**After changing equil defaults**, regenerate both legs and the figure, e.g. `bash run_rpmd_32x32.sh`, then inspect `pipeline_out/rpmd_comparison_32x32.png`. Order-parameter code uses ``scipy.special.sph_harm_y`` on SciPy ≥ 1.15 (and ``sph_harm`` on older releases); see ``ice_order_parameters._sph_harm_batch``.

## Related scripts

- [`sweep_openmm_rpmd_thermostat.sh`](../sweep_openmm_rpmd_thermostat.sh) — short OpenMM-only runs varying `--rpmd-friction` / `--rpmd-centroid-friction`.
- [`plot_ipi_thermo.py`](../plot_ipi_thermo.py) — i-PI extended `temperature` vs `temperature(nm=0)` vs time.
