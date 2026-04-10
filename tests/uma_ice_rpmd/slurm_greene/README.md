# SLURM RPMD Benchmarks for NYU HPC Greene

Submit RPMD simulations of TIP4P, OpenMM UMA, and i-PI+LAMMPS UMA across bead counts 8–32 and 100 replicas.

## Configuration

| Parameter | Value |
|-----------|-------|
| Timestep | 0.1 fs |
| Production | 100 ps |
| Beads | 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32 |
| Models | TIP4P, OpenMM UMA, LAMMPS UMA |
| Replicas | 100 |
| Molecules | 128 |
| Order params | Q6, q_tet from centroid bead positions |

Total jobs: 3 × 13 × 100 = 3900 (array 0–3899).

## Prerequisites

1. **Partition/GPU**: NYU Greene `gpu` partition with 1 GPU. Check with `sinfo -p gpu`. Use `rtx8000` or another GPU partition if required.
2. **Modules**: `anaconda3`, `cuda`. Adjust in `submit_slurm.sh` if your cluster uses different names.
3. **Conda env**: Environment `openmm` (or similar) with:
   - openmm, openmmml
   - fairchem, fairchem-lammps, fairchem-core
   - ipi
   - lammps (built with MISC for fix ipi)
   - numpy, scipy, ase

4. **Data**: Build LAMMPS ice data once:
   ```bash
   cd tests/uma_ice_rpmd
   python lammps/build_lammps_ice_data.py -n 128 -o lammps/data.ice_uma_128
   ```

## Usage

```bash
cd tests/uma_ice_rpmd/slurm_greene
./submit_all.sh
```

Or submit manually:
```bash
mkdir -p logs slurm_out
sbatch --job-name=rpmd_bench --array=0-3899 --time=48:00:00 submit_slurm.sh
```

## Output Layout

```
slurm_out/
  tip4p/replica_000/bead_08/order.csv
  tip4p/replica_000/bead_10/order.csv
  ...
  openmm_uma/replica_000/bead_08/order.csv
  ...
  lammps_uma/replica_000/bead_08/order.csv
  ...
```

Order CSVs use the same schema (step, time_ps, q6_mean, q_tet_mean, …) for downstream plotting.

## Single-Task Test

```bash
SLURM_ARRAY_TASK_ID=0 bash run_rpmd_single.sh
```

Runs one TIP4P replica 0, bead 8 job locally.

## LAMMPS UMA Isolation

i-PI uses a socket port; array tasks can land on the same node. Each LAMMPS UMA task:
- Uses `--port $((64511 + (id % 10000)))` for a unique port
- Copies the project to `$TMPDIR/rpmd_${JOB_ID}_${TASK_ID}` and runs there to avoid overwriting `ipi/` outputs

## Tuning

- **Partition**: Edit `#SBATCH --partition=gpu` in `submit_slurm.sh`.
- **Walltime**: Default 48 h; increase if 100 ps × 128 mol × 32 beads is slow.
- **Memory**: 48 GB; reduce for smaller bead counts if needed.
- **Conda env**: Set in `submit_slurm.sh` (`conda activate openmm`).
