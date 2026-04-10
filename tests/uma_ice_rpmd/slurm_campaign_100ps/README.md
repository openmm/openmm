# SLURM campaign: 100 ps RPMD, size × bead × replica grid

**Agent / handoff summary:** see [`AGENT_SLURM_CAMPAIGN.md`](AGENT_SLURM_CAMPAIGN.md) for a compact description for automation or the next assistant.

Independent **three-way** submissions (TIP4P OpenMM, OpenMM UMA, i-PI+LAMMPS UMA), each as a **job array 0–2499** (2500 tasks).

## Parameter grid

| Dimension | Values |
|-----------|--------|
| Cubic supercell | **1×1×1 … 5×5×5** → molecules \(N = 8 n^3\) → **8, 64, 216, 512, 1000** |
| RPMD beads | **4, 8, 16, 32, 64** |
| Replicas | **100** (indices 0–99) |
| Production | **100 ps** at **dt = 0.1 fs** → **1,000,000 steps** (override with `RPMD_STEPS` or `PROD_PS`; see below) |

**Tasks per engine:** \(5 \times 5 \times 100 = 2500\). Submitting all three engines yields **7500** jobs total.

### Array index decode (`SLURM_ARRAY_TASK_ID` = `ID`)

```text
replica   = ID % 100              # 0..99
rest      = ID / 100              # 0..24
bead_idx  = rest % 5              # 0..4  →  beads = (4, 8, 16, 32, 64)
size_idx  = rest / 5              # 0..4  →  edge n = 1..5
N         = 8 * n³
```

RNG seed: **`SEED = CAMPAIGN_SEED_BASE + ID`** (default base **284759**).

## Prerequisites

1. **Conda env** with OpenMM, (for UMA) openmmml/fairchem, (for i-PI path) `i-pi`, LAMMPS + MISC, fairchem-lammps — same as [../slurm_greene/README.md](../slurm_greene/README.md).
2. **LAMMPS ice data** for all five sizes (once per machine):

   ```bash
   cd tests/uma_ice_rpmd
   bash slurm_campaign_100ps/prepare_lammps_data.sh
   ```

   Writes `lammps/data.ice_uma_8`, `_64`, `_216`, `_512`, `_1000`.

3. Edit **`#SBATCH`** lines in `submit_*.slurm` (**partition**, **mem**, **`--time`**, module names) for your cluster. Default **48 h** may be insufficient for the largest **(1000 molecules × 64 beads)** UMA jobs — increase walltime or split the array by size.

## Resource warnings

- **OpenMM UMA** at large \(N\) and **64 beads** can **OOM** on consumer GPUs. Use `OPENMMML_UMA_RPMD_CHUNK` (set in `submit_openmm_uma.slurm`, default **4**) and `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`.
- **100 ps** with **1e6 steps** is expensive for all three engines; benchmark one `(n, beads)` pair locally before full campaign.
- **i-PI + LAMMPS** tasks copy the `tests/uma_ice_rpmd` tree to `$TMPDIR` (excluding `slurm_campaign_100ps/slurm_out/` when `rsync` is available) to avoid socket/`ipi/` collisions between array tasks.

## Submit (three separate jobs)

From `tests/uma_ice_rpmd/slurm_campaign_100ps/`:

```bash
mkdir -p logs slurm_out
sbatch submit_tip4p.slurm
sbatch submit_openmm_uma.slurm
sbatch submit_lammps_uma.slurm
```

## Output layout

```text
slurm_out/
  tip4p/n{n}_b{beads}/replica_{rrr}/order.csv
  tip4p/.../task_meta.txt
  openmm_uma/...
  lammps_uma/...
```

`order.csv` uses the path-integral ice order schema (see `ice_order_parameters.ICE_ORDER_PI_CSV_HEADER`).

## Environment overrides

| Variable | Meaning |
|----------|---------|
| `RPMD_STEPS` | If set, passed as `--steps` (overrides `PROD_PS` length). Use for short tests (e.g. `1000`). |
| `PROD_PS` | Production time (ps) if `RPMD_STEPS` unset (default **100**). |
| `RPMD_DT_FS` | Timestep in fs (default **0.1**). |
| `ORDER_EVERY` | Order CSV stride for TIP4P / OpenMM UMA (default **100**). |
| `CAMPAIGN_SEED_BASE` | Base for `SEED` (default **284759**). |

## Local smoke test (no SLURM)

After `prepare_lammps_data.sh` and from `tests/uma_ice_rpmd`:

```bash
export RPMD_STEPS=200
export CAMPAIGN_MODEL=tip4p
SLURM_ARRAY_TASK_ID=0 bash slurm_campaign_100ps/run_campaign_task.sh
```

Uses task 0 → \(n=1\), **4 beads**, replica **0** → needs `lammps/data.ice_uma_8`.

## Files

| File | Role |
|------|------|
| `prepare_lammps_data.sh` | Build all `data.ice_uma_*` files for \(n=1..5\). |
| `run_campaign_task.sh` | Dispatcher: decode `SLURM_ARRAY_TASK_ID`, run one model. |
| `submit_tip4p.slurm` | Array **0–2499**, `CAMPAIGN_MODEL=tip4p`. |
| `submit_openmm_uma.slurm` | Array **0–2499**, `CAMPAIGN_MODEL=openmm_uma`. |
| `submit_lammps_uma.slurm` | Array **0–2499**, `CAMPAIGN_MODEL=lammps_uma`. |
