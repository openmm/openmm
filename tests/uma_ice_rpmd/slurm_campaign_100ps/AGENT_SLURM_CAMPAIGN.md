# AI agent note: `slurm_campaign_100ps` RPMD campaign

**Purpose:** This file gives a future assistant enough context to **operate on or extend** the SLURM campaign **without re-reading the whole repo**. Human-oriented detail lives in [`README.md`](README.md).

## What this campaign is

- **Location:** `tests/uma_ice_rpmd/slurm_campaign_100ps/`
- **Goal:** Run **path-integral RPMD** ice benchmarks at **100 ps** (default **0.1 fs** → **1,000,000 steps**) on a fixed grid of **system sizes**, **bead counts**, and **replicas**, using **three separate physics stacks** so they do not block each other on the scheduler.

## Fixed parameter grid (do not assume different defaults)

| Axis | Values |
|------|--------|
| Cubic supercell edge | `n ∈ {1,2,3,4,5}` → molecules `N = 8 n³` → **8, 64, 216, 512, 1000** |
| RPMD beads | **4, 8, 16, 32, 64** |
| Replicas | **100** per `(n, beads)` combination |
| SLURM array | **0–2499** per engine (**2500 tasks** each) |

**Total jobs if all three engines are submitted:** `3 × 2500 = 7500`.

## Three independent submissions (critical)

Engines **must** be submitted separately (three `sbatch` calls). Each uses the same array index logic but **`CAMPAIGN_MODEL`** selects the driver:

| Submit script | `CAMPAIGN_MODEL` | Driver (from `tests/uma_ice_rpmd/`) |
|----------------|------------------|--------------------------------------|
| `submit_tip4p.slurm` | `tip4p` | `run_openmm_tip4p_rpmd.py` |
| `submit_openmm_uma.slurm` | `openmm_uma` | `run_openmm_rpmd_reference.py` → `test_uma_ice_rpmd.py` |
| `submit_lammps_uma.slurm` | `lammps_uma` | `run_ipi_lammps_uma_rpmd.py` + `ipi_order_from_traj.py` |

**Dispatcher:** `run_campaign_task.sh` (sourced by each `.slurm` via `bash run_campaign_task.sh`).

## Array index → physics parameters

For `ID = SLURM_ARRAY_TASK_ID`:

```text
replica   = ID % 100
rest      = ID // 100          # integer division, 0..24
bead_idx  = rest % 5           → beads ∈ {4,8,16,32,64}
size_idx  = rest // 5          → edge n = size_idx + 1 ∈ {1..5}
N_mol     = 8 * n³
SEED      = CAMPAIGN_SEED_BASE + ID   (default base 284759)
```

**Output path pattern:**

```text
slurm_out/<tip4p|openmm_uma|lammps_uma>/n<n>_b<beads>/replica_<rrr>/order.csv
slurm_out/.../task_meta.txt
```

## One-time data prep

Ice LAMMPS data files must exist **before** array tasks run:

```bash
cd tests/uma_ice_rpmd
bash slurm_campaign_100ps/prepare_lammps_data.sh
```

Produces `lammps/data.ice_uma_{8,64,216,512,1000}`. The OpenMM/UMA reference path expects this naming (`data.ice_uma_<N>`).

## Environment overrides (for tests or partial reruns)

| Variable | Effect |
|----------|--------|
| `RPMD_STEPS` | If set → `--steps` (short dry runs, e.g. `200`) |
| `PROD_PS` | Production length (ps) if `RPMD_STEPS` unset (default **100**) |
| `RPMD_DT_FS` | Timestep fs (default **0.1**) |
| `ORDER_EVERY` | Order CSV stride for TIP4P / OpenMM UMA (default **100**) |
| `CAMPAIGN_SEED_BASE` | Seed offset base (default **284759**) |
| `OPENMMML_UMA_RPMD_CHUNK` | UMA batched bead chunk size (set in `submit_openmm_uma.slurm`, default **4**) |

**Local smoke test (no SLURM):** see [`README.md`](README.md) “Local smoke test”.

## Operational pitfalls (check these first)

1. **Walltime:** Default `#SBATCH --time` may be too short for **large N × 64 beads**, especially **OpenMM UMA**. Increase time or split arrays by `ID` range.
2. **GPU OOM (UMA):** Raise `OPENMMML_UMA_RPMD_CHUNK`, set `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` (already exported in `submit_openmm_uma.slurm`).
3. **LAMMPS/i-PI:** Each task uses a **unique port** and a **fresh copy** of the tree under `$TMPDIR` so `ipi/` sockets and trajectories do not collide. If `rsync` is missing, the script falls back to `cp` and deletes `slurm_out` inside the copy target—ensure `rsync` exists on compute nodes for large trees.
4. **Trajectory name:** Post-processing looks for `ipi/ice__i-pi.traj_0.xyz`, then `traj_00`, then `ice_traj.xyz` (same idea as `slurm_greene/run_rpmd_single.sh`).

## How to extend or change the campaign

- **Different sizes/beads/replicas:** Edit the decode logic and array size in `run_campaign_task.sh` and the `#SBATCH --array=` line in all three `submit_*.slurm` files; update this note and `README.md`.
- **Different production length:** Prefer `PROD_PS` / `RPMD_STEPS` via export in `.slurm` or user env; document for users.
- **Do not confuse** with `tests/uma_ice_rpmd/slurm_greene/` — that is an **older** benchmark (bead sweep 8–32, 128 molecules, different array mapping). The **100 ps grid campaign** is **only** under `slurm_campaign_100ps/`.

## File map (quick)

| File | Role |
|------|------|
| `prepare_lammps_data.sh` | Build all five `data.ice_uma_*` files |
| `run_campaign_task.sh` | Single-task logic; **all engines** |
| `submit_tip4p.slurm` | SLURM wrapper, TIP4P only |
| `submit_openmm_uma.slurm` | SLURM wrapper, OpenMM UMA only |
| `submit_lammps_uma.slurm` | SLURM wrapper, i-PI+LAMMPS only |
| `README.md` | Human runbook |
| `AGENT_SLURM_CAMPAIGN.md` | **This** agent-oriented summary |

---

*Last aligned with the campaign layout as of creation: three arrays × 2500 tasks, 100 ps default, 5×5×100 grid.*
