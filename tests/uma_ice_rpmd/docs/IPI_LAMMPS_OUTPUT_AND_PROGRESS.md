# i-PI + LAMMPS: what “progress” means (and why LAMMPS looks stuck)

This note closes the loop on **misleading LAMMPS output** versus **real simulation progress** when using [`run_ipi_lammps_uma_rpmd.py`](../run_ipi_lammps_uma_rpmd.py) with `fix ipi` + `fix external` (UMA).

## What the manuals say

### LAMMPS `fix ipi`

From the [LAMMPS `fix ipi` documentation](https://docs.lammps.org/fix_ipi.html):

- i-PI is the **driver** for **advanced molecular simulations** (including PIMD/RPMD).
- **Philosophy:** separate **energy/force evaluation** (client) from **dynamics** (i-PI).
- **Communication** is minimal (sockets); **all parameters of the dynamics are specified in i-PI’s input**.
- **Force-field parameters** must be in the LAMMPS input **before** `fix ipi`.
- **The initial configuration is ignored** — it is **replaced by coordinates received from i-PI** before forces are evaluated.

So LAMMPS is **not** running a standalone Verlet trajectory in the usual sense; **thermo lines are not a direct “i-PI step” counter.**

### i-PI overview

From the [i-PI introduction](https://docs.ipi-code.org/introduction.html):

- i-PI uses a **client–server** model: i-PI **commands the evolution** of nuclear positions; clients compute properties (forces, etc.).
- Communication is via **sockets** (Internet or Unix domain).
- **Output** (trajectories, restarts) is handled by i-PI’s **output** machinery.

**Implication:** **Step counts, timing, and trajectory strides** should be read from **i-PI** (`input.xml`, console / `i-pi_run.log`, `ice__i-pi.traj_*.xyz`), not from LAMMPS thermo.

## How our wiring behaves

| Item | Where it lives | Consistency check |
|------|----------------|-------------------|
| Total MD steps | `<total_steps>` in `ipi/input.xml` (written by `_write_ipi_input`) | Must match LAMMPS `run N` — the script replaces the template `run` line with `--steps`. |
| Integration timestep | `<timestep units="femtosecond">` in XML | **Authoritative** for the nuclear/RPMD integration in i-PI. |
| Trajectory frequency | `<trajectory ... stride="S"/>` | Expect about `floor(total_steps / S)` **frames per bead file** (e.g. 2000 steps, stride 10 → ~200 frames if the run completes). |
| LAMMPS `Time step` / thermo | Printed by LAMMPS (`metal` default e.g. 0.001 ps) | **Not** the same object as i-PI’s timestep; LAMMPS docs state **dynamics parameters are in i-PI**. Do not infer i-PI dt from this line. |
| `WARNING: No fixes with time integration` ([LAMMPS err0028](https://docs.lammps.org/err0028.html)) | LAMMPS `run` with `fix ipi` | **Expected** in spirit: there is no local `fix nve`/`fix nvt` **driving** atoms — **i-PI** drives positions via the socket protocol. |

## Where to look for “real” progress

1. **`ipi/i-pi_run.log`** (default; see `--ipi-log`) — i-PI banner, input echo, **socket bind**, **client handshake** (`Handshaking was successful`), **force-field polling** (`Starting the polling thread main loop`). This is the right place to confirm the driver is alive and the LAMMPS client connected.
2. **`ipi/ice__i-pi.traj_*.xyz`** — **Bead trajectories** written by i-PI; count **frames** (lines that are only the atom count, e.g. `96`) to measure progress **independent of LAMMPS thermo**.
3. **`ipi/i-pi.RESTART`** — restart data from i-PI.
4. **LAMMPS screen** — useful for **debugging the deck** (`read_data`, `pair_style`, etc.), **not** as the primary progress meter under `fix ipi`.

## Why LAMMPS thermo can sit at “Step 0”

- Under `fix ipi`, LAMMPS is blocked in the **client loop** servicing i-PI’s force requests; **standard thermo may not advance** in a way that tracks **i-PI’s** internal step index.
- **Cost per i-PI step** here is dominated by **UMA** (`fix external`) over **many beads** — wall time can be **large** before the **next trajectory frame** (e.g. stride 10) appears, so **one frame** in `traj_0.xyz` does **not** mean “nothing happened”; it may mean **steps 1–9** completed without a dump.

## Operational pitfalls we observed

- **Port already in use** — a leftover `i-pi` process can leave `OSError: [Errno 98] Address already in use`. Use a free **`--port`** or stop old i-PI processes before re-running.
- **Subprocess timeout / SIGTERM** — killing the wrapper sends **SOFTEXIT** in i-PI’s log (`Kill signal received`); trajectories may show **only early frames**.

## References

- LAMMPS: [fix ipi](https://docs.lammps.org/fix_ipi.html)
- i-PI: [Program overview](https://docs.ipi-code.org/introduction.html), [FAQ](https://docs.ipi-code.org/faq.html)
- Ceriotti et al., *Comput. Phys. Commun.* **185**, 1019 (2014); Kapil et al., *Comput. Phys. Commun.* **236**, 214 (2019)
