# i-PI `pile_g` vs OpenMM PILE-G (reference driver)

This note summarizes **what is wired in the repo** and what you should verify when comparing `ice_order_ipi_rpmd.csv` to `ice_order_openmm_rpmd.csv`.

## Dynamics parity checklist (before trusting Q6 / T_K plots)

Static UMA forces can agree while trajectories differ. Work through:

1. **Log what you ran** — copy the template in [`RPMD_DYNAMICS_PARITY_RUN_METADATA.md`](RPMD_DYNAMICS_PARITY_RUN_METADATA.md) (equil, dt, beads, friction, `OPENMMML_UMA_RPMD_CHUNK`, i-PI `tau`, seed).
2. **Match OpenMM equilibration policy** — [`run_rpmd_32x32.sh`](../run_rpmd_32x32.sh) keeps the same explicit `RPMD_EQUIL_PS` policy as the parity run you request; setting **`RPMD_EQUIL_PS=0`** prints a stderr warning because early centroid `T_K` vs i-PI often diverges.
3. **Same centroid temperature definition** — OpenMM `T_K` is centroid kinetic T. i-PI’s merged CSV should use **`temperature(nm=0)`** for `T_K` when [`ipi_order_from_traj.py`](../ipi_order_from_traj.py) is run with thermo (see below); extended `temperature` is not comparable.
4. **Thermostat knobs** — parity runners now treat i-PI `tau` as the reference input (`RPMD_PILE_TAU_FS`) and derive the OpenMM centroid coupling from it. Internal PILE damping remains the normal-mode critical damping used by `pile_g`.
5. **Optional diagnostics** — [`plot_ipi_thermo.py`](../plot_ipi_thermo.py) overlays i-PI extended vs centroid temperature vs time.

## Equilibration

| Path | Default behaviour |
|------|-------------------|
| [`run_openmm_rpmd_reference.py`](../run_openmm_rpmd_reference.py) | **`--equil 2.0` ps** of ring-polymer NVT before production order rows (unless overridden). |
| [`run_rpmd_32x32.sh`](../run_rpmd_32x32.sh) | Uses the explicit `RPMD_EQUIL_PS` you request. `RPMD_EQUIL_PS=0` is still supported for smoke/parity starts and prints a warning. |
| i-PI + LAMMPS | No separate equil block in the generated XML; dynamics start at `init.xyz`. |

**Action:** Check which `--equil` value produced your OpenMM CSV (`head` the file’s implied timeline: `time_ps` at the first production row vs `step`).

## Thermostat mapping (friction ↔ τ)

OpenMM ([`test_uma_ice_rpmd.py`](../test_uma_ice_rpmd.py)) now treats the centroid thermostat as the exact i-PI `ThermoSVR` update used inside `pile_g`, with the API still expressed as:

- **PILE internal modes:** critical damping from the ring-polymer normal-mode frequencies.
- **Centroid (PILE-G):** `--rpmd-centroid-friction` in 1/ps, interpreted as **`1/tau_ps`** for the i-PI-equivalent centroid SVR step.

i-PI `pile_g` uses a single **`<tau>`** in femtoseconds ([`_write_ipi_input`](../run_ipi_lammps_uma_rpmd.py), default **1000 fs**). The parity shell runners pass the same `RPMD_PILE_TAU_FS` to i-PI and convert it to the OpenMM centroid friction.

Parameter mapping:

- centroid damping time τ (ps) = **1 / (centroid friction in 1/ps)**.

So i-PI `tau = 1000 fs` corresponds to OpenMM centroid friction **1.0 /ps** in the parity runners.

## Kinetic temperature column

OpenMM `T_K` is **centroid kinetic temperature** (massive DoF only); see [`KINETIC_TEMPERATURE_UMA_RPMD.md`](../KINETIC_TEMPERATURE_UMA_RPMD.md).

i-PI exports in [`run_ipi_lammps_uma_rpmd.py`](../run_ipi_lammps_uma_rpmd.py) `<properties>`:

- `temperature` — extended ring-polymer MD temperature (typically **much** larger than 243 K).
- `temperature(nm=0)` — normal-mode 0 / centroid kinetic temperature — **use this** when filling `T_K` next to OpenMM in [`plot_rpmd_comparison.py`](../plot_rpmd_comparison.py).

**Confirmed wiring:** [`ipi_thermo_utils.build_thermo_lookup`](../ipi_thermo_utils.py) sets `T_K` from **`temperature(nm=0)`** when that column exists, otherwise falls back to `temperature`. [`ipi_order_from_traj.py`](../ipi_order_from_traj.py) merges thermo into the order CSV using that helper (default `--thermo` / `ice__i-pi.md` when present), so [`plot_rpmd_comparison.py`](../plot_rpmd_comparison.py)’s third panel compares like-with-like **as long as** the properties file includes `temperature(nm=0)` (the template in `run_ipi_lammps_uma_rpmd.py` requests it).

[`plot_ipi_thermo.py`](../plot_ipi_thermo.py) plots both i-PI temperatures for diagnostics when extended vs centroid must be inspected side-by-side.
