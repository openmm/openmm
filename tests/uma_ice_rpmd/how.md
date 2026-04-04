# RPMD simulations with UMA potential
_Exported on 3/26/2026 at 13:38:39 EDT from Cursor (2.6.19)_

---

**User**

We are in the mlip branch, running RPMD simulations using the UMA potential. Review the LAMMPS implementation and local implementation in this repo and report to me how it works. I think we even have like an independent driver with i-pi?

---

**Cursor**


Exploring the repository for LAMMPS RPMD/UMA setup, local RPMD implementations, and i-pi driver integration.




Here is a concise report of how UMA + RPMD (and related benchmarks) are wired in this repo.

## 1. OpenMM (“local” implementation)

**Dynamics:** `tests/uma_ice_rpmd/test_uma_ice_rpmd.py` builds an ice Ih system and runs **`openmm.RPMDIntegrator`** (PILE / barostat options as in that script) at the chosen temperature and bead count.

**Potential:** **`openmmml.MLPotential`** with the small batched UMA model (README calls out `uma-s-1p1-pythonforce-batch`). Forces come from the Python/Torch path in `wrappers/python/openmmml/models/umapotential_pythonforce_batch.py`, which batches **all RPMD beads in one (or chunked) ML call(s)** so each MD step does one UMA evaluation over the polymer, not one independent classical step.

**C++ side:** The RPMD plugin (`plugins/rpmd/…`) implements ring-polymer evolution; there is separate work (see `FIXES_SUMMARY.md`) on **additive** force handling and hybrid RPMD so UMA + other forces stay consistent.

So: **OpenMM = native RPMD integrator + UMA as a custom force that knows about bead batching.**

---

## 2. LAMMPS implementation (ice test folder)

**Role in the README:** The standard LAMMPS decks under `tests/uma_ice_rpmd/lammps/` are documented as **classical MD with one configuration per step**—the **cost analogue of OpenMM with `--beads 1`**, not full RPMD inside LAMMPS.

**How UMA is attached:** `run_lammps_uma_ice.py` uses **FAIRChem’s `fairchem-lammps`** driver:

- Loads the same ML stack: `pretrained_mlip.get_predict_unit(model, device=...)`.
- Registers LAMMPS **`fix external`** via `FIX_EXTERNAL_CMD` and `FixExternalCallback` (charge/spin for `omol`).
- Splits the input so **`minimize` runs only after** the external fix is active (UMA forces present during minimization).

```265:274:tests/uma_ice_rpmd/lammps/run_lammps_uma_ice.py
    machine = os.environ.get("LAMMPS_MACHINE_NAME")
    lmp = lammps(name=machine, cmdargs=["-nocite", "-log", log_arg, "-echo", "screen"])
    lmp._predictor = pretrained_mlip.get_predict_unit(args.model, device=args.device)
    lmp._task_name = "omol"

    lmp.commands_list([x for x in pre if x.strip()])
    lmp.command(FIX_EXTERNAL_CMD)
    lmp.set_fix_external_callback(FIX_EXT_ID, FixExternalCallback(charge=args.charge, spin=args.spin), lmp)
    print("--- UMA fix external on → minimize + MD ---")
    lmp.commands_list([x for x in post if x.strip() or x == ""])
```

Input scripts use **`pair_style zero`** (no classical FF); all interatomic forces are from the callback—same idea as README: **LAMMPS and OpenMM both end up calling Fairchem’s `predict(AtomicData)`-style path.**

---

## 3. i-PI: yes, and it’s mostly **not** LAMMPS anymore in the orchestrator

There **is** an i-PI workflow, but the main driver script’s **current** design is:

1. **`i-pi`** reads a generated `tests/uma_ice_rpmd/ipi/input.xml` (RPMD / `nbeads`, NVT, PILE-style thermostat, inet socket).
2. **`i-pi-py_driver`** connects as the force client with **`-P tests/uma_ice_rpmd/ipi/uma_ipi_driver.py -m custom`**.

So the **force engine is a small custom i-PI PES driver** (`UMADriver`), which subclasses `ipi.pes.ase.ASEDriver`, builds a **`FAIRChemCalculator`** from `pretrained_mlip.get_predict_unit`, and returns energy/forces (and a zero virial) in i-PI’s internal units:

```97:123:tests/uma_ice_rpmd/ipi/uma_ipi_driver.py
    def check_parameters(self):
        """Load template structure, fix cell, and attach FAIRChemCalculator."""
        from fairchem.core import FAIRChemCalculator
        from fairchem.core.calculate import pretrained_mlip

        super().check_parameters()

        # ASE's xyz reader does not parse i-PI's CELL(abcABC) comment; fix it.
        cell = _read_ipi_cell_angstrom(self.template)
        if cell is not None:
            self.template_ase.set_cell(cell)
            print(f"[UMADriver] Cell set from i-PI comment: {np.diag(cell)} Å", flush=True)
        else:
            print("[UMADriver] WARNING: could not parse cell from template xyz", flush=True)

        self.template_ase.pbc = [True, True, True]
        self.template_ase.info["charge"] = 0
        self.template_ase.info["spin"] = 1

        print(f"[UMADriver] Loading model={self.model} on {self.device} ...", flush=True)
        predict_unit = pretrained_mlip.get_predict_unit(
            self.model, device=self.device
        )
        self.ase_calculator = FAIRChemCalculator(
            predict_unit, task_name=self.task_name
        )
        print(f"[UMADriver] FAIRChemCalculator ready (task={self.task_name})", flush=True)
```

**Geometry source:** `run_ipi_lammps_uma_rpmd.py` still uses **`build_lammps_ice_data.py`-style data** and `convert_lammps_to_ipi_xyz.py` only to produce **`init.xyz`** for i-PI—the name **`run_ipi_lammps_uma_rpmd.py`** is legacy; it does **not** spawn LAMMPS in the current code path.

**XML naming:** The socket block is still labeled `name="lammps"` / `forcefield="lammps"` in the generated XML—**historical label**; the connected client is the Python UMA driver, not LAMMPS.

**Alternate (unused by that script today):** `tests/uma_ice_rpmd/lammps/in.ice_uma_ipi_client.lmp` shows **LAMMPS as i-PI client** with `fix ipi` + commented injection of **`fix external` (UMA)**—the pattern i-PI docs describe for a socket client, but your orchestrator was switched to **inet + `i-pi-py_driver`** to avoid Unix-socket / handshake issues (as noted in `chat.md`).

---

## 4. How the three compare conceptually

| Piece | Who integrates RPMD? | Who evaluates UMA? |
|--------|----------------------|---------------------|
| **OpenMM** | `RPMDIntegrator` + RPMD plugin | `openmmml` batched UMA (`TorchForce` / PythonForce path) |
| **LAMMPS (default ice script)** | None (1 bead / classical) | `fix external` + `FixExternalCallback` (fairchem-lammps) |
| **i-PI + `uma_ipi_driver.py`** | i-PI (ring polymer, thermostat in XML) | `FAIRChemCalculator` in the custom ASE driver |

So: **same underlying UMA stack (Fairchem)** in all three; **RPMD** is native in OpenMM, external in i-PI, and **not** the default in the main LAMMPS ice driver (that’s intentional for a classical force/geometry match and cost comparison).

If you want, we can next map **thermostat / timestep / bead conventions** line-by-line between `input.xml` / `_write_ipi_input` and `RPMDIntegrator` settings in `test_uma_ice_rpmd.py` (that’s where subtle “they don’t match” issues usually hide).

---

**User**

How to run a case that reproduces @tests/uma_ice_rpmd/pipeline_out/rpmd_comparison.png  but for the i+PI version?

---

**Cursor**


Tracing how `rpmd_comparison.png` is produced and which i-PI steps mirror it.




Here’s how to reproduce **`tests/uma_ice_rpmd/pipeline_out/rpmd_comparison.png`**, but with the **i‑PI / UMA RPMD** side filled in the way this repo expects.

### What your PNG actually is

The default **`plot_rpmd_comparison.py`** title talks about **i‑PI+LAMMPS UMA**, but each curve is chosen by **which CSVs you pass**. Your figure’s **green** series is labeled **TIP4P/2005f RPMD** in the script’s legend—that means that plot was built with **`--tip4p pipeline_out/ice_order_tip4p_rpmd.csv`**, not with an i‑PI CSV. To get the **i‑PI** version, you pass **`--ipi pipeline_out/ice_order_ipi_rpmd.csv`** and omit **`--tip4p`** (unless you want all three).

**`run_ipi_lammps_uma_rpmd.py`** now supports two force clients via `--client`: **`lammps`** (default; LAMMPS Python API + `fix ipi` + `fix external` UMA) and **`python`** (`i-pi-py_driver` + `ipi/uma_ipi_driver.py` ASE driver).

### Match the case in the figure (64 H₂O, 32 beads, 1 ps, 243 K, PILE‑G, Δt = 0.1 fs)

From `tests/uma_ice_rpmd`:

**1. One-shot (OpenMM UMA + i‑PI UMA, no TIP4P)** — same physics targets as `run_rpmd_32x32.sh`; it writes **`pipeline_out/rpmd_comparison_32x32.png`** (not `rpmd_comparison.png`):

```bash
cd tests/uma_ice_rpmd
export OPENMMML_UMA_RPMD_CHUNK=4   # if VRAM is tight; optional
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
bash run_rpmd_32x32.sh
```

Then either open `pipeline_out/rpmd_comparison_32x32.png`, or rerun only the plotter:

```bash
python plot_rpmd_comparison.py \
  --openmm pipeline_out/ice_order_openmm_rpmd.csv \
  --ipi pipeline_out/ice_order_ipi_rpmd.csv \
  -o pipeline_out/rpmd_comparison.png \
  --title "RPMD ice order: OpenMM UMA vs i-PI UMA (243 K, PILE-G)" \
  --molecules 64 --beads 32
```

**2. i‑PI leg only** (if you already have `ice_order_openmm_rpmd.csv`):

```bash
cd tests/uma_ice_rpmd
python lammps/build_lammps_ice_data.py --nx 2 --ny 2 --nz 2 -o lammps/data.ice_uma_64
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
python run_ipi_lammps_uma_rpmd.py \
  --molecules 64 --beads 32 --dt-fs 0.1 --steps 10000 \
  --device cuda --ipi-thermostat pile_g
python ipi_order_from_traj.py \
  --traj ipi/ice__i-pi.traj_00.xyz --beads 32 --dt-fs 0.1 \
  -o pipeline_out/ice_order_ipi_rpmd.csv
```

Use **`ice__i-pi.traj_00.xyz`** if your i‑PI install zero-pads indices; `ipi_order_from_traj.py` infers the other bead files from that name (see `tests/uma_ice_rpmd/ipi_order_from_traj.py`).

**3. OpenMM UMA leg** (if you need to regenerate the blue curve with the same settings):

```bash
export OPENMMML_UMA_RPMD_CHUNK=4
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
python run_openmm_rpmd_reference.py \
  --nx 2 --ny 2 --nz 2 --beads 32 --dt 0.1 --steps 10000 \
  --equil 2 \
  --platform cuda --rpmd-thermostat pile-g
```

(`run_rpmd_32x32.sh` and `run_openmm_uma_rpmd_only.sh` default **`RPMD_EQUIL_PS=2`** for dynamics parity with the reference driver; use **`RPMD_EQUIL_PS=0`** only for fast smoke tests.)

### Three curves like the screenshot (OpenMM UMA + TIP4P + i‑PI)

Run **`run_openmm_tip4p_rpmd.py`** with **64 molecules, 32 beads, 0.1 fs, 10000 steps**, write **`pipeline_out/ice_order_tip4p_rpmd.csv`**, then:

```bash
python plot_rpmd_comparison.py \
  --openmm pipeline_out/ice_order_openmm_rpmd.csv \
  --tip4p pipeline_out/ice_order_tip4p_rpmd.csv \
  --ipi pipeline_out/ice_order_ipi_rpmd.csv \
  -o pipeline_out/rpmd_comparison.png
```

### Kinetic temperature panel

**OpenMM** (and **TIP4P**) order CSVs carry **`T_K`** from the simulation. When **`ipi_order_from_traj.py`** is run with thermo (default: merge from `ice__i-pi.md`), **`T_K`** uses i‑PI **`temperature(nm=0)`** (centroid), matching OpenMM’s centroid definition; see [`ipi_thermo_utils.py`](ipi_thermo_utils.py) and [`plot_ipi_thermo.py`](plot_ipi_thermo.py) for extended vs centroid diagnostics.

That’s the main functional gap if you want the bottom panel to look exactly like your current PNG for the i‑PI run.


