# UMA RPMD: centroid kinetic temperature (`T_K`) and high-T artifacts

## What `T_K` is in `ice_order_openmm_rpmd.csv`

`test_uma_ice_rpmd.py` reports **centroid kinetic temperature**: for each atom *i*, centroid velocity  
\(\bar{\mathbf v}_i = \frac{1}{P}\sum_{b=1}^P \mathbf v_i^{(b)}\), then  
\(\mathrm{KE}_{\mathrm{cent}} = \sum_i \frac{1}{2} m_i |\bar{\mathbf v}_i|^2\), and  
\(T = 2\,\mathrm{KE}_{\mathrm{cent}} / (f\,R)\) with \(f = 3\times(\#\text{ particles with }m>0)\) and \(R = 8.314\times10^{-3}\) kJ/(mol·K).  
This matches the degree-of-freedom counting used in OpenMM’s PILE-G centroid thermostat (`ReferenceRpmdKernels::applyBussiCentroidThermostat`).

It is **not** the ring-polymer “quantum” kinetic estimator used in some PIMD analyses; it is the **classical** centroid kinetic temperature used consistently with PILE-G.

## PBC wrapping and UMA RPMD

**Force evaluation (dynamics):** [`test_uma_ice_rpmd.py`](test_uma_ice_rpmd.py) defaults to **`use_atom_wrap_for_lammps_parity=True`** (per-atom ASE-style wrap in the batched UMA path), matching LAMMPS `fix external` / Fairchem. Use **`--use-molecular-wrap-for-uma`** for molecular water translation (`_wrap_water_molecules_per_bead`) instead—for example when per-atom wrapping of unwrapped RPMD coordinates produces unrealistic O–H separations in the NN graph and you observe runaway centroid `T_K`. See [`docs/LAMMPS_OPENMM_UMA_PARITY.md`](docs/LAMMPS_OPENMM_UMA_PARITY.md).

**Order-parameter postprocessing:** Q6 / q_tet CSV rows apply **`wrap_cartesian_orthorhombic`** to each bead’s Cartesian positions before `ice_order_metrics_path_integral`, consistent with [`ipi_order_from_traj.py`](ipi_order_from_traj.py) (default per-atom wrap).

## Why you might see `T_K` ≫ 243 K (e.g. &gt; 1000 K)

1. **Wrap / NN input geometry** — if **per-atom** force wrapping interacts badly with your coordinates, try **`--use-molecular-wrap-for-uma`** as a diagnostic, or reduce timestep / equilibrate.

2. **Real drift / heating**  
   If centroid KE grows in time (e.g. within the first picoseconds), the reported `T_K` will grow with it. That usually points to **integrator + forces + thermostat** balance, not the CSV writer:
   - **Timestep** too large for UMA + many RPMD beads (try smaller `--dt` in `run_openmm_rpmd_reference.py` / `test_uma_ice_rpmd.py`).
   - **No equilibration** (`--equil 0` in the reference driver): starting production immediately after minimization can give a short transient; for serious runs, use a short NVT equilibration (e.g. 0.5–2 ps).
   - **Thermostat coupling**: increase `--rpmd-centroid-friction` (Bussi on centroid) and/or `--rpmd-friction` (PILE internal modes) if the bath is too weak for your settings.
   - **ML potential / batched path**: inconsistent energies or forces (see below) can drive unstable dynamics.

3. **Check potential energy scale**  
   For a small ice box (~64 H₂O), **total** potential energy in OpenMM should be on the order of **10³–10⁵ kJ/mol** (extensive but not ~10⁷). If `PE_kj_mol` in the order CSV looks like **−10⁷**, treat energy/force evaluation as **suspect** (model output interpretation, unit conversion, or batched vs single-path parity) and debug UMA/openmmml before trusting `T_K` or structural observables.

4. **After code changes**  
   Centroid `T_K` uses **`ndof = 3 × (# particles with mass > 0)`**. If your topology ever includes massless sites, older code that used `3 × n_atoms` could bias `T` slightly; that correction does **not** by itself explain runaway heating to thousands of kelvin.

## Practical checklist

- [ ] Confirm **`use_atom_wrap_for_lammps_parity`** choice (default **True** for LAMMPS parity; `--use-molecular-wrap-for-uma` for molecular).  
- [ ] Reduce `dt` (e.g. 0.05 fs) and/or beads if instability persists.  
- [ ] Add short `--equil` production in `test_uma_ice_rpmd.py` via `run_openmm_rpmd_reference.py` (default is 2 ps; use smaller values only for smoke tests).  
- [ ] Compare early `T_K` trend with centroid friction and internal friction.  
- [ ] Validate **PE** magnitude and force norms against a trusted single-bead or smaller-bead run.
