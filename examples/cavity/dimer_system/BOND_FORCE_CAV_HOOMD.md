# Bond force in cav-hoomd: where it comes from

This note traces where the bond force in cav-hoomd is defined and documents the **root cause** of the OpenMM–HOOMD bond force mismatch and its fix.

---

## Root cause and fix (position wrapping)

**Harmonic bond is F = k(r0 − r)r̂.** Same k, r0, positions, box ⇒ same F.

**Root cause:** The function `make_positions_minimage_consistent()` in `compare_openmm_cav_hoomd_forces.py` wrapped all particles relative to particle 0 so that pairwise distances matched minimum image. That **broke molecular bonds**: bonded atoms could end up on opposite sides of the periodic box. **OpenMM's `HarmonicBondForce` uses direct Euclidean distance, not minimum image.** So when bonded atoms were ~39 Bohr apart (across the box), OpenMM computed a huge restoring force even though the minimum-image distance was ~2 Bohr.

**Fix:** We use `unwrap_molecules(cfg)` instead. It keeps each bond's two atoms close by moving the second atom into the same image as the first, then applies a single global translation so the minimum coordinate is in [0, box). OpenMM and the compare script build from these unwrapped positions; bond forces then match the by-hand formula (and typically agree with HOOMD to ~1e-2 Ha/Bohr when HOOMD runs on the original GSD). Full agreement below 1e-5 would require HOOMD to run on the same unwrapped, in-box configuration.

**Verification:** After the fix, the single-bond diagnostic shows by-hand and OpenMM matching for the same bond (e.g. bond 147: |F| ~ 0.036 Ha/Bohr for both). Remaining bond max error ~1e-2 Ha/Bohr when comparing OpenMM (unwrapped) to HOOMD (original GSD) is from comparing different configs. Coulomb can differ due to PPPM vs Ewald; LJ and Cavity are at ~1e-4 to 1e-3 Ha/Bohr.

---

## Step back: the two real puzzles

Assume HOOMD is fine — any scale difference is at most a prefactor (e.g. 0.5). The puzzling things are:

### 1. Why is “within tolerance” bond error ~1e-3 when LJ is ~1e-4–1e-5?

We report *absolute* errors (|ΔF| in Ha/Bohr). Bond forces are **much larger** than LJ (|F_bond| ~ tens of Ha/Bohr, |F_LJ| ~ 0.01–0.1). So the same *relative* error in both codes would produce a much larger *absolute* bond error. A typical bond error of ~1e-3 with |F_bond| ~ 30 would be a relative error ~3e-5; LJ error ~3e-5 with |F_LJ| ~ 0.05 would be relative ~6e-4. So it’s not obvious that bond is “worse” — we need **relative** error (|ΔF|/|F|) or at least typical |F| per component to compare apples to apples.

**What the numbers show:** Once you add p50(|F|) and relative error (|ΔF|/|F|) per component, you see that **bond and Coulomb have large relative errors** (p50 rel_err ~0.25–0.5) while **LJ has small relative error** (p50 rel_err ~0.02). So the “crazy” 1e-3 bond vs 1e-5 LJ is **not** just magnitude — bond and Coulomb are *relatively* wrong (~25–50% typical); LJ agrees to ~2%. The puzzle becomes: why do bond and Coulomb disagree in relative terms while LJ agrees?

### 2. Why are there a few huge outliers?

A few particles have bond error ~50 Ha/Bohr (same order as |F_bond| itself); most are ~1e-3. That suggests for those particles we’re comparing the right OpenMM force to something that is either (a) the force on the *other* atom of the same bond (i.e. −F, so error ~ 2|F|), or (b) the force from a *different* bond, or (c) zero/wrong row. The fact that worst particles are **consecutive pairs** (294,295), (483,482) is consistent with (a): if for some bonds we swapped which particle we’re reading, we’d see F vs −F on that pair → huge error on both. So the outlier pattern points at **bond-pair assignment** (which particle is “first” in the bond, or which row maps to which tag) for a subset of bonds, not a global scale/unit bug.

**Repos (GitRepos):**

- **cav-hoomd:** `~/GitRepos/cav-hoomd` (or `cav-hoomd` under your clone root)
- **HOOMD-blue:** `~/GitRepos/hoomd-blue`

## cav-hoomd does not implement bond forces in C++

- **cav-hoomd source** (`~/GitRepos/cav-hoomd/src` or similar) contains:
  - `CavityForceCompute.{cc,h}` — cavity–molecule coupling
  - `DipoleResponseForceCompute`, `PerturbationForceCompute` — other physics
  - **No** bond/Harmonic force implementation

- Bond forces come from **HOOMD-blue** via `hoomd.md.bond.Harmonic`.

## Where cav-hoomd sets up bonds

| File | Location | What it does |
|------|----------|--------------|
| `cavitymd/simulation/core.py` | `setup_force_parameters()` ~L2202 | `harmonic = hoomd.md.bond.Harmonic()`, then sets params and appends to forces |
| `cavitymd/setup.py` | `create_harmonic_bonds()` ~L358 | Same: `harmonic.params['O-O']`, `harmonic.params['N-N']` |

Relevant snippet from `core.py`:

```python
# Setup harmonic bonds
harmonic = hoomd.md.bond.Harmonic()
harmonic.params['O-O'] = dict(k=2*0.36602, r0=2.281655158)
harmonic.params['N-N'] = dict(k=2*0.71625, r0=2.0743522177)
forces.append(harmonic)
```

- **O-O:** `k = 0.73204` Ha/Bohr², `r0 = 2.281655158` Bohr  
- **N-N:** `k = 1.4325` Ha/Bohr², `r0 = 2.0743522177` Bohr  

These match OpenMM `diamer_forcefield.xml` (same potential, Ha/Bohr² and Bohr). cav-hoomd uses Hartree/a.u. for this setup.

## Where the per-particle bond force array comes from

In `debug_force_components.py` we do:

```python
with force_obj.cpu_local_force_arrays as fa:
    arr = np.array(fa.force, copy=True)[:, :3]
```

For the force we match by name `"harmonic"` or `"bond"`, that `force_obj` is the same `hoomd.md.bond.Harmonic()` instance. So:

- **`cpu_local_force_arrays`** is provided by HOOMD-blue’s MD backend for that Harmonic bond force.
- The layout (local/tag/bond order) and units of `fa.force` are defined entirely in **HOOMD-blue**, not in cav-hoomd.

The bond force layout in HOOMD-blue (local index, same as net_force) is correct; the large OpenMM–HOOMD bond mismatch was caused by **position wrapping** (see "Root cause and fix" above), not by HOOMD's bond implementation.

## What to inspect next (HOOMD-blue)

To understand layout and scale of the bond force array:

1. **HOOMD-blue source**  
   - `hoomd.md.bond.Harmonic` is a Python wrapper; the real work is in the C++ bond force compute (often something like `BondForceCompute` or `HarmonicBondForceCompute`).
   - Find where that compute writes into the per-particle force array and what index convention it uses (e.g. local particle index, bond index, etc.).

2. **Docs / API**  
   - HOOMD’s “ParticleLocalAccessBase” and “Accessing System Configurations With MPI” describe local vs global indexing and `rtag`.
   - Check whether bond forces use the same local layout as `net_force` or a different one (e.g. bond list order).

3. **Installed HOOMD**  
   - If HOOMD is from conda/pip, search the installed package for bond-related sources (e.g. under `site-packages/hoomd` or the compiled extension).  
   - Alternatively, use a HOOMD-blue repo/source build and search for “bond” and “force” in the MD and bond modules.

## Summary

| Component | Implemented in | Notes |
|----------|----------------|--------|
| Bond potential (k, r0) | cav-hoomd Python | `core.py` / `setup.py`; matches OpenMM |
| Bond force *evaluation* | HOOMD-blue | `hoomd.md.bond.Harmonic` → HOOMD backend |
| Per-particle bond array we read | HOOMD-blue | `cpu_local_force_arrays` on that Harmonic force |

So “going deeper” on the bond force means debugging **HOOMD-blue’s** bond force compute (buffer layout and units), not cav-hoomd’s own source.

---

## HOOMD-blue source (e.g. `~/GitRepos/hoomd-blue`)

If HOOMD-blue is checked out under `~/GitRepos/hoomd-blue`, the bond force chain is:

| Path | Role |
|------|------|
| `hoomd/md/bond.py` | Python: `Harmonic` class, `_cpp_class_name = "PotentialBondHarmonic"` |
| `hoomd/md/PotentialBond.h` | C++: `PotentialBond<EvaluatorBondHarmonic, BondData>`, `computeForces()` |
| `hoomd/md/EvaluatorBondHarmonic.h` | C++: `force_divr = K*(r_0/r - 1)`, same as standard harmonic |

### Layout: bond force is written in **local index** order

In `PotentialBond.h::computeForces()`:

- Bond list is iterated; for each bond, `bond.tag[0]` and `bond.tag[1]` are the particle tags.
- **Local indices:** `idx_a = h_rtag.data[bond.tag[0]]`, `idx_b = h_rtag.data[bond.tag[1]]` (rtag = tag → local index).
- Forces are written to **`h_force.data[idx_a]`** and **`h_force.data[idx_b]`** (local index).

So the bond force array is in **local particle index** order, same as `net_force`. Indexing it with `arr[rtag]` gives force per **tag**, so our pipeline (arr[rtag] → tag order, then F[gsd_idx] → OpenMM order) is consistent with the C++ implementation.

### Units in the evaluator

`EvaluatorBondHarmonic` uses `K` and `r_0` from params; the force written is `force_divr * dx` (in the same length/energy units as `K` and `r_0`). cav-hoomd passes k in Ha/Bohr² and r0 in Bohr, so the C++ output is in Ha/Bohr. No extra conversion is done in the evaluator.

### Implications (root cause was position wrapping)

From the HOOMD-blue source:

1. **Layout:** Bond forces are stored by local index; our use of `rtag` to get tag order is correct.
2. **Formula:** Harmonic force matches the usual \(K(r_0/r - 1)\hat{r}\), so no factor-of-2 or sign mistake in the formula.
3. **Units:** Params are passed through; cav-hoomd uses Ha/Bohr² and Bohr, so forces should be Ha/Bohr.

The remaining OpenMM–HOOMD bond discrepancy was due to **position wrapping** (see "Root cause and fix" above). After switching to `unwrap_molecules()` for OpenMM, by-hand and OpenMM match; the residual ~1e-2 Ha/Bohr when comparing to HOOMD on the original GSD is from comparing different configs (OpenMM unwrapped vs HOOMD on raw GSD).
