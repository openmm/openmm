# OpenMM / cav-hoomd parity run

This document describes how to run the **exact same** simulation in cav-hoomd as the OpenMM reference, then compare F(k,t) output to interpret relaxation discrepancies.

## OpenMM reference command

```bash
python run_simulation.py --dimers 250 --g 0.0 --temp 100.0 --dt 0.001 --equil 200 --prod 3000 \
  --cavity-freq 1555 --report-interval 10000 --pdb-interval 10000 --enable-fkt \
  --fkt-output-period-ps 1.0 --fkt-ref-interval-ps 100.0 --constant-density --no-dipole --pdb ./test.pdb
```

Parameters: 250 dimers, g=0 (no cavity), 100 K, 1 fs, Bussi (τ 5 ps), constant density (40 Bohr box), 200 ps equil + 3000 ps prod (3200 ps total), F(k,t) with k=113.4 nm⁻¹, 50 wavevectors, ref interval 100 ps, max refs 10, output period 1 ps.

## Cav-hoomd workflow (two steps)

**Requirement:** cav-hoomd (HOOMD-blue with cavitymd plugin) must be installed in the environment where you run the cav-hoomd script.

### Step 1: Create initial GSD (250 dimers, constant density)

From `examples/cavity/dimer_system/`:

```bash
python initlattice_equilibrium.py --job-dir . --replica 0 --nmol 250 --temperature 100 --seed 42
```

This creates `molecular-0.gsd` (500 particles, box 40 Bohr, 80% O–O / 20% N–N). No cavity particle for g=0 parity.

### Step 2: Run cav-hoomd simulation (exact parity)

From `examples/cavity/dimer_system/`:

```bash
python run_cav_hoomd_advanced.py --no-cavity --temperature 100 --runtime 3200 \
  --fixed-timestep --timestep 1.0 --molecular-bath bussi --molecular-tau 5.0 \
  --enable-fkt --fkt-kmag 113.4 --fkt-wavevectors 50 --fkt-ref-interval 100.0 \
  --fkt-max-refs 10 --fkt-output-period-ps 1.0 --input-gsd molecular-0.gsd \
  --device GPU --seed 42 --console-output-period-ps 10
```

Cav-hoomd writes F(k,t) files under `no_cavity/`. Compare those F(k,t) curves (and relaxation time) to OpenMM’s `fkt_lambda0.0000_fkt_ref_*.txt` to see whether dynamics or the F(k,t) implementation explains the difference.

## Known accepted differences

### Long-range electrostatics: PME (OpenMM) vs PPPM (cav-hoomd)

OpenMM uses Particle Mesh Ewald (PME) for long-range Coulomb interactions, while
cav-hoomd uses Particle–Particle Particle–Mesh (PPPM, resolution=[32,32,32],
order=6, alpha=0.0). Both approximate the Ewald sum but use different
interpolation schemes (B-spline vs cardinal B-spline). This produces a per-
particle force error of approximately 0.25% relative to the exact Ewald result.

This difference is inherent to the two codes: OpenMM does not offer PPPM, and
HOOMD does not offer PME. The error is small enough for production use and does
not affect equilibrium thermodynamic averages.

**Tuning guidance for minimal discrepancy:**
- OpenMM PME error tolerance: `ewaldErrorTolerance=0.0005` (current)
- cav-hoomd PPPM: `resolution=[32,32,32]`, `order=6`
- Increasing OpenMM's PME grid density (via `pmeGridSpacing`) or lowering the
  error tolerance will not eliminate the algorithmic difference, only reduce
  OpenMM's internal error.

### RNG streams

HOOMD uses a counter-based Philox RNG seeded by `(RNGIdentifier, timestep, seed)`.
OpenMM uses `SimTKOpenMMUtilities`. Even with identical seeds, individual
trajectories will diverge stochastically. Statistical averages (temperature,
F(k,t), etc.) should converge to the same values within sampling error.

## Parameter mapping (OpenMM → cav-hoomd)

| OpenMM                                   | cav-hoomd                                          |
| ---------------------------------------- | --------------------------------------------------- |
| 250 dimers, constant density             | GSD from initlattice with `--nmol 250` (box 40 Bohr) |
| g=0                                      | `--no-cavity`                                       |
| temp 100 K                               | `--temperature 100`                                 |
| dt 0.001 ps = 1 fs                       | `--fixed-timestep --timestep 1.0`                    |
| equil 200 + prod 3000 = 3200 ps          | `--runtime 3200`                                    |
| Bussi τ 5 ps                             | `--molecular-tau 5.0`                               |
| F(k,t): k=113.4, 50 wavevectors          | `--fkt-kmag 113.4 --fkt-wavevectors 50`             |
| F(k,t): ref interval 100 ps, max refs 10 | `--fkt-ref-interval 100.0 --fkt-max-refs 10`         |
| F(k,t): output period 1 ps               | `--fkt-output-period-ps 1.0`                        |
| seed 42                                  | `--seed 42`                                         |
| report every 10000 steps ≈ 10 ps         | `--console-output-period-ps 10`                     |
