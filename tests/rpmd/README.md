# UMA/RPMD Test Suite

Validates UMA potential integration with RPMD in OpenMM: force accumulation, hybrid RPMD, and numerical stability.

## Running tests

```bash
cd tests/rpmd
python run_all_uma_rpmd_tests.py --quick   # ~5 min
python run_all_uma_rpmd_tests.py           # Full suite ~15-20 min
python run_all_uma_rpmd_tests.py --category regression   # After kernel changes
python run_all_uma_rpmd_tests.py --category forces       # After force changes
python run_all_uma_rpmd_tests.py --category hybrid      # After hybrid RPMD changes
```

| File | Purpose | ~Runtime |
|------|---------|----------|
| test_bug_regression.py | Force overwrite, classical centroid regression | ~30s |
| test_uma_force_accumulation.py | UMA + bond force accumulation | ~20s |
| test_classical_centroid_force.py | Classical particle centroid averaging | ~30s |
| test_uma_hybrid_rpmd_cavity.py | UMA + hybrid RPMD + cavity | ~40s |
| test_batch_consistency.py | Batch vs sequential consistency | ~60s |
| test_uma_rpmd_stress.py | P=32 beads, multi-group, 500 steps | ~90s |
| test_quantum_statistics.py | Quantum harmonic oscillator vs theory | ~120s |
| test_edge_cases.py | Single particle, extremes, exclusion, NaN/Inf | ~80s |
| test_end_to_end_pipeline.py | Full pipeline validation | ~50s |

## When to run

- **After kernel changes (C++/GLSL):** `test_bug_regression.py`, `test_uma_force_accumulation.py`, `test_classical_centroid_force.py`
- **Before commit:** `run_all_uma_rpmd_tests.py --quick`
- **Before release:** Full suite

## Bug coverage

| Bug | Primary tests |
|-----|---------------|
| Force overwrite (batched path) | test_uma_force_accumulation, test_bug_regression |
| Classical centroid wrong (bead 0 vs centroid) | test_classical_centroid_force, test_bug_regression |
| Batch/sequential inconsistency | test_batch_consistency |
| Buffer overflow (many beads) | test_uma_rpmd_stress |
| Particle exclusion | test_uma_hybrid_rpmd_cavity, test_edge_cases |
| Quantum statistics wrong | test_quantum_statistics |

## Common failures

- **"Force accumulation error too large"** → Check force accumulation (+= vs =) in kernels
- **"Velocity not aligned with centroid force"** → Check syncClassicalBeads, integrateStepHybrid
- **"Cavity has non-zero UMA force"** → Check atom_indices in PythonForce
- **"Energy became NaN or Inf"** → Check force computation, timestep, precision

## CI

```bash
cd tests/rpmd && python run_all_uma_rpmd_tests.py --quick
```

Exit 0 = pass, non-zero = failures.

## Limits

- Hybrid RPMD full support: Common (GPU) only; Reference may skip hybrid tests.
- UMA models required; tests skip if unavailable. Use `uma-s-1p1-pythonforce-batch` for fast runs.
