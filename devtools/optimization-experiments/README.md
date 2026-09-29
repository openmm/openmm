# CUDA and common-compute optimization experiments

This contribution exposes the complete implemented optimization suite from a local OpenMM 8.2 experiment as a current-master source port. It targets CPU launch preparation, PME, bonded forces, constrained integration, neighbor-list preparation, reordering, and coordinate transfers. It is submitted for review; the validation listed below remains required before merging or default enablement. No performance improvement or scientific equivalence is established for the current-master port.

The starting upstream commit is `0c4bcaba734ea574f52f5e09ffaaf4d6fe6109c2`. The earlier CUDA argument-pointer cache commit remains in the branch with its regression test. Unlike the historical local build, all newly added experimental paths require explicit environment opt-in and their runtime eligibility checks. The argument-pointer cache is the already-contributed unconditional host optimization. Existing upstream optimizations stay enabled as upstream defines them.

## Complete scope and disposition

- [Common, bonded, coordinate-transfer and constrained-integration accounting](common-and-integration.md).
- [PME, sorting, energy-reduction and CUDA-context accounting](pme-and-cuda.md).
- [Reorder, neighbor-list and direct-launch accounting](reorder-and-neighbors.md).
- [Historical requested options](legacy-requested-options.json): 37 performance requests plus five disabled trace options. This is provenance, not a recommended configuration for this branch; some options no longer exist.
- Three initially ungated optimizations now have explicit switches: `OPENMM_EXPERIMENT_BAROSTAT_POTENTIAL_ONLY`, `OPENMM_EXPERIMENT_PME_ENERGY_ONLY_SKIP_FORCE`, and `OPENMM_EXPERIMENT_PME_REAL_GRID_CLEAR`.

Boolean switches accept exactly `1`. SETTLE block selection accepts `128`; eligible direct-force geometry accepts `5` or `6`; the guarded neighbor-skin policy accepts the values documented in its accounting table. Set options before creating a Context. Flags express requests, and do not override platform, precision, force, size, or topology restrictions. Most experiments target single-context mixed-precision CUDA PME with the exact Langevin Middle integrator. Experiments are internal development switches, not proposed stable user API.

Three old requests are intentionally not duplicated in executable current-master code:

| Historical request | Current-master disposition |
| --- | --- |
| `OPENMM_EXPERIMENT_PME_SORT_CADENCE` | Equivalent large-system cadence already exists upstream ([#5305](https://github.com/openmm/openmm/pull/5305)). |
| `OPENMM_EXPERIMENT_PARALLEL_BAROSTAT` | Upstream already has parallel large-molecule scaling and newer component support ([#5432](https://github.com/openmm/openmm/pull/5432)). |
| `OPENMM_EXPERIMENT_REUSE_BOX_POSQ_DOWNLOAD` | The newer packed-coordinate setter eliminated that download entirely ([#4945](https://github.com/openmm/openmm/pull/4945)). The fused/parallel setter experiments were adapted to this layout. |

The old `OPENMM_EXPERIMENT_CUDA_ARGUMENT_POINTER_CACHE` switch is also unnecessary: the independently contributed cache is unconditional. Trace switches remain optional diagnostics and do not activate their associated optimization.

## Historical full patch

[`legacy-openmm82.patch`](legacy-openmm82.patch) preserves all 33 changed source files against public tag `8.2.0`, commit `53770948682c40bd460b39830d4e0f0fd3a4b868`. It is review/provenance material, **not a patch to apply to this branch**, and is not used by the build. It includes superseded paths and the historical default-on behavior; those should not be mistaken for the current implementation.

[`legacy-manifest.json`](legacy-manifest.json) binds each original and experimental file by SHA-256, plus LF-normalized comparison hashes. Applying the patch to the pristine 8.2 source and comparing all resulting source texts to the cumulative snapshot passed. No Folding work-unit files, account configuration, private input data, or binaries are included.

## Validation

The full port was configured and built using the repository's CMake project on Windows x64, Visual Studio 2022 / MSVC 19.44, Release, with CUDA 11.8 headers/import libraries. `OpenMM`, `OpenMMCUDA`, and `TestCudaKernel` built successfully. Existing upstream linker/conversion/encoding warnings remain. Building a CUDA plugin does not compile every runtime-generated device kernel or run its tests.

Fresh CPU checks on the actual current-port sources passed:

- Hilbert specialization: all 16,777,216 3D/8-bit coordinates equal the bundled generic implementation; 1,267 outside-domain and 1,331 disabled-path cases also pass.
- Neighbor-skin selection and floating-point padding: 598 checks pass.
- Packed-coordinate setter: 314 comparisons, including 153 fused and 161 fallback cases (six failure injections are included). Full position/correction bytes, charges, padding, upload/copy/reorder ordering, rounding/FTZ/DAZ controls and fallback selection match the actual baseline function.

The setter uses deterministic CPU mock workers and a CPU interpretation of the actual xyz-copy kernel. It does not validate native ThreadPool concurrency, CUDA transfers, stream lifetime, or GPU arithmetic. The prior argument-cache host mock tests remain applicable to its unchanged implementation. GPU simulation tests and benchmarks were **not run for this port**, and the local Folding installation was not changed by preparing this contribution.

Representative runtime device-source compilation also passed without loading the CUDA driver: 13 common/bonded/integration variants and the PME/sort/reduction variants recorded in [validation.json](validation.json), using NVRTC 11.8 and `compute_90` PTX. The mixed-precision envelope was reconstructed from current source; this is not a capture of a live Context, driver JIT or GPU execution. Four representative bonded host-generation cases also passed.

Reproduce the self-contained CPU checks from a Git checkout containing the baseline commit (no OpenMM or CUDA library loading):

```sh
cmake -S devtools/optimization-experiments/cpu-checks -B build-experiment-cpu
cmake --build build-experiment-cpu --config Release
ctest --test-dir build-experiment-cpu -C Release --output-on-failure
```

The setter check is Windows/MSVC x64-specific and requires Python 3 and Git to extract actual source bodies. Hilbert and skin checks have no CUDA dependency. Changes to source layout may require updating the extraction script; extraction failures are not silently ignored.

## Required before merge or performance claims

Run force/energy/derivative comparisons, constrained trajectories and RNG/checkpoint continuation; exercise reorders, box changes, parameter updates, sort fallback/flag transitions, multi-context and unsupported-force/precision/platform cases. Check default-OFF behavior with the normal OpenMM suite and other common-compute backends. Validate all relevant option combinations on GPU, particularly compiler-sensitive bonded fence removal, fused SETTLE, floating-point PME spreading/permutations, and the smaller real FFT grid. Benchmark matched workloads with each candidate isolated and then in combination. Lower launch counts, smaller allocations, or high GPU utilization alone are not performance evidence.

## AI assistance and provenance

AI assistance (OpenAI Codex) was used for implementation, porting, review, checks and this description. This contribution makes that explicit under the repository's AI policy and does not assert that CPU checks replace scientific validation. Existing OpenMM and bundled Rice University license notices are retained. The large suite is presented together for review of its complete scope; it can be split into focused changes before any merge decision.
