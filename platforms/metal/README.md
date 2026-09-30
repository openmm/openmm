# Experimental Metal Platform

Metal 3 host API and runtime-compiled MSL 3.0; one Apple-silicon GPU and single precision.

This series introduces the OpenCL/Common baseline in component-sized changes before adding independently selectable alternatives. New component sources are wired into the library by the Common-runtime integration change. Optional force plugins, multi-GPU and CPU fallbacks are out of scope.

## Components in this revision

- Common/OpenCL source adaptation and MSL primitives.
- Reflected argument buffers and queue-ordered dispatch.
- OpenCL-derived neighbor lists, sorting and integration utilities.
- Bundled VkFFT with a private C++17/metal-cpp bridge.
- Scoped floating accumulator variants for GPU-only minimization.
- Core force/integrator/state Platform registration; no optional force plugins.

## Independent comparison switches

| Build option | OFF | ON |
| --- | --- | --- |
| `OPENMM_METAL_RECORD_AND_COMMIT` | Batch up to 64 operations or a synchronization boundary | Commit each operation (default) |
| `OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS` | Q32.32 minimization with limited range diagnostics; no CPU fallback | Scoped GPU floating accumulators for large-force minimization (default) |
| `OPENMM_METAL_NATIVE_FLOAT_ATOMICS` | OpenCL-style float-add CAS loop | Native Metal float atomic add, independently of accumulator representation |
| `OPENMM_METAL_FAST_MATH` | Strict Metal compilation | Broad fast math where compatible with numerical safety requirements |
| `OPENMM_METAL_FAST_SHORT_LIST_SORT` | OpenCL sorting selection | CUDA short-list selection using the shared alternative kernel |
| `OPENMM_METAL_FAST_MINIMIZE_SHUFFLE` | Common local-memory reduction | Common CUDA/HIP shuffle reduction in the minimizer only |
| `OPENMM_METAL_FAST_CONSTANT_POTENTIAL_REDUCTION` | Base ConstantPotential local-memory reduction | Common shuffle reduction in the base ConstantPotential kernels only |
| `OPENMM_METAL_FAST_CONSTANT_POTENTIAL_CG_REDUCTION` | CG local-memory reductions | Common scalar and compensated `float2` shuffle reductions in CG only |
| `OPENMM_METAL_FAST_CONSTANT_POTENTIAL_MATRIX_REDUCTION` | Matrix-solver local-memory reduction | Common shuffle reduction in the matrix solver only |
| `OPENMM_METAL_FAST_CONSTANT_POTENTIAL_MATRIX_BROADCAST` | Matrix-solver local-memory broadcasts | Common lane-shuffle broadcasts in the matrix solver only |
| `OPENMM_METAL_FAST_CUSTOM_NONBONDED_GROUPS_SHUFFLE` | Interaction-group local-memory maximum | Common CUDA xor-shuffle maximum for CustomNonbonded interaction groups |
| `OPENMM_METAL_FAST_CUSTOM_MANY_PARTICLE_BALLOT` | CustomManyParticle local neighbor-block flags | Common CUDA ballot-based neighbor-block scan |
| `OPENMM_METAL_FAST_LCPO_BALLOT` | LCPO local neighbor-block flags | Common CUDA ballot-based neighbor-block scan |
| `OPENMM_METAL_FAST_PME_FLOAT_SPREAD` | Fixed-point PME grid | Common floating-grid PME spreading, without the conversion dispatch |
| `OPENMM_METAL_FAST_LJPME_FLOAT_SPREAD` | Fixed-point grids in LJPME systems | Floating electrostatic and dispersion grids in LJPME systems |
| `OPENMM_METAL_FAST_CONSTANT_POTENTIAL_FLOAT_SPREAD` | Fixed-point ConstantPotential PME grid | Common floating-grid spreading |
| `OPENMM_METAL_FAST_BLOCK_BOUNDS` | Serial atom-block bounds | SIMD-cooperative bounds and radius reduction |
| `OPENMM_METAL_FAST_FP16_BOUNDS` | Float sorted/large bounding boxes | Conservatively rounded half4 storage; public Common bounds remain float |
| `OPENMM_METAL_FAST_NEIGHBOR_BALLOT` | OpenCL local flags and atom prefix sums | SIMD ballot/popcount block iteration and atom compaction |
| `OPENMM_METAL_FAST_SPARSE_PAIRS` | All interactions in tiles | Separate sparse pairs and the CUDA sparse-pair force loop |
| `OPENMM_METAL_FAST_NONBONDED_SHUFFLE` | OpenCL local-memory force tiles | CUDA register/shuffle force template with Metal spelling adaptations |
| `OPENMM_METAL_FAST_CUSTOM_GB_VALUE_SHUFFLE` | Common CustomGB value local arrays | Register exchange for value, parameters, and secondary accumulators |
| `OPENMM_METAL_FAST_CUSTOM_GB_ENERGY_SHUFFLE` | Common CustomGB energy local arrays | Register exchange for force and parameter-derivative state |
| `OPENMM_METAL_FAST_GBSA_BORN_SHUFFLE` | Common Born-sum local structs | Register exchange of Born data and secondary sums |
| `OPENMM_METAL_FAST_GBSA_FORCE_SHUFFLE` | Common GBSA force local structs | Register exchange of force/Born state |
| `OPENMM_METAL_FAST_DPD_PARTICLE_SHUFFLE` | Common DPD local particle arrays | Register broadcasts preserving original pair/RNG visitation order |
| `OPENMM_METAL_FAST_DPD_TILE_BROADCAST` | Local-memory tile-counter broadcast | SIMD lane-zero broadcast |
| `OPENMM_METAL_FAST_CUSTOM_HBOND_SHUFFLE` | Common local acceptor structs | Register rotation of acceptor positions and forces |
| `OPENMM_METAL_FAST_CENTROID_REDUCTION` | Common centroid local reduction | SIMD reductions with per-group shared partials |
| `OPENMM_METAL_FAST_RG_REDUCTION` | Common radius-of-gyration local reduction | SIMD reductions with per-group shared partials |
| `OPENMM_METAL_FAST_RMSD_REDUCTION` | Common RMSD local reduction | SIMD center and correlation reductions |
| `OPENMM_METAL_FAST_ORIENTATION_REDUCTION` | Common orientation local reduction | SIMD center and correlation reductions |

All performance alternatives default to OFF. The GPU-only minimizer mode defaults to ON. Q32.32 force accumulation retains the OpenCL two-word algorithm; there is no native 64-bit atomic add. Compilation or a skipped GPU test is not a performance or numerical-parity result.
