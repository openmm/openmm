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
| `OPENMM_METAL_FAST_MINIMIZE_SHUFFLE` | Common local-memory reduction | Common CUDA/HIP shuffle reduction in the minimizer only |
| `OPENMM_METAL_FAST_CONSTANT_POTENTIAL_REDUCTION` | Base ConstantPotential local-memory reduction | Common shuffle reduction in the base ConstantPotential kernels only |
| `OPENMM_METAL_FAST_PME_FLOAT_SPREAD` | Fixed-point PME grid | Common floating-grid PME spreading, without the conversion dispatch |
| `OPENMM_METAL_FAST_LJPME_FLOAT_SPREAD` | Fixed-point grids in LJPME systems | Floating electrostatic and dispersion grids in LJPME systems |
| `OPENMM_METAL_FAST_CONSTANT_POTENTIAL_FLOAT_SPREAD` | Fixed-point ConstantPotential PME grid | Common floating-grid spreading |

All performance alternatives default to OFF. The GPU-only minimizer mode defaults to ON. Q32.32 force accumulation retains the OpenCL two-word algorithm; there is no native 64-bit atomic add. Compilation or a skipped GPU test is not a performance or numerical-parity result.
