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

All performance alternatives default to OFF. The GPU-only minimizer mode defaults to ON. Q32.32 force accumulation retains the OpenCL two-word algorithm; there is no native 64-bit atomic add. Compilation or a skipped GPU test is not a performance or numerical-parity result.
