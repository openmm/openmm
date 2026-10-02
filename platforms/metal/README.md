# Experimental Metal Platform

- Metal 3 host API and runtime-compiled MSL 3.0.
- One Apple-silicon GPU; single precision only.
- Existing Common Compute host classes and kernel algorithms.
- OpenCL-derived neighbor-list, sorting, and constraint orchestration.
- Bundled VkFFT for GPU FFT/PME.
- Independent build-time switches for submission policy, GPU-only large-force
  minimization, native float atomics, fast math, and each implemented fast path.

This is an experimental core-platform port, not an assertion of numerical or
performance parity. It registers the core force, integrator, constraint, virtual
site, state/checkpoint, and minimization interfaces. Optional force plugins,
multiple GPUs, mixed/double precision, and CPU PME are outside this version.
Unsupported device/precision properties are rejected, not silently replaced
with a CPU implementation. User-defined CustomCPP/Python callbacks retain their
normal Common semantics; they are not fallback algorithms for built-in forces.

## Reuse and Metal-specific code

Common host classes and kernel sources are compiled directly from
`platforms/common`. OpenCL kernel sources are embedded directly from
`platforms/opencl/src/kernels`; they are not copied into independent Metal
force implementations. Narrow Common hooks allow disabling CPU minimizer
recovery, selecting its existing overflow-safe reductions, and decoding a
backend's accumulator buffer. Existing platforms retain their original defaults.

The backend's small source adapter handles entry-point arguments, address
spaces, vector constructors, and explicit execution builtins. It preserves
preprocessor branches and shared mathematical bodies; it is not a general
OpenCL compiler. Reflected argument buffers avoid the limit of 31 direct Metal
buffer slots. Each dispatch retains its own argument snapshot and declares its
indirect resources. Native test MSL still uses direct buffer slots.

Forces use the existing 64-bit component-plane layout. The OpenCL two-word
atomic accumulation algorithm supplies fixed-point addition on Metal; ordinary
64-bit integers do not imply native 64-bit atomic addition. Shared reductions
are read only after their accumulation dispatches finish. A small MSL preamble
provides the required math and synchronization operations.

With `OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS=ON` (the default), minimization
temporarily uses floating accumulators in the original Context,
avoiding fixed-point overflow for severely overlapping particles. It reuses the
Common force algorithms and L-BFGS solver; no CPU/Reference or replacement force
Context is created. Initialized force parameters and callback identity are
preserved. The backend adapts accumulator types/conversions while retaining
64-bit integer indices and random state. PME and intermediate accumulators use
matching floating producers and consumers. Linked contexts switch together at
a completion boundary; accumulators are cleared and cached forces invalidated
on entry/exit. Allocations keep their eight-byte capacity, with packed floats
used temporarily. Programs cache both pipeline variants and compile the second
only when needed. This is still single precision: energies and forces must
remain representable as finite floats.

When this option is OFF, minimization keeps the OpenCL-style Q32.32 accumulators
and still cannot fall back to CPU/Reference. A GPU diagnostic rejects
out-of-range contributions passing through the fixed-point conversion helper.
This is not a complete overflow detector: direct casts and overflow of an
accumulated sum are not all covered. The OFF path is a limited-range comparison
baseline, not a substitute for the large-force path. Its conversion helpers
also read the diagnostic-enable flag outside minimization, so this build is
not a zero-diagnostic-cost replica of OpenCL.

Float atomic addition independently defaults to an OpenCL-style 32-bit CAS
loop; `OPENMM_METAL_NATIVE_FLOAT_ATOMICS` selects native Metal float addition.
This does not change the accumulator representation or enable float PME grids
for ordinary MD. Separate PME, LJPME, and ConstantPotential switches select
Common's existing floating-grid producer/FFT path, omitting fixed-point spread
conversion and the finishing dispatch. These grid-format switches work with
either CAS or native float atomics; force-buffer storage is unchanged. An LJPME
system uses its own switch for both electrostatic and dispersion grids.
Broad Metal fast math is also OFF by default. Even when
requested, it stays disabled for Common's compensated ConstantPotential CG
solver and Common-source floating-accumulator programs, which require
finite-value checks and scaled-reduction semantics. FP16 neighbor-bound programs
also request precise math: half overflow must remain conservative infinity,
not be optimized under finite-only assumptions. Register-bitonic sort programs
also use precise math to preserve NaN classification and signed-zero payloads.
VkFFT compilation follows
the build-time fast-math setting independently of the force-buffer mode.
API selection follows the running OS, not the SDK or deployment target:
macOS 15+ uses `mathMode` and `mathFloatingPointFunctions`, while macOS 13/14
uses the guarded legacy setter. OFF selects `Safe`/`Precise`, ON selects
`Fast`/`Fast` on the modern API. Runtime version detection does not enable
optimizations or change the build-time switches. Initial compilation, lazy
accumulator variants, and VkFFT all follow this policy. MSL remains 3.0 and
the minimum deployment target remains macOS 13.

The runtime uses private Objective-C++ and C++11 public headers. Only the
private VkFFT translation unit requires C++17/metal-cpp. Backend-5-only fixes
in the bundled VkFFT address ownership, per-dispatch push constants,
threadgroup barriers, and in-place resource binding. Unit-length axes are
handled without a CPU FFT. Apple metal-cpp headers are an external build
dependency; they are not a replacement for the runtime's host API.

Linked CustomCV/ATM contexts share their parent's command queue. Other
same-device contexts and independent queues require explicit synchronization
when sharing arrays. Encoding operations are serialized to accommodate Common
worker-thread uploads. Nonblocking transfers use `getPinnedBuffer()`; do not
read or reuse a pending range until its queue/event completes.

The OFF host path mirrors OpenCL's Apple-device scheduling where possible:

- Neighbor counts are copied into dedicated reusable host storage before force
  dispatches. The later host check waits only for that copy, not the force work.
- Autoclears submit groups of up to six buffers, like OpenCL's fused clears.

These baseline changes do not enable optional algorithm switches or general
command batching. Metal submission costs, VkFFT, and shader compilation still
differ from OpenCL; matching the orchestration does not establish equal speed.

## Build and test

Requires macOS 13+, an arm64 shared build, a macOS 15+ SDK (Xcode 16 or newer),
and Apple's metal-cpp headers exposing the modern math options. The Metal
build is disabled by default. The official
[macOS15.2/iOS18.2 metal-cpp release](https://github.com/apple/metal-cpp/tree/9b4029d993648571595451f3b965c7c438fb7dc4)
provides these APIs; supply its root directory with `METAL_CPP_INCLUDE_DIR`.
Older headers are rejected during configuration instead of failing in VkFFT.

```sh
cmake -S . -B build/metal-common \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_ARCHITECTURES=arm64 \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=13.0 \
  -DOPENMM_BUILD_METAL_LIB=ON \
  -DOPENMM_BUILD_SHARED_LIB=ON \
  -DOPENMM_BUILD_STATIC_LIB=OFF \
  -DMETAL_CPP_INCLUDE_DIR=/path/to/metal-cpp \
  -DBUILD_TESTING=ON
cmake --build build/metal-common -j 8
MTL_DEBUG_LAYER=1 ctest --test-dir build/metal-common -R '^TestMetal' --output-on-failure
```

The core Platform library installs to OpenMM's normal `lib/plugins` directory;
this does not add optional force plugins such as AMOEBA, Drude, or RPMD.

Tests reuse the existing single-precision core assertions and tolerances.
Checkpoint coverage is restricted explicitly to one device; the existing
multi-device test is not applicable. Additional tests cover dynamic platform
loading, Common-source adaptation, arrays/queues, sorting, dense neighbor-list
overflow, FFTs, and fixed-point/floating-accumulator switching. Separate
floating-accumulator force suites exercise the minimizer's shader ABI, including
PME, GB, CustomCV/ATM, callbacks, and virtual sites. No-GPU skips are not successful GPU validation. Existing
tolerances remain provisional regression checks, not a newly agreed numerical
accuracy contract. Test execution times are not performance benchmarks.

## Independent comparison switches

All optional performance fast paths default to `OFF`. GPU-only large-force
minimization and per-operation submission default to `ON`. No unused placeholder
switches are provided. Necessary Metal API/MSL compatibility and correctness
fixes (including the fixed-point atomic emulation and VkFFT integration) are
part of the backend, not optional optimizations.

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
| `OPENMM_METAL_FAST_SPARSE_FORCE_AGGREGATION` | One Q32.32 write per sparse-pair component | Experimental SIMD aggregation of equal-target, already-quantized contributions; floating mode unchanged |
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
| `OPENMM_METAL_FAST_LCPO_NEIGHBOR_SCAN` | Common LCPO shared prefix scan | Hierarchical SIMD prefix scan |
| `OPENMM_METAL_FAST_CUSTOM_MANY_PARTICLE_NEIGHBOR_SCAN` | Common CustomManyParticle shared prefix scan | Hierarchical SIMD prefix scan |
| `OPENMM_METAL_FAST_SORT_BUCKET_SCAN` | OpenCL bucket shared prefix scan | Hierarchical SIMD prefix scan |
| `OPENMM_METAL_FAST_SORT_REGISTER_BITONIC` | Existing short-list sort | Register bitonic sort for supported records of up to 32 elements |

Source transformations match reviewed Common/OpenCL templates, including only
the existing generated-parameter holes. They do not globally define CUDA/HIP
macros. Independent switches preserve the baseline algorithm when OFF; sparse
pairs do not require the shuffle or ballot options. Ballot masks use unsigned
32-bit arithmetic, including the lane-31-only case. FP16 bounds round outward,
including subnormal and overflow boundaries. Periodic block bounds retain
OpenCL's ordered image selection before reducing radii. Register paths preserve
the shared force formulas and Q32.32 global accumulation; their performance and
register pressure must still be measured on each GPU family.
When both short-list sort switches are enabled, register bitonic takes precedence
for supported records with at most 32 elements; the CUDA-style scan selection
still applies to eligible larger lists. Register bitonic places NaN keys after
numeric keys and preserves input order for equal keys and NaNs.

The local MSL preamble exposes SIMD votes, shuffle/rotate, scalar reductions,
prefix sums, and separate int32/float32 atomic interfaces. Uint64 min/max are
enabled only on macOS GPUs advertising Apple8 or newer; their MSL operations
return no previous value. `atomicAddUInt64()` is explicitly unsupported. The
two-word fixed-point accumulation helper is a reduction, not a native fetch-add.
No general HAL API was added to Common Compute or the other backends.
The sparse-force aggregation experiment preserves unsigned modulo-2^64 sums
and has exact-bit carry/sign/wrap tests. It does not enable sparse pairs itself,
change Q32.32 conversion, or redesign the buffer ABI. Unique targets may make it
slower, and it does not supply a general threadgroup force-staging architecture.

Standalone OpenCL stream compaction has no Metal caller; the active neighbor
compaction is covered by the neighbor switches. FFT remains VkFFT, so no
in-house radix optimization is added. Optional AMOEBA/HIPPO kernels remain
outside this core-platform scope. Matrix-hardware neighbor screening and a
general replacement Q32.32 accumulation architecture are not production capabilities.
`OPENMM_METAL_EXPERIMENTAL_MATRIX_SCREEN=ON` builds a separate matrix-screening
test/profiling target; it does not change the MD neighbor list. The prototype
packs each 32-by-32 tile into sixteen 8-by-8 float matrix operations, records
candidate masks, and always performs the scalar minimum-image recheck for the
authoritative result. The candidates alone are not established as conservative,
especially near cutoff boundaries or for periodic images. This experiment is
not an acceleration claim and cannot yet replace exact neighbor screening.
The Apple M4 prototype tests observed missed candidates for wide-coordinate and
periodic tiles; the mandatory scalar recheck restores correctness. Including
that recheck, the prototype was slower in the measured synthetic workloads.

Submission policy does not change the kernels or add a host wait per dispatch.
An FFT is submitted as one logical operation; internal VkFFT stages remain
inside its compute encoder. Events, flushes, blocking transfers, and destruction
preserve completion boundaries. `MetalQueue::flush()` submits without waiting;
Common's `flushQueue()` waits, following the CUDA/HIP interface contract.

For an A/B comparison, configure otherwise identical Release builds, change
one option at a time, warm up compilation, and time the same workload through
its final completion wait. Measure host elapsed time and GPU execution
separately. A CUDA-derived path is an experimental alternative, not a claim
that it is faster on Apple silicon. Other CUDA-specific tuning remains
separate porting work, not an enabled capability of these switches.

### Reproducible timing

`BenchmarkMetal` is an opt-in build target, not a correctness test. It saves a
checkpoint before warmup, warms up compilation, restores the original checkpoint
for each repeat, and reports synchronized wall time plus device/OS/build-switch
metadata, canonical input fingerprints, and integrator type/step size as JSON lines:

```sh
cmake --build build/metal-common --target BenchmarkMetal
build/metal-common/BenchmarkMetal --case pme --particles 4096 --steps 100 --warmup 20 --repeats 5
```

The synthetic cases `cutoff`, `pme`, `ljpme`, and `gbsa` are deterministic
micro-workloads, not representative proteins or cross-generation performance
claims. For real systems, export uncompressed System, State, and optionally
Integrator XML from the existing `examples/benchmarks/benchmark.py` fixtures:

```sh
build/metal-common/BenchmarkMetal --system system.xml --state state.xml \
  --integrator integrator.xml --steps 1000 --warmup 100 --repeats 5
```

Use the same inputs and otherwise identical builds, change one switch, and
repeat both runs. Serialized stochastic integrators and forces must have explicit
nonzero random seeds; automatic seeds are rejected. This fixes input selection,
not bitwise trajectory equivalence across different algorithms or platforms.
Deterministic `CustomIntegrator` inputs are supported. The benchmark rejects
`CustomIntegrator` expressions containing the reserved `uniform` or `gaussian`
variables, including kinetic-energy expressions and nested `CompoundIntegrator`
children: not all of Common's custom random-generator state is checkpointed.
This restriction belongs to the timing harness, not the Metal backend.
Keep Metal validation/capture enabled for correctness runs
but disabled for timing. Wall time includes host dispatch and the final energy
readback; it is not isolated GPU kernel time. Use a Metal GPU capture/counters
for register spills, occupancy, atomic contention, and individual kernel costs.
Measurements on one Apple GPU do not establish results on other generations.
