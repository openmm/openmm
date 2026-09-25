# Experimental Metal Platform

- Implements `ComputeContext`, `ArrayInterface`, `ComputeQueue`, `ComputeEvent`,
  `ComputeProgram`, and `ComputeKernel` for one Apple-silicon GPU.
- Follows the CUDA/HIP class responsibilities, array operations, argument
  rebinding, and block-based launch flow.
- Reuses the existing Common context, array wrapper, and force-info sources.
- Uses private Objective-C++ for Metal API calls and ordinary C++11 headers.
- Compiles native MSL at runtime with Metal 3.0; no offline library or capability
  probe is required.

The Metal Platform currently contains only the Common runtime foundation and
is not yet registered as an OpenMM simulation Platform. It does not yet execute
existing Common kernel sources or implement force,
integrator, sorting, neighbor-list, or FFT/PME paths. Unsupported interfaces
throw explicitly. There is no separate Metal force algorithm, fixed-point
emulation, or shader translator.

## Porting boundaries

`MetalArray`, `MetalKernel`, `MetalProgram`, `MetalQueue`, and `MetalEvent`
follow their `Cuda*`/`Hip*` counterparts. `MetalContext` includes the runtime and
single-precision state-buffer portion, without the force-dependent platform
initialization. The force buffer uses ordinary 8-byte signed integer elements
and the existing three-component-plane layout; allocating storage does not
enable 64-bit atomic accumulation.

The required Metal adaptations are resource ownership, blit transfers, command
buffers/events, runtime compilation, and encoder bindings. Nonblocking uploads
and downloads use the context's `getPinnedBuffer()` (or a range within it), like
CUDA's page-locked-memory requirement. Do not read or reuse that memory until
the relevant queue/event completes. Blocking transfers accept ordinary host
memory. Operations select the current queue at execution time, not at array
creation. Submission to each queue must be serialized by the caller.

Native kernel buffer arguments use sequential `[[buffer(i)]]` slots in Common
argument order; scalars use constant-buffer references. Stage builtins are
explicit MSL parameters. The caller supplies the matching signature and legal,
non-overlapping buffer bindings. Launches use complete threadgroups (64 threads
by default), capped at a conservative 128 groups, so kernels must support the
existing grid-stride convention. This cap is not performance tuning.

## Build and test

Configure an arm64 shared-library build on macOS 13 or newer with
`OPENMM_BUILD_METAL_LIB=ON` and `BUILD_TESTING=ON`. The Metal option is off by
default; static builds and installation as a simulation plugin are not part
of this initial target.

```sh
cmake -S . -B build/metal-common \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_ARCHITECTURES=arm64 \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=13.0 \
  -DOPENMM_BUILD_METAL_LIB=ON \
  -DOPENMM_BUILD_SHARED_LIB=ON \
  -DOPENMM_BUILD_STATIC_LIB=OFF \
  -DBUILD_TESTING=ON
cmake --build build/metal-common --target TestMetalComputeContext
ctest --test-dir build/metal-common -R '^TestMetalComputeContext$' --output-on-failure
```

The test exercises the Common C++ interfaces using standalone, test-only MSL
smoke kernels. It covers transfers, resize/rebinding, full and partial logical
blocks, grid-stride execution, clearing, queue/event ordering, pinned-memory
readback, and invalid inputs. No GPU returns skip code 77, not a successful GPU
validation. Passing this test does not establish simulation or Common shader
compatibility.
