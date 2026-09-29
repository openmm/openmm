# Reorder, neighbor-list and direct-launch port

Ported onto the current PR checkout at `773ca815`, retaining upstream source changes. The implementation is in eight owned files; all new experimental behavior is opt-in. After implementation, two CPU regression programs were compiled and run against the actual ported headers and current upstream reference source. See README.md for combined validation and its limitations.

## Feature accounting

| Live feature | Current-master disposition | Gate and fallback |
| --- | --- | --- |
| `OPENMM_EXPERIMENT_REORDER_GATHER` | Common reorder dispatch and two protected virtual hooks ported in `ComputeContext.h/.cpp`; CUDA hook implementation is included in CudaContext. New-to-old physical permutation and padded `-1` sentinels are retained. | Common default hook returns false, so unsupported platforms preserve all original host downloads/copies/uploads. CUDA preparation must succeed before any downloads are skipped. |
| `OPENMM_EXPERIMENT_REORDER_RADIX_SORT` | Four stable byte passes retained, with reusable per-reorder scratch and the signed-key transform. Original `(bin, i)` tie ordering is preserved because `i` is constructed in increasing order. | Exact `1`, at least 4096 molecules; default `std::sort`, including allocation failure before mutation. |
| `OPENMM_EXPERIMENT_REORDER_HILBERT_LUT` | Exact 3D/8-bit lookup specialization retained in the new `ReorderHilbert.h`; bundled Rice University notices preserved. | Exact `1`, mixed precision, one context, device gather active. Coordinates outside unsigned 8-bit range and disabled path call the bundled generic function. |
| `OPENMM_EXPERIMENT_NEIGHBOR_REUSE` | Separate reset/displacement-precheck kernels, conditional sort and neighbor-list metadata work, periodic-box representation tracking, and reorder/restore listener retained. Current direct-kernel xyz bounds remain refreshed even while the list is reused. | Changed live default-on to explicit `1` for this contribution. Requires mixed precision, one context, more than 3000 atom blocks, cutoff + periodic + padding + neighbor list, and standard nonbonded kernel source. OFF uses the upstream call sequence. |
| `OPENMM_EXPERIMENT_NEIGHBOR_FORCED_REBUILD` | Mandatory-rebuild branch skips the two precheck launches; a by-value host flag makes all bounds CTAs take the same branch and one thread publishes the GPU flag for subsequent kernels. | Exact `1` plus active neighbor reuse. Otherwise reset + precheck sequence is retained. |
| `OPENMM_EXPERIMENT_NEIGHBOR_SKIN_FRACTION` | Fixed-value policy retained in new `NeighborSkinPolicy.h`. | Active neighbor reuse and exactly one ordinary PME `NonbondedForce` with only the explicitly recognized companion forces. `0` and `0.08` mean baseline 0.08; `0.10` and `0.12` select experimental fractions. Other values throw only when the experiment is eligible. Unsupported systems keep 0.08. |
| `OPENMM_EXPERIMENT_DIRECT_BLOCKS_PER_SM` | Separate direct-force block count, larger energy-buffer requirement, and widened exclusion-tile partition arithmetic retained. | Only exact `5` or `6`, compute capability 12.0, mixed precision, more than 90000 atoms. Default force block count stays 4 per SM; custom kernel sources use the original block count. |
| `OPENMM_EXPERIMENT_DIRECT_CUTOFF_GUARD` | Implemented in CommonCalcNonbondedForce, the new location of the legacy CudaKernels code. | Preserve the legacy force/source/precision guards; verify in the combined integration review. |

## Current-master adaptations

- `CudaNonbondedUtilities` now owns a `ComputeSort` shared pointer and uses the common `ComputeSortImpl::SortTrait`. The opt-in conditional sorter uses `ComputeSort(new CudaSort(..., &rebuildNeighborList))`; the disabled path remains `context.createSort(...)`. CudaSort provides the optional execution-flag constructor.
- Shared arrays now expose `ArrayInterface`; the new displacement-precheck argument uses `context.unwrap(context.getPosq()).getDevicePointer()`.
- Current master computes neighbor lists at the global `maxCutoff`, rather than per-group cutoff with a `lastCutoff` transition. The obsolete legacy `lastCutoff` comparison and assignment were deliberately not resurrected. All newly generated precheck/list kernels use the same current `maxCutoff` and padding.
- Current master's lazy force parameter argument initialization, half-precision bounding-box storage, shared ownership, and newer Common initialization methods were retained.
- `CudaContext` already caps launches at six blocks per SM for the eligible architecture, and sizes energy storage using `nonbonded->getNumEnergyBuffers()`. The revised getter therefore supplies the required larger direct-force capacity without changing that launch cap.
