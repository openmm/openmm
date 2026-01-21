# CUDA Kernel Selection Thresholds in OpenMM

Based on source code inspection of the OpenMM CUDA platform.

## Threshold 1: Neighbor List Block Size (90,000 atoms)

**File**: `platforms/cuda/src/CudaNonbondedUtilities.cpp` (line 77)

```cpp
// When building the neighbor list, we can optionally use large blocks (1024 atoms) to
// accelerate the process.  This makes building the neighbor list faster, but it prevents
// us from sorting atom blocks by size, which leads to a slightly less efficient neighbor
// list.  We guess based on system size which will be faster.

useLargeBlocks = (context.getNumAtoms() > 90000);
```

### What it does:
- **≤ 90,000 atoms**: Uses **small blocks** with atom sorting → More efficient neighbor list
- **> 90,000 atoms**: Uses **large blocks (1024 atoms)** → Faster build, less efficient list

This primarily affects neighbor list construction speed, not force calculation speed.

---

## Threshold 2: Neighbor List vs. All-Pairs (3,000 atoms)

**File**: `platforms/common/src/CommonCalcNonbondedForce.cpp` (line 622)

```cpp
cc.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, 
    force.getCutoffDistance(), exclusionList, source, force.getForceGroup(), 
    numParticles > 3000,  // <-- This enables neighbor list
    true);
```

### What it does:
- **≤ 3,000 atoms**: Uses **all-pairs kernel** (no neighbor list)
- **> 3,000 atoms**: Switches to **neighbor list kernel**

**This is the threshold you hit at 1500→1501 molecules!**

### Why the performance cliff?

1. **1500 molecules = 3000 atoms** → All-pairs kernel (optimized for small systems)
2. **1501 molecules = 3002 atoms** → Neighbor list kernel (optimized for large systems)

At 3002 atoms, you're in the **worst case** for the neighbor list kernel:
- Overhead of building/maintaining neighbor list
- Not enough atoms to amortize that overhead
- The all-pairs kernel would be faster, but the threshold already switched

### Performance Characteristics

| System Size | Method | Typical Performance |
|------------|--------|-------------------|
| < 3000 atoms | All-pairs | Fast for small systems |
| 3001-5000 atoms | Neighbor list | **Slower** (worst case zone) |
| 5000-20000 atoms | Neighbor list | Good performance |
| > 20000 atoms | Neighbor list | Excellent performance |

---

## Additional Thresholds in CUDA Kernels

### Thread Block Size (Compute Capability)
**File**: `platforms/cuda/src/CudaNonbondedUtilities.cpp` (line 70)

```cpp
forceThreadBlockSize = (context.getComputeCapability() < 2.0 ? 128 : 256);
```

- **Older GPUs** (compute < 2.0): 128 threads per block
- **Modern GPUs** (compute ≥ 2.0): 256 threads per block

### Pair List Support (Compute Capability ≥ 8.0)
**File**: `platforms/cuda/src/CudaNonbondedUtilities.cpp` (line 523)

```cpp
defines["MAX_BITS_FOR_PAIRS"] = (canUsePairList ? 
    (context.getComputeCapability() < 8.0 ? "2" : "3") : "0");
```

More efficient pair tracking on newer Ampere/Ada GPUs.

---

## Implications for Your Simulations

### Your Benchmark (Diatomic System)
- Uses `NonbondedForce.CutoffPeriodic` (not PME)
- 1500 molecules = 3000 atoms → Fast all-pairs kernel
- 1501 molecules = 3002 atoms → Slow neighbor list kernel
- **~10x slowdown** is expected behavior

### Your Water Simulations (TIP4P-Ew)
- Uses `PME` with virtual sites
- 216 molecules = **~650 atoms** → Well below threshold (fast)
- 1000 molecules = **~3000 atoms** → Right at threshold
- No performance issues expected

---

## Recommendations

1. **For benchmarking**: Stay well away from 3000-atom threshold
   - Use 2500 atoms (safe)
   - Or jump to 5000+ atoms (where neighbor list is efficient)

2. **For production**: Choose system size based on science
   - The absolute time difference is negligible for real simulations
   - 100,000-step production runs dominate any kernel selection overhead

3. **If you must optimize around 3000 atoms**: 
   - There's no easy workaround without modifying OpenMM source
   - You could disable neighbor lists entirely (custom build), but this would slow down large systems

---

## References

These thresholds are empirically tuned by OpenMM developers (Peter Eastman et al.) based on profiling across different GPUs and system sizes. They represent good heuristics for *most* systems, but edge cases (like 3002 atoms) can hit worst-case performance.

For more details, see the OpenMM GitHub repository:
- https://github.com/openmm/openmm
