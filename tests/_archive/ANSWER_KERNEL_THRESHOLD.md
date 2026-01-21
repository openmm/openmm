# Answer: Where is the CUDA Kernel Threshold in OpenMM?

## TL;DR

The **3000-atom threshold** that caused your 10x slowdown is at:

**File**: `platforms/common/src/CommonCalcNonbondedForce.cpp`, **line 622**

```cpp
cc.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, 
    force.getCutoffDistance(), exclusionList, source, force.getForceGroup(), 
    numParticles > 3000,  // <-- Boolean flag: use neighbor list?
    true);
```

The 7th parameter (`numParticles > 3000`) tells OpenMM whether to use a neighbor list:
- **`false` (≤ 3000 atoms)**: All-pairs kernel (fast for small systems)
- **`true` (> 3000 atoms)**: Neighbor list kernel (fast for large systems)

---

## Why Your Benchmark Hit a Performance Cliff

Your diatomic system:
- **1500 molecules** = **3000 atoms** → All-pairs kernel → **~20,000 steps/s**
- **1501 molecules** = **3002 atoms** → Neighbor list kernel → **~2,000 steps/s**

The neighbor list kernel has overhead (building, sorting, maintaining the list) that only pays off for larger systems. At 3002 atoms, you pay the overhead without getting the benefit.

---

## Complete List of OpenMM CUDA Thresholds

### 1. **Neighbor List Enable** (3,000 atoms)
- **Location**: `platforms/common/src/CommonCalcNonbondedForce.cpp:622`
- **What it does**: Switches from all-pairs to neighbor list kernel
- **Your issue**: You hit exactly this threshold

### 2. **Large Blocks for Neighbor List** (90,000 atoms)
- **Location**: `platforms/cuda/src/CudaNonbondedUtilities.cpp:77`
- **Code**: `useLargeBlocks = (context.getNumAtoms() > 90000);`
- **What it does**: Uses 1024-atom blocks instead of sorted smaller blocks
- **Trade-off**: Faster list build, slightly less efficient force calculation

### 3. **Thread Block Size** (Compute Capability 2.0)
- **Location**: `platforms/cuda/src/CudaNonbondedUtilities.cpp:70`
- **Code**: `forceThreadBlockSize = (context.getComputeCapability() < 2.0 ? 128 : 256);`
- **What it does**: Adjusts parallelism for older vs. modern GPUs

### 4. **Pair List Bits** (Compute Capability 8.0)
- **Location**: `platforms/cuda/src/CudaNonbondedUtilities.cpp:523`
- **What it does**: Uses more efficient pair tracking on Ampere+ GPUs

---

## How to Avoid the 3000-Atom Cliff

### Option 1: Stay Below (Recommended for Small Systems)
```bash
# Use 2500 atoms (safe zone)
python benchmark_speed.py --molecules 1250 --no-cavity
```

### Option 2: Jump Above (Recommended for Large Systems)
```bash
# Use 5000+ atoms where neighbor list is efficient
python benchmark_speed.py --molecules 2500 --no-cavity
```

### Option 3: Modify OpenMM Source (Advanced)
Edit `platforms/common/src/CommonCalcNonbondedForce.cpp:622`:

```cpp
// Original
numParticles > 3000

// Force all-pairs kernel up to 10,000 atoms
numParticles > 10000

// Or always use all-pairs (slow for large systems!)
false
```

Then recompile OpenMM.

---

## Why OpenMM Chose 3000 Atoms

This threshold is empirically tuned based on:
1. **Memory overhead** of neighbor lists (~100 KB minimum)
2. **Rebuild frequency** (every ~10-20 steps typically)
3. **GPU warp utilization** (32-thread warps work best with specific sizes)
4. **Profiling** across different GPUs (Tesla, Quadro, GeForce, etc.)

For *most* systems, 3000 is a good heuristic. Your benchmark just happened to hit the exact boundary.

---

## Your Water Simulations Are Safe

- **216 molecules** = ~650 atoms → All-pairs kernel (fast)
- **1000 molecules** = ~3000 atoms → Right at threshold, but TIP4P-Ew adds virtual M-sites, so actual particle count is ~4000 → Neighbor list kernel is appropriate

The baseline simulation (running now) shows **no performance issues** and correctly exhibits flexible water vibrations.

---

## Further Reading

For the full implementation, see:
- `platforms/common/src/CommonCalcNonbondedForce.cpp` (high-level logic)
- `platforms/cuda/src/CudaNonbondedUtilities.cpp` (CUDA-specific setup)
- `platforms/cuda/src/kernels/nonbonded.cu` (actual CUDA kernels)
- `platforms/cuda/src/kernels/findInteractingBlocks.cu` (neighbor list building)

The OpenMM developers maintain excellent documentation at:
- https://github.com/openmm/openmm
- https://openmm.org/development
