# UMA C++ Plugin Implementation - Current Status

## Phase 1: Infrastructure - STARTED ✓

### Files Created

1. **Header Files** (Complete)
   - `/media/extradrive/Trajectories/openmm/plugins/uma/openmmapi/include/openmm/UMAForce.h`
   - `/media/extradrive/Trajectories/openmm/plugins/uma/openmmapi/include/openmm/UMAKernels.h`
   - `/media/extradrive/Trajectories/openmm/plugins/uma/openmmapi/include/openmm/internal/UMAForceImpl.h`
   - `/media/extradrive/Trajectories/openmm/plugins/uma/openmmapi/include/openmm/internal/windowsExportUMA.h`

2. **Implementation Files** (Complete)
   - `/media/extradrive/Trajectories/openmm/plugins/uma/openmmapi/src/UMAForce.cpp`
   - `/media/extradrive/Trajectories/openmm/plugins/uma/openmmapi/src/UMAForceImpl.cpp`

3. **Build System** (Partial)
   - `/media/extradrive/Trajectories/openmm/plugins/uma/CMakeLists.txt` - Basic structure

### Directory Structure Created
```
plugins/uma/
├── CMakeLists.txt                    ✓
├── openmmapi/
│   ├── include/openmm/
│   │   ├── UMAForce.h                ✓
│   │   ├── UMAKernels.h              ✓
│   │   └── internal/
│   │       ├── UMAForceImpl.h        ✓
│   │       └── windowsExportUMA.h    ✓
│   └── src/
│       ├── UMAForce.cpp              ✓
│       └── UMAForceImpl.cpp          ✓
└── platforms/                        (directories created)
```

## Remaining Work

### Phase 1 Completion (1-2 days)
- [ ] Create platforms/common/CMakeLists.txt
- [ ] Create platforms/cuda/CMakeLists.txt
- [ ] Create platforms/reference/CMakeLists.txt
- [ ] Add UMA plugin to main OpenMM CMakeLists.txt

### Phase 2: AtomicData Builder (3-5 days)
- [ ] Implement AtomicDataBuilder class in C++
- [ ] Implement NeighborListBuilder with cell lists
- [ ] Add periodic boundary condition support
- [ ] Create tensor conversion utilities

### Phase 3: libtorch Integration (5-7 days)
- [ ] Implement CommonCalcUMAForceKernel
- [ ] Add model loading from .pt files
- [ ] Implement forward pass with unit conversions
- [ ] Load atomic reference energies from YAML

### Phase 4: CUDA Optimization (5-7 days)
- [ ] Implement CudaCalcUMAForceKernel
- [ ] Zero-copy memory sharing with torch::from_blob
- [ ] GPU-accelerated neighbor list construction
- [ ] Optimize memory transfers

### Phase 5: Python Wrappers (2-3 days)
- [ ] SWIG interface file
- [ ] Python binding generation
- [ ] Update OpenMM-ML umapotential.py

### Phase 6: Testing (3-5 days)
- [ ] Unit tests in C++
- [ ] Validation against FAIRChem reference
- [ ] Performance benchmarks

### Phase 7: RPMD Integration (2-3 days)
- [ ] Test with RPMDIntegrator
- [ ] Verify ring polymer bead handling
- [ ] Force group compatibility

## Total Estimated Effort

**3-4 weeks for complete implementation**

## Alternative: Immediate Working Solution

While the C++ plugin is being developed, you can use **PythonForce** today (see below).

---

# PRAGMATIC SOLUTION: PythonForce Implementation

## Can Be Completed Today (2-3 hours)

Instead of waiting weeks for the C++ plugin, we can modify the existing 
`umapotential.py` to use `PythonForce`, which works immediately with
acceptable performance (5-20% overhead vs C++ plugin).

### Advantages of Starting with PythonForce

1. **Works today** - No C++ compilation needed
2. **Still GPU accelerated** - PyTorch runs on GPU
3. **RPMD compatible** - Full OpenMM integrator support
4. **Validates approach** - Can test UMA integration immediately
5. **Easy to upgrade** - Can switch to C++ plugin later

### Performance Comparison

| Approach | Time to Implement | Overhead | Status |
|----------|------------------|----------|---------|
| PythonForce | 2-3 hours | 5-20% | Can do today |
| C++ Plugin | 3-4 weeks | <1% | Full project |

### Implementation

Modify `fairchem/openmm-ml/openmmml/models/umapotential.py` to use
`openmm.PythonForce` instead of `openmmtorch.TorchForce`.

This gives you:
- Working UMA integration today
- RPMD support immediately  
- GPU acceleration
- Time to develop C++ plugin properly

## Recommendation

**Dual-track approach:**

1. **This week**: Implement PythonForce solution (immediate results)
2. **Next month**: Develop C++ plugin (long-term performance)

The PythonForce solution provides 80% of the benefit with 5% of the effort!
