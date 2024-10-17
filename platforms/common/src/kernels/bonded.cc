inline DEVICE void saveBondedForce(GLOBAL mm_ulong* RESTRICT forceBuffer, unsigned int atom, real3 force) {
    ATOMIC_ADD(&forceBuffer[atom], (mm_ulong) realToFixedPoint(force.x));
    ATOMIC_ADD(&forceBuffer[atom+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force.y));
    ATOMIC_ADD(&forceBuffer[atom+PADDED_NUM_ATOMS*2], (mm_ulong) realToFixedPoint(force.z));
}

// Several threads of one warp write forces of the same atom quite often.
// Head segmented reduction decreases the number of such intra-warp atomic conflicts and the
// total number of atomic adds.

#if defined(USE_HIP) && !defined(USE_DOUBLE_PRECISION)

inline DEVICE inline real3 headSegmentedSum(real3 input, unsigned int key, bool* isLeader) {
    const unsigned int laneId = (threadIdx.x & (warpSize - 1));
    const unsigned int prevLaneKey = __shfl(key, laneId - 1);
    const bool head = laneId == 0 || key != prevLaneKey;

#if defined(AMD_RDNA)
    const unsigned int n = __popc(__ballot(1));
    unsigned int flags = __ballot(head);
    flags >>= 1;
    flags &= ~((static_cast<unsigned int>(1) << laneId) - 1);
    flags |= static_cast<unsigned int>(1) << (n - 1);
    const unsigned int nextSegmentStart = __ffs(flags);
#else
    const unsigned int n = __popcll(__ballot(1));
    unsigned long long flags = __ballot(head);
    flags >>= 1;
    flags &= ~((static_cast<unsigned long long>(1) << laneId) - 1);
    flags |= static_cast<unsigned long long>(1) << (n - 1);
    const unsigned int nextSegmentStart = __ffsll(flags);
#endif

    real3 output = input;
    *isLeader = head;

    unsigned int offset = 1;
    while (__ballot(laneId + offset < nextSegmentStart)) {
        const real3 other = warpShuffle(output, laneId + offset);
        const real k = real(laneId + offset < nextSegmentStart);
        output += k * other;
        offset *= 2;
    }
    return output;
}

inline DEVICE void saveBondedForceWithConflicts(unsigned int bond, unsigned int numBonds, GLOBAL mm_ulong* RESTRICT forceBuffer, unsigned int atom, real3 force) {
    (void)bond;
    (void)numBonds;
    bool isLeader;
    force = headSegmentedSum(force, atom, &isLeader);
    if (isLeader) {
        saveBondedForce(forceBuffer, atom, force);
    }
}

#elif defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 700 && !defined(USE_DOUBLE_PRECISION)

inline DEVICE real3 headSegmentedSum(unsigned int mask, real3 input, unsigned int key, bool* outHead) {
    const unsigned int laneId = (threadIdx.x & (warpSize - 1));
    const unsigned int prevLaneKey = __shfl_up_sync(mask, key, 1);
    const bool head = laneId == 0 || key != prevLaneKey;

    const unsigned int n = __popc(mask);
    unsigned int flags = __ballot_sync(mask, head);
    flags >>= 1;
    flags &= ~((static_cast<unsigned int>(1) << laneId) - 1);
    flags |= static_cast<unsigned int>(1) << (n - 1);
    const unsigned int nextSegmentStart = __ffs(flags);

    real3 output = input;
    *outHead = head;

    unsigned int offset = 1;
    while (__any_sync(mask, laneId + offset < nextSegmentStart)) {
        const real k = real(laneId + offset < nextSegmentStart);
        output.x += k * __shfl_down_sync(mask, output.x, offset),
        output.y += k * __shfl_down_sync(mask, output.y, offset),
        output.z += k * __shfl_down_sync(mask, output.z, offset);
        offset *= 2;
    }
    return output;
}

inline DEVICE void saveBondedForceWithConflicts(unsigned int bond, unsigned int numBonds, GLOBAL mm_ulong* RESTRICT forceBuffer, unsigned int atom, real3 force) {
    const unsigned int activeLanes = numBonds - bond / warpSize * warpSize;
    const unsigned int activeMask = activeLanes >= warpSize ? 0xffffffff : ((1u << activeLanes) - 1);
    bool isLeader;
    force = headSegmentedSum(activeMask, force, atom, &isLeader);
    if (isLeader) {
        saveBondedForce(forceBuffer, atom, force);
    }
}

#else

inline DEVICE void saveBondedForceWithConflicts(unsigned int bond, unsigned int numBonds, GLOBAL mm_ulong* RESTRICT forceBuffer, unsigned int atom, real3 force) {
    (void)bond;
    (void)numBonds;
    saveBondedForce(forceBuffer, atom, force);
}

#endif
