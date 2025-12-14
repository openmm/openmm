/**
 * This file contains the device function for using cross-lane operations (ballot and shuffle)
 */

template<class T>
static __inline__ __device__
T warpShuffle(const T& input, const int src_lane) {
    static_assert(sizeof(T) % sizeof(int) == 0, "incorrect type size");
    constexpr int words_no = sizeof(T) / sizeof(int);

    T output;
    #pragma unroll
    for (int i = 0; i < words_no; i++) {
        int word;
        __builtin_memcpy(&word, reinterpret_cast<const char*>(&input) + i * sizeof(int), sizeof(int));
        word = __builtin_amdgcn_ds_bpermute(src_lane << 2, word);
        __builtin_memcpy(reinterpret_cast<char*>(&output) + i * sizeof(int), &word, sizeof(int));
    }

    return output;
}

#if defined(AMD_RDNA)

// RDNA: device warp size = tile size = 32

#define SHFL(var, srcLane) __shfl(var, srcLane)
#define BALLOT(var) __ballot(var)

template<class T>
static __inline__ __device__
T warpRotateLeft(const T& input) {
    return warpShuffle(input, threadIdx.x + 1);
}

#else

// CDNA: device warp size = 64, tile size = 32

#define SHFL(var, srcLane) __shfl(var, (srcLane) & (32 - 1), 32)
#define BALLOT(var) (unsigned int)(__ballot(var) >> (threadIdx.x & ((64 - 1) ^ (32 - 1))))

template<class T>
static __inline__ __device__
T warpRotateLeft(const T& input) {
    return warpShuffle(input, ((threadIdx.x + 1) & (32 - 1)) | (threadIdx.x & ~(32 - 1)));
}

#endif
