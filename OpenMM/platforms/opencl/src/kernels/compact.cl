/* Code for CUDA stream compaction. Roughly based on:
    Billeter M, Olsson O, Assarsson U. Efficient Stream Compaction on Wide SIMD Many-Core Architectures.
        High Performance Graphics 2009.

    Notes:
        - paper recommends 128 threads/block, so this is hard coded.
        - I only implement the prefix-sum based compact primitive, and not the POPC one, as that is more
          complicated and performs poorly on current hardware
        - I only implement the scattered- and staged-write variant of phase III as it they have reasonable
          performance across most of the tested workloads in the paper. The selective variant is not
          implemented.
        - The prefix sum of per-block element counts (phase II) is not done in a particularly efficient
          manner. It is, however, done in a very easy to program manner, and integrated into the top of
          phase III, reducing the number of kernel invocations required. If one wanted to use existing code,
          it'd be easy to take the CUDA SDK scanLargeArray sample, and do a prefix sum over dgBlockCounts in
          a phase II kernel. You could also adapt the existing prescan128 to take an initial value, and scan
          dgBlockCounts in stages.

  Date:         23 Aug 2009
  Author:       CUDA version by Imran Haque (ihaque@cs.stanford.edu), converted to OpenCL by Peter Eastman
  Affiliation:  Stanford University
  License:      Public Domain
*/

// Phase 1: Count valid elements per thread block
// Hard-code 128 thd/blk
unsigned int sumReduce128(__local unsigned int* arr) {
    // Parallel reduce element counts
    // Assumes 128 thd/block
    int thread = get_local_id(0);
    if (thread < 64) arr[thread] += arr[thread+64];
    barrier(CLK_LOCAL_MEM_FENCE);
#ifdef WARPS_ARE_ATOMIC
    if (thread < 32) {
        arr[thread] += arr[thread+32];
        if (thread < 16) arr[thread] += arr[thread+16];
        if (thread < 8) arr[thread] += arr[thread+8];
        if (thread < 4) arr[thread] += arr[thread+4];
        if (thread < 2) arr[thread] += arr[thread+2];
        if (thread < 1) arr[thread] += arr[thread+1];
    }
#else
    if (thread < 32) arr[thread] += arr[thread+32];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 16) arr[thread] += arr[thread+16];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 8) arr[thread] += arr[thread+8];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 4) arr[thread] += arr[thread+4];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 2) arr[thread] += arr[thread+2];
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread < 1) arr[thread] += arr[thread+1];
#endif
    barrier(CLK_LOCAL_MEM_FENCE);
    return arr[0];
}

__kernel void countElts(__global unsigned int* restrict dgBlockCounts, __global const unsigned int* restrict dgValid, const unsigned int len, __local unsigned int* restrict dsCount) {
    dsCount[get_local_id(0)] = 0;
    unsigned int ub;
    const unsigned int eltsPerBlock = len/get_num_groups(0) + ((len % get_num_groups(0)) ? 1 : 0);
    ub = (len < (get_group_id(0)+1)*eltsPerBlock) ? len : ((get_group_id(0) + 1)*eltsPerBlock);
    for (int base = get_group_id(0) * eltsPerBlock; base < (get_group_id(0)+1)*eltsPerBlock; base += get_local_size(0)) {
        if ((base + get_local_id(0)) < ub && dgValid[base+get_local_id(0)])
            dsCount[get_local_id(0)]++;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    unsigned int blockCount = sumReduce128(dsCount);
    if (get_local_id(0) == 0) dgBlockCounts[get_group_id(0)] = blockCount;
    return;
}

// Phase 2/3: Move valid elements using SIMD compaction (phase 2 is done implicitly at top of __global__ method)
// Exclusive prefix scan over 128 elements
// Assumes 128 threads
// Taken from cuda SDK "scan" sample for naive scan, with small modifications
int exclusivePrescan128(__local const unsigned int* in, __local unsigned int* outAndTemp) {
    const int n=128;
    //TODO: this temp storage could be reduced since we write to shared memory in out anyway, and n is hardcoded
    //__shared__ int temp[2*n];
    __local unsigned int* temp = outAndTemp;
    int pout = 1, pin = 0;

    // load input into temp
    // This is exclusive scan, so shift right by one and set first elt to 0
    int thread = get_local_id(0);
    temp[pout*n + get_local_id(0)] = (get_local_id(0) > 0) ? in[get_local_id(0)-1] : 0;
    barrier(CLK_LOCAL_MEM_FENCE);

    for (int offset = 1; offset < n; offset *= 2)
    {
        pout = 1 - pout; // swap double buffer indices
        pin  = 1 - pout;
        barrier(CLK_LOCAL_MEM_FENCE);
        temp[pout*n+get_local_id(0)] = temp[pin*n+get_local_id(0)];
        if (get_local_id(0) >= offset)
            temp[pout*n+get_local_id(0)] += temp[pin*n+get_local_id(0) - offset];
    }

    //out[get_local_id(0)] = temp[pout*n+get_local_id(0)]; // write output
    barrier(CLK_LOCAL_MEM_FENCE);
    return outAndTemp[127]+in[127]; // Return sum of all elements
}

int compactSIMDPrefixSum(__local const unsigned int* dsData, __local const unsigned int* dsValid, __local unsigned int* dsCompact, __local unsigned int* dsLocalIndex) {
    int numValid = exclusivePrescan128(dsValid,dsLocalIndex);
    int thread = get_local_id(0);
    if (dsValid[get_local_id(0)]) dsCompact[dsLocalIndex[get_local_id(0)]] = dsData[get_local_id(0)];
    return numValid;
}

__kernel void moveValidElementsStaged(__global const unsigned int* restrict dgData, __global unsigned int* restrict dgCompact, __global const unsigned int* restrict dgValid,
            __global const unsigned int* restrict dgBlockCounts, unsigned int len, __global unsigned int* restrict dNumValidElements,
            __local unsigned int* restrict inBlock, __local unsigned int* restrict validBlock, __local unsigned int* restrict compactBlock) {
    __local unsigned int dsLocalIndex[256];
    int blockOutOffset=0;
    // Sum up the blockCounts before us to find our offset
    // This is totally inefficient - lots of repeated work b/w blocks, and uneven balancing.
    // Paper implements this as a prefix sum kernel in phase II
    // May still be faster than an extra kernel invocation?
    int thread = get_local_id(0);
    for (int base = 0; base < get_group_id(0); base += get_local_size(0)) {
        // Load up the count of valid elements for each block before us in batches of 128
        if ((base + get_local_id(0)) < get_group_id(0)) {
            validBlock[get_local_id(0)] = dgBlockCounts[base+get_local_id(0)];
        } else {
            validBlock[get_local_id(0)] = 0;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        // Parallel reduce these counts
        // Accumulate in the final offset variable
        blockOutOffset += sumReduce128(validBlock);
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    unsigned int ub;
    const unsigned int eltsPerBlock = len/get_num_groups(0) + ((len % get_num_groups(0)) ? 1 : 0);
    ub = (len < (get_group_id(0)+1)*eltsPerBlock) ? len : ((get_group_id(0) + 1)*eltsPerBlock);
    for (int base = get_group_id(0) * eltsPerBlock; base < (get_group_id(0)+1)*eltsPerBlock; base += get_local_size(0)) {
        if ((base + get_local_id(0)) < ub) {
            validBlock[get_local_id(0)] = dgValid[base+get_local_id(0)];
            inBlock[get_local_id(0)] = dgData[base+get_local_id(0)];
        } else {
            validBlock[get_local_id(0)] = 0;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        int numValidBlock = compactSIMDPrefixSum(inBlock,validBlock,compactBlock,dsLocalIndex);
        barrier(CLK_LOCAL_MEM_FENCE);
        if (get_local_id(0) < numValidBlock) {
            dgCompact[blockOutOffset + get_local_id(0)] = compactBlock[get_local_id(0)];
        }
        blockOutOffset += numValidBlock;
    }
    if (get_group_id(0) == (get_num_groups(0)-1) && get_local_id(0) == 0) {
        *dNumValidElements = blockOutOffset;
    }
}
