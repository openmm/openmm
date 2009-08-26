
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
  Author:       Imran Haque (ihaque@cs.stanford.edu)
  Affiliation:  Stanford University
  License:      Public Domain
*/

#include "cudaCompact.h"

typedef unsigned int T;

// Phase 1: Count valid elements per thread block
// Hard-code 128 thd/blk
__device__ unsigned int sumReduce128(unsigned int* arr) {
    // Parallel reduce element counts
    // Assumes 128 thd/block
    if (threadIdx.x < 64) arr[threadIdx.x] += arr[threadIdx.x+64];
    __syncthreads();
    if (threadIdx.x < 32) {
        arr[threadIdx.x] += arr[threadIdx.x+32];
        if (threadIdx.x < 16) arr[threadIdx.x] += arr[threadIdx.x+16];
        if (threadIdx.x < 8) arr[threadIdx.x] += arr[threadIdx.x+8];
        if (threadIdx.x < 4) arr[threadIdx.x] += arr[threadIdx.x+4];
        if (threadIdx.x < 2) arr[threadIdx.x] += arr[threadIdx.x+2];
        if (threadIdx.x < 1) arr[threadIdx.x] += arr[threadIdx.x+1];
    }
    __syncthreads();
    return arr[0];
}

__global__ void countElts(unsigned int* dgBlockCounts,const unsigned int* dgValid,const size_t eltsPerBlock,const size_t len) {
    __shared__ unsigned int dsCount[128];
    dsCount[threadIdx.x] = 0;
    size_t ub;
    ub = (len < (blockIdx.x+1)*eltsPerBlock) ? len : ((blockIdx.x + 1)*eltsPerBlock);
    for (int base = blockIdx.x * eltsPerBlock; base < (blockIdx.x+1)*eltsPerBlock; base += blockDim.x) {
        if ((base + threadIdx.x) < ub && dgValid[base+threadIdx.x])
            dsCount[threadIdx.x]++;
    }
    __syncthreads();
    unsigned int blockCount = sumReduce128(dsCount);
    if (threadIdx.x == 0) dgBlockCounts[blockIdx.x] = blockCount;
    return;
}

// Phase 2/3: Move valid elements using SIMD compaction (phase 2 is done implicitly at top of __global__ method)
// Exclusive prefix scan over 128 elements
// Assumes 128 threads
// Taken from cuda SDK "scan" sample for naive scan, with small modifications
__device__ int exclusivePrescan128(const unsigned int* in,unsigned int* outAndTemp) {
    const int n=128;
    //TODO: this temp storage could be reduced since we write to shared memory in out anyway, and n is hardcoded
    //__shared__ int temp[2*n];
    unsigned int* temp = outAndTemp;
    int pout = 1, pin = 0;

    // load input into temp
    // This is exclusive scan, so shift right by one and set first elt to 0
    temp[pout*n + threadIdx.x] = (threadIdx.x > 0) ? in[threadIdx.x-1] : 0;
    __syncthreads();

    for (int offset = 1; offset < n; offset *= 2)
    {
        pout = 1 - pout; // swap double buffer indices
        pin  = 1 - pout;
        __syncthreads();
        temp[pout*n+threadIdx.x] = temp[pin*n+threadIdx.x];
        if (threadIdx.x >= offset)
            temp[pout*n+threadIdx.x] += temp[pin*n+threadIdx.x - offset];
    }

    //out[threadIdx.x] = temp[pout*n+threadIdx.x]; // write output
    __syncthreads();
    return outAndTemp[127]+in[127]; // Return sum of all elements
}
__device__ int compactSIMDPrefixSum(const T* dsData,const unsigned int* dsValid,T* dsCompact) {
    __shared__ unsigned int dsLocalIndex[256];
    int numValid = exclusivePrescan128(dsValid,dsLocalIndex);
    if (dsValid[threadIdx.x]) dsCompact[dsLocalIndex[threadIdx.x]] = dsData[threadIdx.x];
    return numValid;
}

__global__ void moveValidElementsStaged(const T* dgData,T* dgCompact,const unsigned int* dgValid,const unsigned int* dgBlockCounts,size_t eltsPerBlock,size_t len,size_t* dNumValidElements) {
    __shared__ T inBlock[128];
    __shared__ unsigned int validBlock[128];
    __shared__ T compactBlock[128];
    int blockOutOffset=0;
    // Sum up the blockCounts before us to find our offset
    // This is totally inefficient - lots of repeated work b/w blocks, and uneven balancing.
    // Paper implements this as a prefix sum kernel in phase II
    // May still be faster than an extra kernel invocation?
    for (int base = 0; base < blockIdx.x; base += blockDim.x) {
        // Load up the count of valid elements for each block before us in batches of 128
        if ((base + threadIdx.x) < blockIdx.x) {
            validBlock[threadIdx.x] = dgBlockCounts[base+threadIdx.x];
        } else {
            validBlock[threadIdx.x] = 0;
        }
        __syncthreads();
        // Parallel reduce these counts
        // Accumulate in the final offset variable
        blockOutOffset += sumReduce128(validBlock);
    }

    size_t ub;
    ub = (len < (blockIdx.x+1)*eltsPerBlock) ? len : ((blockIdx.x + 1)*eltsPerBlock);
    for (int base = blockIdx.x * eltsPerBlock; base < (blockIdx.x+1)*eltsPerBlock; base += blockDim.x) {
        if ((base + threadIdx.x) < ub) {
            validBlock[threadIdx.x] = dgValid[base+threadIdx.x];
            inBlock[threadIdx.x] = dgData[base+threadIdx.x];
        } else {
            validBlock[threadIdx.x] = 0;
        }
        __syncthreads();
        int numValidBlock = compactSIMDPrefixSum(inBlock,validBlock,compactBlock);
        __syncthreads();
        if (threadIdx.x < numValidBlock) {
            dgCompact[blockOutOffset + threadIdx.x] = compactBlock[threadIdx.x];
        }
        blockOutOffset += numValidBlock;
    }
    if (blockIdx.x == (gridDim.x-1) && threadIdx.x == 0) {
        *dNumValidElements = blockOutOffset;
    }
}

__global__ void moveValidElementsScattered(const T* dgData,T* dgCompact,const unsigned int* dgValid,const unsigned int* dgBlockCounts,size_t eltsPerBlock,size_t len,size_t* dNumValidElements) {
    __shared__ T inBlock[128];
    __shared__ unsigned int validBlock[128];
    T* compactBlock=dgCompact;
    size_t blockOutOffset = 0;
    // Sum up the blockCounts before us to find our offset
    // This is totally inefficient - lots of repeated work b/w blocks, and uneven balancing.
    // Paper implements this as a prefix sum kernel in phase II
    // May still be faster than an extra kernel invocation?
    for (int base = 0; base < blockIdx.x; base += blockDim.x) {
        // Load up the count of valid elements for each block before us in batches of 128
        if ((base + threadIdx.x) < blockIdx.x) {
            validBlock[threadIdx.x] = dgBlockCounts[base+threadIdx.x];
        } else {
            validBlock[threadIdx.x] = 0;
        }
        __syncthreads();
        // Parallel reduce these counts
        // Accumulate in the final offset variable
        blockOutOffset += sumReduce128(validBlock);
    }
    compactBlock += blockOutOffset;
    size_t ub;
    ub = (len < (blockIdx.x+1)*eltsPerBlock) ? len : ((blockIdx.x + 1)*eltsPerBlock);
    for (int base = blockIdx.x * eltsPerBlock; base < (blockIdx.x+1)*eltsPerBlock; base += blockDim.x) {
        if ((base + threadIdx.x) < ub) {
            validBlock[threadIdx.x] = dgValid[base+threadIdx.x];
            inBlock[threadIdx.x] = dgData[base+threadIdx.x];
        } else {
            validBlock[threadIdx.x] = 0;
        }
        __syncthreads();
        int numValidBlock = compactSIMDPrefixSum(inBlock,validBlock,compactBlock);
        blockOutOffset += numValidBlock;
        compactBlock += numValidBlock;
    }
    if (blockIdx.x == (gridDim.x-1) && threadIdx.x == 0) {
        *dNumValidElements = blockOutOffset;
    }
}

void planCompaction(compactionPlan& d,bool stageOutput) {
    int device;
    cudaGetDevice(&device);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    d.nThreadBlocks = 16*deviceProp.multiProcessorCount;
    cudaMalloc((void**)&(d.dgBlockCounts), d.nThreadBlocks*sizeof(unsigned int));
    d.stageOutput = stageOutput;
    // TODO: make sure allocation worked
    d.valid = true;
}

void destroyCompactionPlan(compactionPlan& d) {
    if (d.valid) cudaFree(d.dgBlockCounts);
}

int compactStream(const compactionPlan& d,T* dOut,const T* dIn,const unsigned int* dValid,size_t len,size_t* dNumValid) {
    if (!d.valid) {
        return -1;
    }
    // Figure out # elements per block
    unsigned int numBlocks = d.nThreadBlocks;
    if (numBlocks*128 > len)
        numBlocks = (len+127)/128;
    const size_t eltsPerBlock = len/numBlocks + ((len % numBlocks) ? 1 : 0);

    // TODO: implement loop over blocks of 10M
    // Phase 1: Calculate number of valid elements per thread block
    countElts<<<numBlocks,128>>>(d.dgBlockCounts,dValid,eltsPerBlock,len);

    // Phase 2/3: Move valid elements using SIMD compaction
    if (d.stageOutput) {
        moveValidElementsStaged<<<numBlocks,128>>>(dIn,dOut,dValid,d.dgBlockCounts,eltsPerBlock,len,dNumValid);
    } else {
        moveValidElementsScattered<<<numBlocks,128>>>(dIn,dOut,dValid,d.dgBlockCounts,eltsPerBlock,len,dNumValid);
    }
    return 0;
}
