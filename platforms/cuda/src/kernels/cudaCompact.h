#ifndef __OPENMM_CUDACOMPACT_H__
#define __OPENMM_CUDACOMPACT_H__

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

struct compactionPlan {
    bool valid;
    unsigned int* dgBlockCounts;
    unsigned int nThreadBlocks;
    bool stageOutput;
};

extern "C"
void planCompaction(compactionPlan& d,bool stageOutput=true);

extern "C"
void destroyCompactionPlan(compactionPlan& d);

extern "C"
int compactStream(const compactionPlan& d,unsigned int* dOut,const unsigned int* dIn,const unsigned int* dValid,size_t len,size_t* dNumValid);

#endif // __OPENMM_CUDACOMPACT_H__