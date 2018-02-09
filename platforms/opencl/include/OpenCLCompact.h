#ifndef __OPENMM_OPENCLCOMPACT_H__
#define __OPENMM_OPENCLCOMPACT_H__

/* Code for OPENCL stream compaction. Roughly based on:
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

#include "OpenCLArray.h"
#include "OpenCLContext.h"

namespace OpenMM {

class OPENMM_EXPORT_OPENCL OpenCLCompact {
public:
    OpenCLCompact(OpenCLContext& context);
    void compactStream(OpenCLArray& dOut, OpenCLArray& dIn, OpenCLArray& dValid, OpenCLArray& numValid);
private:
    OpenCLContext& context;
    OpenCLArray dgBlockCounts;
    cl::Kernel countKernel;
    cl::Kernel moveValidKernel;
};

} // namespace OpenMM

#endif // __OPENMM_OPENCLCOMPACT_H__
