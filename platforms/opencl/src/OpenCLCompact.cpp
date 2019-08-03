

/* Code for OpenCL stream compaction. Roughly based on:
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

#include "OpenCLCompact.h"
#include "OpenCLKernelSources.h"

using namespace OpenMM;

OpenCLCompact::OpenCLCompact(OpenCLContext& context) : context(context) {
    dgBlockCounts.initialize<cl_uint>(context, context.getNumThreadBlocks(), "dgBlockCounts");
    cl::Program program = context.createProgram(OpenCLKernelSources::compact);
    countKernel = cl::Kernel(program, "countElts");
    moveValidKernel = cl::Kernel(program, "moveValidElementsStaged");
}

void OpenCLCompact::compactStream(OpenCLArray& dOut, OpenCLArray& dIn, OpenCLArray& dValid, OpenCLArray& numValid) {
    // Figure out # elements per block
    unsigned int len = dIn.getSize();
    unsigned int numBlocks = context.getNumThreadBlocks();
    if (numBlocks*128 > len)
        numBlocks = (len+127)/128;
    const size_t eltsPerBlock = len/numBlocks + ((len % numBlocks) ? 1 : 0);

    // TODO: implement loop over blocks of 10M
    // Phase 1: Calculate number of valid elements per thread block
    countKernel.setArg<cl::Buffer>(0, dgBlockCounts.getDeviceBuffer());
    countKernel.setArg<cl::Buffer>(1, dValid.getDeviceBuffer());
    countKernel.setArg<cl_uint>(2, len);
    countKernel.setArg(3, 128*sizeof(cl_uint), NULL);
    context.executeKernel(countKernel, len, 128);

    // Phase 2/3: Move valid elements using SIMD compaction
    moveValidKernel.setArg<cl::Buffer>(0, dIn.getDeviceBuffer());
    moveValidKernel.setArg<cl::Buffer>(1, dOut.getDeviceBuffer());
    moveValidKernel.setArg<cl::Buffer>(2, dValid.getDeviceBuffer());
    moveValidKernel.setArg<cl::Buffer>(3, dgBlockCounts.getDeviceBuffer());
    moveValidKernel.setArg<cl_uint>(4, len);
    moveValidKernel.setArg<cl::Buffer>(5, numValid.getDeviceBuffer());
    moveValidKernel.setArg(6, 128*sizeof(cl_uint), NULL);
    moveValidKernel.setArg(7, 128*sizeof(cl_uint), NULL);
    moveValidKernel.setArg(8, 128*sizeof(cl_uint), NULL);
    context.executeKernel(moveValidKernel, len, 128);
}
