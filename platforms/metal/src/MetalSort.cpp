/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform.
 * Source: platforms/opencl/src/OpenCLSort.cpp
 * Optional parallel short-list selection follows platforms/cuda/src/CudaSort.cpp.
 *
 * Original OpenCL Platform code:
 * Portions copyright (c) 2010-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 *
 * Metal Platform code:
 * Portions copyright (c) 2026 Chun-Chi Hung.
 * Authors: Chun-Chi Hung
 *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#ifdef _MSC_VER
    // Prevent Windows from defining macros that interfere with other code.
    #define NOMINMAX
#endif
#include "MetalSort.h"
#include "MetalKernel.h"
#include "MetalOpenCLKernelSources.h"
#include "MetalReductionOptimizations.h"
#include <algorithm>
#include <map>
#include <string>

#ifndef OPENMM_METAL_FAST_SHORT_LIST_SORT
#define OPENMM_METAL_FAST_SHORT_LIST_SORT 0
#endif

using namespace OpenMM;
using namespace std;

MetalSort::MetalSort(MetalContext& context, ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform) :
        context(context), trait(trait), dataLength(length), useShortList2(false), uniform(uniform) {
    if (trait == nullptr || trait->getDataSize() <= 0 || trait->getKeySize() <= 0)
        throw OpenMMException("MetalSort requires a valid nonempty sort trait");
    if (length < 2)
        return;
    // Create kernels.

    std::map<std::string, std::string> replacements;
    replacements["DATA_TYPE"] = trait->getDataType();
    replacements["KEY_TYPE"] =  trait->getKeyType();
    replacements["SORT_KEY"] = trait->getSortKey();
    replacements["MIN_KEY"] = trait->getMinKey();
    replacements["MAX_KEY"] = trait->getMaxKey();
    replacements["MAX_VALUE"] = trait->getMaxValue();
    replacements["UNIFORM"] = (uniform ? "1" : "0");
    MetalReductionOptimizations::Settings optimizations = MetalReductionOptimizations::getBuildSettings();
    // Apply before substituting DATA_TYPE/KEY_TYPE so the exact OpenCL
    // fingerprint remains available to the narrowly scoped transformer.
    string source = MetalReductionOptimizations::apply(MetalOpenCLKernelSources::sort, optimizations);
    map<string, string> defines;
    ComputeProgram program = context.compileProgram(context.replaceStrings(source, replacements), defines);
    shortListKernel = program->createKernel("sortShortList");
    for (int i = 0; i < 3; i++)
        shortListKernel->addArg();
    computeRangeKernel = program->createKernel("computeRange");
    for (int i = 0; i < 7; i++)
        computeRangeKernel->addArg();
    assignElementsKernel = program->createKernel(uniform ? "assignElementsToBuckets" : "assignElementsToBuckets2");
    for (int i = 0; i < 7; i++)
        assignElementsKernel->addArg();
    computeBucketPositionsKernel = program->createKernel("computeBucketPositions");
    for (int i = 0; i < 3; i++)
        computeBucketPositionsKernel->addArg();
    copyToBucketsKernel = program->createKernel("copyDataToBuckets");
    for (int i = 0; i < 6; i++)
        copyToBucketsKernel->addArg();
    sortBucketsKernel = program->createKernel("sortBuckets");
    for (int i = 0; i < 5; i++)
        sortBucketsKernel->addArg();

    // Work out the work group sizes for various kernels.

    unsigned int maxGroupSize = std::min(256, (int) context.getMaxThreadBlockSize());
    int maxSharedMem = context.getMaxThreadgroupMemory();
    unsigned int maxRangeSize = std::min(maxGroupSize, (unsigned int) computeRangeKernel->getMaxBlockSize());
    unsigned int maxPositionsSize = std::min(maxGroupSize, (unsigned int) computeBucketPositionsKernel->getMaxBlockSize());
    int maxLocalBuffer = (maxSharedMem/trait->getDataSize())/2;
    maxRangeSize = min(maxRangeSize, (unsigned int) (maxSharedMem/(2*trait->getKeySize())));
    if (maxLocalBuffer < 2 || maxRangeSize == 0 || maxPositionsSize == 0)
        throw OpenMMException("The Metal sort trait requires too much threadgroup memory");
    int maxShortList = min(1024, maxLocalBuffer);
    isShortList = (length <= maxShortList);
#if OPENMM_METAL_FAST_SHORT_LIST_SORT
    // CUDA's alternate short-list algorithm scans in parallel instead of sorting
    // within one threadgroup. Reuse its existing OpenCL source verbatim. The
    // kernel has a fixed 64-element local tile and no grid-stride outer loop.
    useShortList2 = (length <= min(3000, MetalContext::ThreadBlockSize*context.getNumThreadBlocks()) &&
            64*trait->getDataSize() <= maxSharedMem);
    if (useShortList2) {
        shortList2Kernel = program->createKernel("sortShortList2");
        if (shortList2Kernel->getMaxBlockSize() < 64)
            throw OpenMMException("The Metal short-list scan requires a 64-thread group");
        for (int i = 0; i < 3; i++)
            shortList2Kernel->addArg();
    }
#endif
    for (rangeKernelSize = 1; rangeKernelSize*2 <= maxRangeSize; rangeKernelSize *= 2)
        ;
    positionsKernelSize = std::min(rangeKernelSize, maxPositionsSize);
    sortKernelSize = (isShortList ? rangeKernelSize : rangeKernelSize/2);
    unsigned int maxSortSize = (isShortList ? shortListKernel->getMaxBlockSize() : sortBucketsKernel->getMaxBlockSize());
    while (sortKernelSize > maxSortSize || sortKernelSize > maxLocalBuffer)
        sortKernelSize /= 2;
    if (sortKernelSize == 0)
        throw OpenMMException("The Metal device cannot execute the sorting kernel");
    if (rangeKernelSize > length)
        rangeKernelSize = length;
    unsigned int targetBucketSize = max(1u, sortKernelSize/2);
    unsigned int numBuckets = length/targetBucketSize;
    if (numBuckets < 1)
        numBuckets = 1;
    if (positionsKernelSize > numBuckets)
        positionsKernelSize = numBuckets;

    // Create workspace arrays.

    dataRange.initialize(context, 2, trait->getKeySize(), "sortDataRange");
    bucketOffset.initialize<unsigned int>(context, numBuckets, "bucketOffset");
    bucketOfElement.initialize<unsigned int>(context, length, "bucketOfElement");
    offsetInBucket.initialize<unsigned int>(context, length, "offsetInBucket");
    buckets.initialize(context, length, trait->getDataSize(), "buckets");
}

MetalSort::~MetalSort() {
}

void MetalSort::sort(ArrayInterface& data) {
    if (data.getSize() != dataLength || data.getElementSize() != trait->getDataSize())
        throw OpenMMException("MetalSort called with different data size");
    if (data.getSize() < 2)
        return;
    ArrayInterface& cldata = data;
    if (useShortList2) {
        shortList2Kernel->setArg(0, cldata);
        shortList2Kernel->setArg(1, buckets);
        shortList2Kernel->setArg(2, (int) dataLength);
        shortList2Kernel->execute(dataLength, 64);
        buckets.copyTo(cldata);
        return;
    }
    if (isShortList) {
        shortListKernel->setArg(0, cldata);
        shortListKernel->setArg(1, dataLength);
        static_cast<MetalKernel&>(*shortListKernel).setLocalArg(2, dataLength*trait->getDataSize());
        shortListKernel->execute(sortKernelSize, sortKernelSize);
        return;
    }

    // Compute the range of data values.

    unsigned int numBuckets = bucketOffset.getSize();
    computeRangeKernel->setArg(0, cldata);
    computeRangeKernel->setArg(1, (unsigned int) (cldata.getSize()));
    computeRangeKernel->setArg(2, dataRange);
    static_cast<MetalKernel&>(*computeRangeKernel).setLocalArg(3, rangeKernelSize*trait->getKeySize());
    static_cast<MetalKernel&>(*computeRangeKernel).setLocalArg(4, rangeKernelSize*trait->getKeySize());
    computeRangeKernel->setArg(5, (int) (numBuckets));
    computeRangeKernel->setArg(6, bucketOffset);
    computeRangeKernel->execute(rangeKernelSize, rangeKernelSize);

    // Assign array elements to buckets.

    assignElementsKernel->setArg(0, cldata);
    assignElementsKernel->setArg(1, (int) (cldata.getSize()));
    assignElementsKernel->setArg(2, (int) (numBuckets));
    assignElementsKernel->setArg(3, dataRange);
    assignElementsKernel->setArg(4, bucketOffset);
    assignElementsKernel->setArg(5, bucketOfElement);
    assignElementsKernel->setArg(6, offsetInBucket);
    assignElementsKernel->execute(cldata.getSize());

    // Compute the position of each bucket.

    computeBucketPositionsKernel->setArg(0, (int) (numBuckets));
    computeBucketPositionsKernel->setArg(1, bucketOffset);
    static_cast<MetalKernel&>(*computeBucketPositionsKernel).setLocalArg(2, positionsKernelSize*sizeof(int));
    computeBucketPositionsKernel->execute(positionsKernelSize, positionsKernelSize);

    // Copy the data into the buckets.

    copyToBucketsKernel->setArg(0, cldata);
    copyToBucketsKernel->setArg(1, buckets);
    copyToBucketsKernel->setArg(2, (int) (cldata.getSize()));
    copyToBucketsKernel->setArg(3, bucketOffset);
    copyToBucketsKernel->setArg(4, bucketOfElement);
    copyToBucketsKernel->setArg(5, offsetInBucket);
    copyToBucketsKernel->execute(cldata.getSize());

    // Sort each bucket.

    sortBucketsKernel->setArg(0, cldata);
    sortBucketsKernel->setArg(1, buckets);
    sortBucketsKernel->setArg(2, (int) (numBuckets));
    sortBucketsKernel->setArg(3, bucketOffset);
    static_cast<MetalKernel&>(*sortBucketsKernel).setLocalArg(4, sortKernelSize*trait->getDataSize());
    sortBucketsKernel->execute(((cldata.getSize()+sortKernelSize-1)/sortKernelSize)*sortKernelSize, sortKernelSize);
}
