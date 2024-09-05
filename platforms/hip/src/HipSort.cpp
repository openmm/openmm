/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2021 Stanford University and the Authors.      *
 * Portions copyright (c) 2020-2023 Advanced Micro Devices, Inc.              *
 * Authors: Peter Eastman, Nicholas Curtis                                    *
 * Contributors:                                                              *
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

#include "HipSort.h"
#include "HipKernelSources.h"
#include <algorithm>
#include <map>

using namespace OpenMM;
using namespace std;

HipSort::HipSort(HipContext& context, SortTrait* trait, unsigned int length, bool uniform) :
        context(context), trait(trait), dataLength(length), uniform(uniform) {
    // Create kernels.

    map<string, string> replacements;
    replacements["DATA_TYPE"] = trait->getDataType();
    replacements["KEY_TYPE"] =  trait->getKeyType();
    replacements["SORT_KEY"] = trait->getSortKey();
    replacements["MIN_KEY"] = trait->getMinKey();
    replacements["MAX_KEY"] = trait->getMaxKey();
    replacements["MAX_VALUE"] = trait->getMaxValue();
    replacements["UNIFORM"] = (uniform ? "1" : "0");
    hipModule_t module = context.createModule(context.replaceStrings(HipKernelSources::sort, replacements));
    shortListKernel = context.getKernel(module, "sortShortList");
    shortList2Kernel = context.getKernel(module, "sortShortList2");
    computeRangeKernel = context.getKernel(module, "computeRange");
    assignElementsKernel = context.getKernel(module, uniform ? "assignElementsToBuckets" : "assignElementsToBuckets2");
    computeBucketPositionsKernel = context.getKernel(module, "computeBucketPositions");
    copyToBucketsKernel = context.getKernel(module, "copyDataToBuckets");
    sortBucketsKernel = context.getKernel(module, "sortBuckets");

    // Work out the work group sizes for various kernels.

    int maxSharedMem;
    hipDeviceGetAttribute(&maxSharedMem, hipDeviceAttributeMaxSharedMemoryPerBlock, context.getDevice());
    int maxLocalBuffer = (maxSharedMem/trait->getDataSize())/2;
    int maxShortList = min(1024, max(maxLocalBuffer, HipContext::ThreadBlockSize*context.getNumThreadBlocks()));
    isShortList = (length <= maxShortList);
    sortKernelSize = 256;
    rangeKernelSize = 256;
    if (rangeKernelSize > length)
        rangeKernelSize = length;
    rangeKernelBlocks = (length + rangeKernelSize - 1) / rangeKernelSize;
    if (sortKernelSize > maxLocalBuffer)
        sortKernelSize = maxLocalBuffer;
    unsigned int targetBucketSize = uniform ? sortKernelSize/2 : sortKernelSize/8;
    unsigned int numBuckets = length/targetBucketSize;
    if (numBuckets < 1)
        numBuckets = 1;
    // computeBucketPositions is executed as a single work group so larger block size is faster.
    positionsKernelSize = 1024;
    if (positionsKernelSize > numBuckets)
        positionsKernelSize = numBuckets;

    // Create workspace arrays.

    if (!isShortList) {
        counters.initialize<unsigned int>(context, 1, "counters");
        unsigned int zero = 0;
        counters.upload(&zero);
        dataRange.initialize(context, 2*rangeKernelBlocks, trait->getKeySize(), "sortDataRange");
        bucketOffset.initialize<uint1>(context, numBuckets, "bucketOffset");
        bucketOfElement.initialize<uint1>(context, length, "bucketOfElement");
        offsetInBucket.initialize<uint1>(context, length, "offsetInBucket");
    }
    buckets.initialize(context, length, trait->getDataSize(), "buckets");
}

HipSort::~HipSort() {
    delete trait;
}

void HipSort::sort(HipArray& data) {
    if (data.getSize() != dataLength || data.getElementSize() != trait->getDataSize())
        throw OpenMMException("HipSort called with different data size");
    if (data.getSize() == 0)
        return;
    if (isShortList) {
        // We can use a simpler sort kernel that does the entire operation in one kernel.

        if (dataLength <= HipContext::ThreadBlockSize*context.getNumThreadBlocks()) {
            void* sortArgs[] = {&data.getDevicePointer(), &buckets.getDevicePointer(), &dataLength};
            context.executeKernel(shortList2Kernel, sortArgs, dataLength, HipContext::ThreadBlockSize, HipContext::ThreadBlockSize*trait->getKeySize());
            buckets.copyTo(data);
        }
        else {
            void* sortArgs[] = {&data.getDevicePointer(), &dataLength};
            context.executeKernel(shortListKernel, sortArgs, sortKernelSize, sortKernelSize, dataLength*trait->getDataSize());
        }
    }
    else {
        // Compute the range of data values.

        unsigned int numBuckets = bucketOffset.getSize();
        void* rangeArgs[] = {&data.getDevicePointer(), &dataLength, &dataRange.getDevicePointer(), &numBuckets, &bucketOffset.getDevicePointer(), &counters.getDevicePointer()};
        context.executeKernel(computeRangeKernel, rangeArgs, rangeKernelBlocks*rangeKernelSize, rangeKernelSize, 2*rangeKernelSize*trait->getKeySize());

        // Assign array elements to buckets.

        void* elementsArgs[] = {&data.getDevicePointer(), &dataLength, &numBuckets, &dataRange.getDevicePointer(),
                &bucketOffset.getDevicePointer(), &bucketOfElement.getDevicePointer(), &offsetInBucket.getDevicePointer()};
        context.executeKernel(assignElementsKernel, elementsArgs, data.getSize(), 128);

        // Compute the position of each bucket.

        void* computeArgs[] = {&numBuckets, &bucketOffset.getDevicePointer(), &counters.getDevicePointer()};
        context.executeKernel(computeBucketPositionsKernel, computeArgs, positionsKernelSize, positionsKernelSize, positionsKernelSize*sizeof(int));

        // Copy the data into the buckets.

        void* copyArgs[] = {&data.getDevicePointer(), &buckets.getDevicePointer(), &dataLength, &bucketOffset.getDevicePointer(),
                &bucketOfElement.getDevicePointer(), &offsetInBucket.getDevicePointer()};
        context.executeKernel(copyToBucketsKernel, copyArgs, data.getSize());

        // Sort each bucket.

        void* sortArgs[] = {&data.getDevicePointer(), &buckets.getDevicePointer(), &bucketOffset.getDevicePointer()};
        context.executeKernelFlat(sortBucketsKernel, sortArgs, numBuckets*sortKernelSize, sortKernelSize, sortKernelSize*trait->getDataSize());
    }
}
