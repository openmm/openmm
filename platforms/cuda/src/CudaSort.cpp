/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2010-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "CudaSort.h"
#include "CudaKernelSources.h"
#include <algorithm>
#include <map>

using namespace OpenMM;
using namespace std;

CudaSort::CudaSort(CudaContext& context, ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform) :
        CudaSort(context, trait, length, uniform, -1) {
}

CudaSort::CudaSort(CudaContext& context, ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform,
        int knownIntegerMaximum) :
        context(context), trait(trait), dataLength(length), uniform(uniform),
        hasFixedRange(knownIntegerMaximum >= 0) {
    if (knownIntegerMaximum != -1 && (knownIntegerMaximum < 1 || knownIntegerMaximum > 16777215 ||
            !uniform || string(trait->getKeyType()) != "int"))
        throw OpenMMException("Known-range CudaSort requires uniform int keys and maximum in [1, 16777215]");
    // Create kernels.

    map<string, string> replacements;
    replacements["DATA_TYPE"] = trait->getDataType();
    replacements["KEY_TYPE"] =  trait->getKeyType();
    replacements["SORT_KEY"] = trait->getSortKey();
    replacements["MIN_KEY"] = trait->getMinKey();
    replacements["MAX_KEY"] = trait->getMaxKey();
    replacements["MAX_VALUE"] = trait->getMaxValue();
    replacements["UNIFORM"] = (uniform ? "1" : "0");
    replacements["FIXED_RANGE_ENABLED"] = (knownIntegerMaximum >= 0 ? "1" : "0");
    replacements["FIXED_RANGE_UPPER"] = context.intToString(knownIntegerMaximum);
    CUmodule module = context.createModule(context.replaceStrings(CudaKernelSources::sort, replacements));
    shortListKernel = context.getKernel(module, "sortShortList");
    shortList2Kernel = context.getKernel(module, "sortShortList2");
    computeRangeKernel = context.getKernel(module, "computeRange");
    assignElementsKernel = context.getKernel(module, uniform ? "assignElementsToBuckets" : "assignElementsToBuckets2");
    computeBucketPositionsKernel = context.getKernel(module, "computeBucketPositions");
    copyToBucketsKernel = context.getKernel(module, "copyDataToBuckets");
    sortBucketsKernel = context.getKernel(module, "sortBuckets");
    scatterPhysicalIndicesKernel = (hasFixedRange ? context.getKernel(module, "scatterPhysicalIndices") : NULL);

    // Work out the work group sizes for various kernels.

    int maxBlockSize;
    cuDeviceGetAttribute(&maxBlockSize, CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X, context.getDevice());
    int maxSharedMem;
    cuDeviceGetAttribute(&maxSharedMem, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, context.getDevice());
    int maxLocalBuffer = (maxSharedMem/trait->getDataSize())/2;
    int maxShortList = min(3000, max(maxLocalBuffer, CudaContext::ThreadBlockSize*context.getNumThreadBlocks()));
    isShortList = (length <= maxShortList);
    for (rangeKernelSize = 1; rangeKernelSize*2 <= maxBlockSize; rangeKernelSize *= 2)
        ;
    positionsKernelSize = rangeKernelSize;
    sortKernelSize = (isShortList ? rangeKernelSize/2 : rangeKernelSize/4);
    if (rangeKernelSize > length)
        rangeKernelSize = length;
    if (sortKernelSize > maxLocalBuffer)
        sortKernelSize = maxLocalBuffer;
    unsigned int targetBucketSize = sortKernelSize/2;
    unsigned int numBuckets = length/targetBucketSize;
    if (numBuckets < 1)
        numBuckets = 1;
    if (positionsKernelSize > numBuckets)
        positionsKernelSize = numBuckets;

    // Create workspace arrays.

    if (!isShortList) {
        dataRange.initialize(context, 2, trait->getKeySize(), "sortDataRange");
        bucketOffset.initialize<uint1>(context, numBuckets, "bucketOffset");
        bucketOfElement.initialize<uint1>(context, length, "bucketOfElement");
        offsetInBucket.initialize<uint1>(context, length, "offsetInBucket");
    }
    buckets.initialize(context, length, trait->getDataSize(), "buckets");
}

CudaSort::~CudaSort() {
    delete trait;
}

void CudaSort::sort(ArrayInterface& data) {
    // Preserve validation and zero-length behavior before unwrapping the array.
    if (data.getSize() != dataLength || data.getElementSize() != trait->getDataSize())
        throw OpenMMException("CudaSort called with different data size");
    if (data.getSize() == 0)
        return;
    sortImpl(context.unwrap(data), true);
}

void CudaSort::bucketize(CudaArray& data) {
    sortImpl(data, false);
}

bool CudaSort::sortWithGeneratedKeys(CudaArray& data, ComputeKernel generator,
        bool sortWithinBuckets, bool directPhysicalIndices) {
    if (data.getSize() != dataLength || data.getElementSize() != trait->getDataSize())
        throw OpenMMException("CudaSort called with different data size");
    if (data.getSize() == 0)
        return true;
    if (!hasFixedRange || isShortList || !uniform || !generator ||
            dataLength > 0x7fffffffu || (directPhysicalIndices && (sortWithinBuckets || scatterPhysicalIndicesKernel == NULL)) ||
            trait->getDataSize() != sizeof(int2) || trait->getKeySize() != sizeof(int) ||
            string(trait->getDataType()) != "int2" || string(trait->getKeyType()) != "int" ||
            string(trait->getSortKey()) != "value.y")
        return false;
    // The PME generator reserves the 15-argument contract before this call.
    sortImpl(data, sortWithinBuckets, generator, directPhysicalIndices);
    return true;
}

void CudaSort::sortImpl(CudaArray& data, bool sortWithinBuckets,
        ComputeKernel commonGenerator, bool directPhysicalIndices) {
    if (data.getSize() != dataLength || data.getElementSize() != trait->getDataSize())
        throw OpenMMException("CudaSort called with different data size");
    if (data.getSize() == 0)
        return;
    if (isShortList) {
        // We can use a simpler sort kernel that does the entire operation in one kernel.
        
        if (dataLength <= CudaContext::ThreadBlockSize*context.getNumThreadBlocks()) {
            void* sortArgs[] = {&data.getDevicePointer(), &buckets.getDevicePointer(), &dataLength};
            context.executeKernel(shortList2Kernel, sortArgs, dataLength);
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
        void* rangeArgs[] = {&data.getDevicePointer(), &dataLength, &dataRange.getDevicePointer(), &numBuckets, &bucketOffset.getDevicePointer()};
        context.executeKernel(computeRangeKernel, rangeArgs, rangeKernelSize, rangeKernelSize, 2*rangeKernelSize*trait->getKeySize());

        // Assign array elements to buckets.

        if (commonGenerator) {
            commonGenerator->setArg(10, numBuckets);
            commonGenerator->setArg(11, dataRange);
            commonGenerator->setArg(12, bucketOffset);
            commonGenerator->setArg(13, bucketOfElement);
            commonGenerator->setArg(14, offsetInBucket);
            commonGenerator->execute(data.getSize(), 128);
        }
        else {
            void* elementsArgs[] = {&data.getDevicePointer(), &dataLength, &numBuckets, &dataRange.getDevicePointer(),
                    &bucketOffset.getDevicePointer(), &bucketOfElement.getDevicePointer(), &offsetInBucket.getDevicePointer()};
            context.executeKernel(assignElementsKernel, elementsArgs, data.getSize(), 128);
        }

        // Compute the position of each bucket.

        void* computeArgs[] = {&numBuckets, &bucketOffset.getDevicePointer()};
        context.executeKernel(computeBucketPositionsKernel, computeArgs, positionsKernelSize, positionsKernelSize, positionsKernelSize*sizeof(int));

        if (directPhysicalIndices) {
            // The generator left only per-physical-index metadata. Read no data
            // elements while scattering, so output aliases no scatter input.
            void* directArgs[] = {&data.getDevicePointer(), &dataLength, &bucketOffset.getDevicePointer(),
                    &bucketOfElement.getDevicePointer(), &offsetInBucket.getDevicePointer()};
            context.executeKernel(scatterPhysicalIndicesKernel, directArgs, data.getSize());
            return;
        }

        // Copy the data into the buckets.

        void* copyArgs[] = {&data.getDevicePointer(), &buckets.getDevicePointer(), &dataLength, &bucketOffset.getDevicePointer(),
                &bucketOfElement.getDevicePointer(), &offsetInBucket.getDevicePointer()};
        context.executeKernel(copyToBucketsKernel, copyArgs, data.getSize());

        if (sortWithinBuckets) {
            // Preserve the complete sorting contract for every existing caller.
            void* sortArgs[] = {&data.getDevicePointer(), &buckets.getDevicePointer(), &numBuckets, &bucketOffset.getDevicePointer()};
            context.executeKernel(sortBucketsKernel, sortArgs, ((data.getSize()+sortKernelSize-1)/sortKernelSize)*sortKernelSize, sortKernelSize, sortKernelSize*trait->getDataSize());
        }
        else {
            // Partition already contains every element exactly once, including
            // buckets larger than one block. Copy on the current PME stream;
            // preserve both array ownership and all bound device pointers.
            buckets.copyTo(data);
        }
    }
}
