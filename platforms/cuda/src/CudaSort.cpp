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
        CudaSort(context, trait, length, uniform, NULL) {
}

CudaSort::CudaSort(CudaContext& context, ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform, CudaArray* executionFlag) :
        CudaSort(context, trait, length, uniform, executionFlag, -1) {
}

CudaSort::CudaSort(CudaContext& context, ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform,
        CudaArray* executionFlag, int knownIntegerMaximum) :
        executionFlag(executionFlag), context(context), trait(trait), dataLength(length), uniform(uniform),
        hasFixedRange(knownIntegerMaximum >= 0) {
    if (knownIntegerMaximum != -1 && (knownIntegerMaximum < 1 || knownIntegerMaximum > 16777215 ||
            !uniform || executionFlag != NULL || string(trait->getKeyType()) != "int"))
        throw OpenMMException("Known-range CudaSort requires unguarded uniform int keys and maximum in [1, 16777215]");
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
    replacements["GUARDED_SORT_ARGUMENTS"] = (executionFlag == NULL ? "" : ", const int* __restrict__ executionFlag");
    replacements["GUARDED_SORT_BODY"] = (executionFlag == NULL ? "" : "if (executionFlag[0] == 0) return;");
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
    if (executionFlag != NULL && (isShortList || executionFlag->getSize() < 1 || executionFlag->getElementSize() != sizeof(int) ||
            &executionFlag->getContext() != &context))
        throw OpenMMException("Conditional CudaSort requires a long list and a device int flag");
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
    if (!hasFixedRange || isShortList || executionFlag != NULL || !uniform || !generator ||
            dataLength > 0x7fffffffu || (directPhysicalIndices && (sortWithinBuckets || scatterPhysicalIndicesKernel == NULL)) ||
            trait->getDataSize() != sizeof(int2) || trait->getKeySize() != sizeof(int) ||
            string(trait->getDataType()) != "int2" || string(trait->getKeyType()) != "int" ||
            string(trait->getSortKey()) != "value.y")
        return false;
    // This overload is reserved for the PME generators, whose 15-argument
    // contract is established before invoking it. Preserve the raw API above.
    sortImpl(data, sortWithinBuckets, NULL, NULL, 0, directPhysicalIndices, generator);
    return true;
}

void CudaSort::sortImpl(CudaArray& data, bool sortWithinBuckets) {
    sortImpl(data, sortWithinBuckets, NULL, NULL, 0);
}

bool CudaSort::sortWithGeneratedKeys(CudaArray& data, CUfunction generator, void* const* arguments,
        int numArguments, bool sortWithinBuckets) {
    if (data.getSize() != dataLength || data.getElementSize() != trait->getDataSize())
        throw OpenMMException("CudaSort called with different data size");
    if (data.getSize() == 0)
        return true;
    // Every rejection precedes range/reset: dynamic range would read unwritten keys.
    // Five workspace arguments are appended to the fixed-size host stack array.
    if (!hasFixedRange || isShortList || executionFlag != NULL || !uniform ||
            generator == NULL || arguments == NULL || numArguments < 0 || numArguments > 27 ||
            trait->getDataSize() != sizeof(int2) || trait->getKeySize() != sizeof(int) ||
            string(trait->getDataType()) != "int2" || string(trait->getKeyType()) != "int" ||
            string(trait->getSortKey()) != "value.y")
        return false;
    for (int i = 0; i < numArguments; i++)
        if (arguments[i] == NULL)
            return false;
    sortImpl(data, sortWithinBuckets, generator, arguments, numArguments);
    return true;
}

bool CudaSort::bucketizeGeneratedPhysicalIndices(CudaArray& data, CUfunction generator,
        void* const* arguments, int numArguments) {
    if (data.getSize() != dataLength || data.getElementSize() != trait->getDataSize())
        throw OpenMMException("CudaSort called with different data size");
    if (data.getSize() == 0)
        return true;
    // No reset or generator may run until every physical-index contract is checked.
    if (!hasFixedRange || isShortList || executionFlag != NULL || !uniform || dataLength > 0x7fffffffu ||
            generator == NULL || scatterPhysicalIndicesKernel == NULL || arguments == NULL ||
            numArguments < 0 || numArguments > 27 ||
            trait->getDataSize() != sizeof(int2) || trait->getKeySize() != sizeof(int) ||
            string(trait->getDataType()) != "int2" || string(trait->getKeyType()) != "int" ||
            string(trait->getSortKey()) != "value.y")
        return false;
    for (int i = 0; i < numArguments; i++)
        if (arguments[i] == NULL)
            return false;
    sortImpl(data, false, generator, arguments, numArguments, true);
    return true;
}

void CudaSort::sortImpl(CudaArray& data, bool sortWithinBuckets, CUfunction generator,
        void* const* arguments, int numArguments) {
    // Preserve the existing exported overload and its complete data/key contract.
    sortImpl(data, sortWithinBuckets, generator, arguments, numArguments, false);
}

void CudaSort::sortImpl(CudaArray& data, bool sortWithinBuckets, CUfunction generator,
        void* const* arguments, int numArguments, bool directPhysicalIndices, ComputeKernel commonGenerator) {
    if (data.getSize() != dataLength || data.getElementSize() != trait->getDataSize())
        throw OpenMMException("CudaSort called with different data size");
    if (data.getSize() == 0)
        return;
    if (executionFlag != NULL) {
        sortConditional(data);
        return;
    }
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
        else if (generator != NULL) {
            void* generatedArgs[32];
            std::copy(arguments, arguments+numArguments, generatedArgs);
            generatedArgs[numArguments] = &numBuckets;
            generatedArgs[numArguments+1] = &dataRange.getDevicePointer();
            generatedArgs[numArguments+2] = &bucketOffset.getDevicePointer();
            generatedArgs[numArguments+3] = &bucketOfElement.getDevicePointer();
            generatedArgs[numArguments+4] = &offsetInBucket.getDevicePointer();
            context.executeKernel(generator, generatedArgs, data.getSize(), 128);
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

void CudaSort::sortConditional(CudaArray& data) {
    // All five stages, including bucket reset, carry the same GPU-side guard.
    // No CPU readback is introduced and no unconditional short-list copy is used.
    unsigned int numBuckets = bucketOffset.getSize();
    void* rangeArgs[] = {&data.getDevicePointer(), &dataLength, &dataRange.getDevicePointer(), &numBuckets,
            &bucketOffset.getDevicePointer(), &executionFlag->getDevicePointer()};
    context.executeKernel(computeRangeKernel, rangeArgs, rangeKernelSize, rangeKernelSize, 2*rangeKernelSize*trait->getKeySize());
    void* elementsArgs[] = {&data.getDevicePointer(), &dataLength, &numBuckets, &dataRange.getDevicePointer(),
            &bucketOffset.getDevicePointer(), &bucketOfElement.getDevicePointer(), &offsetInBucket.getDevicePointer(),
            &executionFlag->getDevicePointer()};
    context.executeKernel(assignElementsKernel, elementsArgs, data.getSize(), 128);
    void* computeArgs[] = {&numBuckets, &bucketOffset.getDevicePointer(), &executionFlag->getDevicePointer()};
    context.executeKernel(computeBucketPositionsKernel, computeArgs, positionsKernelSize, positionsKernelSize, positionsKernelSize*sizeof(int));
    void* copyArgs[] = {&data.getDevicePointer(), &buckets.getDevicePointer(), &dataLength, &bucketOffset.getDevicePointer(),
            &bucketOfElement.getDevicePointer(), &offsetInBucket.getDevicePointer(), &executionFlag->getDevicePointer()};
    context.executeKernel(copyToBucketsKernel, copyArgs, data.getSize());
    void* sortArgs[] = {&data.getDevicePointer(), &buckets.getDevicePointer(), &numBuckets, &bucketOffset.getDevicePointer(),
            &executionFlag->getDevicePointer()};
    context.executeKernel(sortBucketsKernel, sortArgs, ((data.getSize()+sortKernelSize-1)/sortKernelSize)*sortKernelSize, sortKernelSize, sortKernelSize*trait->getDataSize());
}
