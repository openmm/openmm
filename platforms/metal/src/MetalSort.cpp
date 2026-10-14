/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2010-2026 Stanford University and the Authors.      *
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

#include "MetalSort.h"
#include "MetalKernel.h"
#include "MetalKernelSources.h"
#include <algorithm>
#include <map>

using namespace OpenMM;
using namespace std;

MetalSort::MetalSort(MetalContext& context, ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform) :
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
    ComputeProgram program = context.compileProgram(context.replaceStrings(MetalKernelSources::sort, replacements));
    shortListKernel = program->createKernel("sortShortList");
    shortList2Kernel = program->createKernel("sortShortList2");
    computeRangeKernel = program->createKernel("computeRange");
    assignElementsKernel = program->createKernel(uniform ? "assignElementsToBuckets" : "assignElementsToBuckets2");
    computeBucketPositionsKernel = program->createKernel("computeBucketPositions");
    copyToBucketsKernel = program->createKernel("copyDataToBuckets");
    sortBucketsKernel = program->createKernel("sortBuckets");
    for (int i = 0; i < 2; i++)
        shortListKernel->addArg();
    for (int i = 0; i < 3; i++)
        shortList2Kernel->addArg();
    for (int i = 0; i < 5; i++)
        computeRangeKernel->addArg();
    for (int i = 0; i < 7; i++)
        assignElementsKernel->addArg();
    for (int i = 0; i < 2; i++)
        computeBucketPositionsKernel->addArg();
    for (int i = 0; i < 6; i++)
        copyToBucketsKernel->addArg();
    for (int i = 0; i < 4; i++)
        sortBucketsKernel->addArg();

    // Work out the work group sizes for various kernels.

    int maxBlockSize = 1024;
    int maxSharedMem = 32768;
    int maxLocalBuffer = (maxSharedMem/trait->getDataSize())/2;
    int maxShortList = min(3000, max(maxLocalBuffer, MetalContext::ThreadBlockSize*context.getNumThreadBlocks()));
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
        bucketOffset.initialize<int>(context, numBuckets, "bucketOffset");
        bucketOfElement.initialize<int>(context, length, "bucketOfElement");
        offsetInBucket.initialize<int>(context, length, "offsetInBucket");
    }
    buckets.initialize(context, length, trait->getDataSize(), "buckets");
}

MetalSort::~MetalSort() {
    delete trait;
}

void MetalSort::sort(ArrayInterface& data) {
    if (data.getSize() != dataLength || data.getElementSize() != trait->getDataSize())
        throw OpenMMException("MetalSort called with different data size");
    if (data.getSize() == 0)
        return;
    if (isShortList) {
        // We can use a simpler sort kernel that does the entire operation in one kernel.

        if (dataLength <= MetalContext::ThreadBlockSize*context.getNumThreadBlocks()) {
            shortList2Kernel->setArg(0, data);
            shortList2Kernel->setArg(1, buckets);
            shortList2Kernel->setArg(2, dataLength);
            shortList2Kernel->execute(dataLength);
            buckets.copyTo(data);
        }
        else {
            shortListKernel->setArg(0, data);
            shortListKernel->setArg(1, dataLength);
            dynamic_cast<MetalKernel*>(shortListKernel.get())->setDynamicLocalMemory(dataLength*trait->getDataSize());
            shortListKernel->execute(dataLength);
        }
    }
    else {
        // Compute the range of data values.

        unsigned int numBuckets = bucketOffset.getSize();
        computeRangeKernel->setArg(0, data);
        computeRangeKernel->setArg(1, dataLength);
        computeRangeKernel->setArg(2, dataRange);
        computeRangeKernel->setArg(3, numBuckets);
        computeRangeKernel->setArg(4, bucketOffset);
        dynamic_cast<MetalKernel*>(computeRangeKernel.get())->setDynamicLocalMemory(2*rangeKernelSize*trait->getKeySize());
        computeRangeKernel->execute(rangeKernelSize, rangeKernelSize);

        // Assign array elements to buckets.

        assignElementsKernel->setArg(0, data);
        assignElementsKernel->setArg(1, dataLength);
        assignElementsKernel->setArg(2, numBuckets);
        assignElementsKernel->setArg(3, dataRange);
        assignElementsKernel->setArg(4, bucketOffset);
        assignElementsKernel->setArg(5, bucketOfElement);
        assignElementsKernel->setArg(6, offsetInBucket);
        assignElementsKernel->execute(data.getSize(), 128);

        // Compute the position of each bucket.

        computeBucketPositionsKernel->setArg(0, numBuckets);
        computeBucketPositionsKernel->setArg(1, bucketOffset);
        dynamic_cast<MetalKernel*>(computeBucketPositionsKernel.get())->setDynamicLocalMemory(positionsKernelSize*sizeof(int));
        computeBucketPositionsKernel->execute(positionsKernelSize, positionsKernelSize);

        // Copy the data into the buckets.

        copyToBucketsKernel->setArg(0, data);
        copyToBucketsKernel->setArg(1, buckets);
        copyToBucketsKernel->setArg(2, dataLength);
        copyToBucketsKernel->setArg(3, bucketOffset);
        copyToBucketsKernel->setArg(4, bucketOfElement);
        copyToBucketsKernel->setArg(5, offsetInBucket);
        copyToBucketsKernel->execute(data.getSize());

        // Sort each bucket.

        sortBucketsKernel->setArg(0, data);
        sortBucketsKernel->setArg(1, buckets);
        sortBucketsKernel->setArg(2, numBuckets);
        sortBucketsKernel->setArg(3, bucketOffset);
        dynamic_cast<MetalKernel*>(sortBucketsKernel.get())->setDynamicLocalMemory(sortKernelSize*trait->getDataSize());
        sortBucketsKernel->execute(((data.getSize()+sortKernelSize-1)/sortKernelSize)*sortKernelSize, sortKernelSize);
    }
}
