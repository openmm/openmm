/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2020 Stanford University and the Authors.      *
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

#ifdef _MSC_VER
    // Prevent Windows from defining macros that interfere with other code.
    #define NOMINMAX
#endif
#include "OpenCLSort.h"
#include "OpenCLKernelSources.h"
#include <algorithm>
#include <map>
#include <string>

using namespace OpenMM;
using namespace std;

OpenCLSort::OpenCLSort(OpenCLContext& context, SortTrait* trait, unsigned int length) : context(context), trait(trait), dataLength(length) {
    // Create kernels.

    std::map<std::string, std::string> replacements;
    replacements["DATA_TYPE"] = trait->getDataType();
    replacements["KEY_TYPE"] =  trait->getKeyType();
    replacements["SORT_KEY"] = trait->getSortKey();
    replacements["MIN_KEY"] = trait->getMinKey();
    replacements["MAX_KEY"] = trait->getMaxKey();
    replacements["MAX_VALUE"] = trait->getMaxValue();
    cl::Program program = context.createProgram(context.replaceStrings(OpenCLKernelSources::sort, replacements));
    shortListKernel = cl::Kernel(program, "sortShortList");
    shortList2Kernel = cl::Kernel(program, "sortShortList2");
    computeRangeKernel = cl::Kernel(program, "computeRange");
    assignElementsKernel = cl::Kernel(program, "assignElementsToBuckets");
    computeBucketPositionsKernel = cl::Kernel(program, "computeBucketPositions");
    copyToBucketsKernel = cl::Kernel(program, "copyDataToBuckets");
    sortBucketsKernel = cl::Kernel(program, "sortBuckets");

    // Work out the work group sizes for various kernels.

    unsigned int maxGroupSize = std::min(256, (int) context.getDevice().getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
    int maxSharedMem = context.getDevice().getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
    unsigned int maxRangeSize = std::min(maxGroupSize, (unsigned int) computeRangeKernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(context.getDevice()));
    unsigned int maxPositionsSize = std::min(maxGroupSize, (unsigned int) computeBucketPositionsKernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(context.getDevice()));
    int maxLocalBuffer = (maxSharedMem/trait->getDataSize())/2;
    unsigned int maxShortList = min(8192, max(maxLocalBuffer, (int) OpenCLContext::ThreadBlockSize*context.getNumThreadBlocks()));
    // The following line checks CL_KERNEL_WORK_GROUP_SIZE to make sure we don't create too large a workgroup.
    // Unfortunately, AMD's OpenCL returns an inappropriately small value for it that is much shorter than the actual
    // maximum, so including the check hurts performance.  For the moment I'm just leaving it commented out.
    // If the workgroup size turns out to be too large, we catch the exception and switch back to the standard
    // sorting kernels.
    //maxShortList = min(maxShortList, shortListKernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(context.getDevice()));
    isShortList = (length <= maxShortList);
    string vendor = context.getDevice().getInfo<CL_DEVICE_VENDOR>();
    if (vendor.size() >= 6 && vendor.substr(0, 6) == "NVIDIA")
        useShortList2 = (dataLength <= OpenCLContext::ThreadBlockSize*context.getNumThreadBlocks());
    else
        useShortList2 = false;
    for (rangeKernelSize = 1; rangeKernelSize*2 <= maxRangeSize; rangeKernelSize *= 2)
        ;
    positionsKernelSize = std::min(rangeKernelSize, maxPositionsSize);
    sortKernelSize = (isShortList ? rangeKernelSize : rangeKernelSize/2);
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

    dataRange.initialize(context, 2, trait->getKeySize(), "sortDataRange");
    bucketOffset.initialize<cl_uint>(context, numBuckets, "bucketOffset");
    bucketOfElement.initialize<cl_uint>(context, length, "bucketOfElement");
    offsetInBucket.initialize<cl_uint>(context, length, "offsetInBucket");
    buckets.initialize(context, length, trait->getDataSize(), "buckets");
}

OpenCLSort::~OpenCLSort() {
    delete trait;
}

void OpenCLSort::sort(OpenCLArray& data) {
    if (data.getSize() != dataLength || data.getElementSize() != trait->getDataSize())
        throw OpenMMException("OpenCLSort called with different data size");
    if (data.getSize() == 0)
        return;
    if (isShortList) {
        // We can use a simpler sort kernel that does the entire operation in one kernel.
        
        try {
            if (useShortList2) {
                shortList2Kernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
                shortList2Kernel.setArg<cl::Buffer>(1, buckets.getDeviceBuffer());
                shortList2Kernel.setArg<cl_int>(2, dataLength);
                context.executeKernel(shortList2Kernel, dataLength);
                buckets.copyTo(data);
            }
            else {
                shortListKernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
                shortListKernel.setArg<cl_uint>(1, dataLength);
                shortListKernel.setArg(2, dataLength*trait->getDataSize(), NULL);
                context.executeKernel(shortListKernel, sortKernelSize, sortKernelSize);
            }
            return;
        }
        catch (exception& ex) {
            // This can happen if we chose too large a size for the kernel.  Switch
            // over to the standard sorting method.
            
            isShortList = false;
        }
    }

    // Compute the range of data values.

    unsigned int numBuckets = bucketOffset.getSize();
    computeRangeKernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
    computeRangeKernel.setArg<cl_uint>(1, data.getSize());
    computeRangeKernel.setArg<cl::Buffer>(2, dataRange.getDeviceBuffer());
    computeRangeKernel.setArg(3, rangeKernelSize*trait->getKeySize(), NULL);
    computeRangeKernel.setArg(4, rangeKernelSize*trait->getKeySize(), NULL);
    computeRangeKernel.setArg<cl_int>(5, numBuckets);
    computeRangeKernel.setArg<cl::Buffer>(6, bucketOffset.getDeviceBuffer());
    context.executeKernel(computeRangeKernel, rangeKernelSize, rangeKernelSize);

    // Assign array elements to buckets.

    assignElementsKernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
    assignElementsKernel.setArg<cl_int>(1, data.getSize());
    assignElementsKernel.setArg<cl_int>(2, numBuckets);
    assignElementsKernel.setArg<cl::Buffer>(3, dataRange.getDeviceBuffer());
    assignElementsKernel.setArg<cl::Buffer>(4, bucketOffset.getDeviceBuffer());
    assignElementsKernel.setArg<cl::Buffer>(5, bucketOfElement.getDeviceBuffer());
    assignElementsKernel.setArg<cl::Buffer>(6, offsetInBucket.getDeviceBuffer());
    context.executeKernel(assignElementsKernel, data.getSize());

    // Compute the position of each bucket.

    computeBucketPositionsKernel.setArg<cl_int>(0, numBuckets);
    computeBucketPositionsKernel.setArg<cl::Buffer>(1, bucketOffset.getDeviceBuffer());
    computeBucketPositionsKernel.setArg(2, positionsKernelSize*sizeof(cl_int), NULL);
    context.executeKernel(computeBucketPositionsKernel, positionsKernelSize, positionsKernelSize);

    // Copy the data into the buckets.

    copyToBucketsKernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
    copyToBucketsKernel.setArg<cl::Buffer>(1, buckets.getDeviceBuffer());
    copyToBucketsKernel.setArg<cl_int>(2, data.getSize());
    copyToBucketsKernel.setArg<cl::Buffer>(3, bucketOffset.getDeviceBuffer());
    copyToBucketsKernel.setArg<cl::Buffer>(4, bucketOfElement.getDeviceBuffer());
    copyToBucketsKernel.setArg<cl::Buffer>(5, offsetInBucket.getDeviceBuffer());
    context.executeKernel(copyToBucketsKernel, data.getSize());

    // Sort each bucket.

    sortBucketsKernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
    sortBucketsKernel.setArg<cl::Buffer>(1, buckets.getDeviceBuffer());
    sortBucketsKernel.setArg<cl_int>(2, numBuckets);
    sortBucketsKernel.setArg<cl::Buffer>(3, bucketOffset.getDeviceBuffer());
    sortBucketsKernel.setArg(4, sortKernelSize*trait->getDataSize(), NULL);
    context.executeKernel(sortBucketsKernel, ((data.getSize()+sortKernelSize-1)/sortKernelSize)*sortKernelSize, sortKernelSize);
}
