#ifndef __OPENMM_OPENCLSORT_H__
#define __OPENMM_OPENCLSORT_H__

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "OpenCLArray.h"
#include "OpenCLKernelSources.h"
#include "openmm/internal/windowsExport.h"
#include <map>

namespace OpenMM {

/**
 * This class sorts arrays of values.  It supports any type of values, not just scalars,
 * so long as an appropriate sorting key can be defined by which to sort them.
 *
 * The algorithm used is a bucket sort, followed by a bitonic sort within each bucket
 * (in local memory when possible, in global memory otherwise).  This is similar to
 * the algorithm described in
 *
 * Shifu Chen, Jing Qin, Yongming Xie, Junping Zhao, and Pheng-Ann Heng.  "An Efficient
 * Sorting Algorithm with CUDA"  Journal of the Chinese Institute of Engineers, 32(7),
 * pp. 915-921 (2009)
 *
 * but with many modifications and simplifications.  In particular, this algorithm
 * involves much less communication between host and device, which is critical to get
 * good performance with the array sizes we typically work with (10,000 to 100,000
 * elements).
 */

template <class TYPE>
class OPENMM_EXPORT OpenCLSort {
public:
    /**
     * Create an OpenCLSort object for sorting data of a particular type.
     *
     * @param context    the context in which to perform calculations
     * @param length     the length of the arrays this object will be used to sort
     * @param typeName   the name of the data type being sorting (e.g. "float")
     * @param sortKey    an expression that returns the value by which the variable "value" should be sorted.
     *                   For primitive types, this will simply be "value".
     */
    OpenCLSort(OpenCLContext& context, int length, const std::string& typeName, const std::string& sortKey) : context(context),
            dataRange(NULL), bucketOfElement(NULL), offsetInBucket(NULL), bucketOffset(NULL), buckets(NULL) {
        // Create kernels.

        std::map<std::string, std::string> replacements;
        replacements["TYPE"] = typeName;
        replacements["SORT_KEY"] = sortKey;
        cl::Program program = context.createProgram(context.replaceStrings(OpenCLKernelSources::sort, replacements));
        computeRangeKernel = cl::Kernel(program, "computeRange");
        assignElementsKernel = cl::Kernel(program, "assignElementsToBuckets");
        computeBucketPositionsKernel = cl::Kernel(program, "computeBucketPositions");
        copyToBucketsKernel = cl::Kernel(program, "copyDataToBuckets");
        sortBucketsKernel = cl::Kernel(program, "sortBuckets");

        // Work out the work group sizes for various kernels.

        int maxGroupSize = std::min(256, (int) context.getDevice().getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
        for (rangeKernelSize = 1; rangeKernelSize*2 <= maxGroupSize; rangeKernelSize *= 2)
            ;
        positionsKernelSize = rangeKernelSize;
        sortKernelSize = rangeKernelSize/2;
        if (rangeKernelSize > length)
            rangeKernelSize = length;
        int maxLocalBuffer = (int) ((context.getDevice().getInfo<CL_DEVICE_LOCAL_MEM_SIZE>()/sizeof(TYPE))/2);
        if (sortKernelSize > maxLocalBuffer)
            sortKernelSize = maxLocalBuffer;
        int targetBucketSize = sortKernelSize/2;
        int numBuckets = length/targetBucketSize;
        if (numBuckets < 1)
            numBuckets = 1;
        if (positionsKernelSize > numBuckets)
            positionsKernelSize = numBuckets;

        // Create workspace arrays.

        dataRange = new OpenCLArray<mm_float2>(context, 1, "sortDataRange");
        bucketOffset = new OpenCLArray<cl_int>(context, numBuckets, "bucketOffset");
        bucketOfElement = new OpenCLArray<cl_int>(context, length, "bucketOfElement");
        offsetInBucket = new OpenCLArray<cl_int>(context, length, "offsetInBucket");
        buckets = new OpenCLArray<TYPE>(context, length, "buckets");
    }
    ~OpenCLSort() {
        if (dataRange != NULL)
            delete dataRange;
        if (bucketOfElement != NULL)
            delete bucketOfElement;
        if (offsetInBucket != NULL)
            delete offsetInBucket;
        if (bucketOffset != NULL)
            delete bucketOffset;
        if (buckets != NULL)
            delete buckets;
    }
    /**
     * Sort an array.
     */
    void sort(OpenCLArray<TYPE>& data) {
        // Compute the range of data values.

        computeRangeKernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
        computeRangeKernel.setArg<cl_int>(1, data.getSize());
        computeRangeKernel.setArg<cl::Buffer>(2, dataRange->getDeviceBuffer());
        computeRangeKernel.setArg(3, rangeKernelSize*sizeof(cl_float), NULL);
        context.executeKernel(computeRangeKernel, rangeKernelSize, rangeKernelSize);

        // Assign array elements to buckets.

        int numBuckets = bucketOffset->getSize();
        context.clearBuffer(bucketOffset->getDeviceBuffer(), numBuckets);
        assignElementsKernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
        assignElementsKernel.setArg<cl_int>(1, data.getSize());
        assignElementsKernel.setArg<cl_int>(2, numBuckets);
        assignElementsKernel.setArg<cl::Buffer>(3, dataRange->getDeviceBuffer());
        assignElementsKernel.setArg<cl::Buffer>(4, bucketOffset->getDeviceBuffer());
        assignElementsKernel.setArg<cl::Buffer>(5, bucketOfElement->getDeviceBuffer());
        assignElementsKernel.setArg<cl::Buffer>(6, offsetInBucket->getDeviceBuffer());
        context.executeKernel(assignElementsKernel, data.getSize());

        // Compute the position of each bucket.

        computeBucketPositionsKernel.setArg<cl_int>(0, numBuckets);
        computeBucketPositionsKernel.setArg<cl::Buffer>(1, bucketOffset->getDeviceBuffer());
        computeBucketPositionsKernel.setArg(2, positionsKernelSize*sizeof(cl_int), NULL);
        context.executeKernel(computeBucketPositionsKernel, positionsKernelSize, positionsKernelSize);

        // Copy the data into the buckets.

        copyToBucketsKernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
        copyToBucketsKernel.setArg<cl::Buffer>(1, buckets->getDeviceBuffer());
        copyToBucketsKernel.setArg<cl_int>(2, data.getSize());
        copyToBucketsKernel.setArg<cl::Buffer>(3, bucketOffset->getDeviceBuffer());
        copyToBucketsKernel.setArg<cl::Buffer>(4, bucketOfElement->getDeviceBuffer());
        copyToBucketsKernel.setArg<cl::Buffer>(5, offsetInBucket->getDeviceBuffer());
        context.executeKernel(copyToBucketsKernel, data.getSize());

        // Sort each bucket.

        sortBucketsKernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
        sortBucketsKernel.setArg<cl::Buffer>(1, buckets->getDeviceBuffer());
        sortBucketsKernel.setArg<cl_int>(2, numBuckets);
        sortBucketsKernel.setArg<cl::Buffer>(3, bucketOffset->getDeviceBuffer());
        sortBucketsKernel.setArg(4, sortKernelSize*sizeof(TYPE), NULL);
        context.executeKernel(sortBucketsKernel, ((data.getSize()+sortKernelSize-1)/sortKernelSize)*sortKernelSize, sortKernelSize);
    }
private:
    OpenCLContext& context;
    OpenCLArray<mm_float2>* dataRange;
    OpenCLArray<cl_int>* bucketOfElement;
    OpenCLArray<cl_int>* offsetInBucket;
    OpenCLArray<cl_int>* bucketOffset;
    OpenCLArray<TYPE>* buckets;
    cl::Kernel computeRangeKernel, assignElementsKernel, computeBucketPositionsKernel, copyToBucketsKernel, sortBucketsKernel;
    int rangeKernelSize, positionsKernelSize, sortKernelSize;
};

} // namespace OpenMM

#endif // __OPENMM_OPENCLSORT_H__
