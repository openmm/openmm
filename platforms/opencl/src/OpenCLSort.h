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
 * The class is templatized by a "trait" class that defines the type of data to
 * sort and the key for sorting it.  Here is an example of a trait class for
 * sorting floats:
 * 
 * struct FloatTrait {
 *     // The name of the data and key types being sorted.
 *     // Both the host type and OpenCL type is required.
 *     // For primitive types they will be the same.
 *     typedef cl_float DataType;
 *     typedef cl_float KeyType;
 *     static const char* clDataType() {return "float";}
 *     static const char* clKeyType() {return "float";}
 *     // The minimum value a key can take.
 *     static const char* clMinKey() {return "-MAXFLOAT";}
 *     // The maximum value a key can take.
 *     static const char* clMaxKey() {return "MAXFLOAT";}
 *     // A value whose key is guaranteed to equal clMaxKey().
 *     static const char* clMaxValue() {return "MAXFLOAT";}
 *     // The OpenCL code to select the key from the data value.
 *     static const char* clSortKey() {return "value";}
 * };
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
    
template <class TRAIT>
class OpenCLSort {
public:
    /**
     * Create an OpenCLSort object for sorting data of a particular type.
     *
     * @param context    the context in which to perform calculations
     * @param length     the length of the arrays this object will be used to sort
     */
    OpenCLSort(OpenCLContext& context, unsigned int length) : context(context),
            dataRange(NULL), bucketOfElement(NULL), offsetInBucket(NULL), bucketOffset(NULL), buckets(NULL) {
        // Create kernels.

        std::map<std::string, std::string> replacements;
        replacements["DATA_TYPE"] = TRAIT::clDataType();
        replacements["KEY_TYPE"] =  TRAIT::clKeyType();
        replacements["SORT_KEY"] = TRAIT::clSortKey();
        replacements["MIN_KEY"] = TRAIT::clMinKey();
        replacements["MAX_KEY"] = TRAIT::clMaxKey();
        replacements["MAX_VALUE"] = TRAIT::clMaxValue();
        replacements["VALUE_IS_INT2"] = (TRAIT::clDataType() == "int2" ? "1" : "0");
        cl::Program program = context.createProgram(context.replaceStrings(OpenCLKernelSources::sort, replacements));
        computeRangeKernel = cl::Kernel(program, "computeRange");
        assignElementsKernel = cl::Kernel(program, "assignElementsToBuckets");
        computeBucketPositionsKernel = cl::Kernel(program, "computeBucketPositions");
        copyToBucketsKernel = cl::Kernel(program, "copyDataToBuckets");
        sortBucketsKernel = cl::Kernel(program, "sortBuckets");

        // Work out the work group sizes for various kernels.

        unsigned int maxGroupSize = std::min(256, (int) context.getDevice().getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
        for (rangeKernelSize = 1; rangeKernelSize*2 <= maxGroupSize; rangeKernelSize *= 2)
            ;
        positionsKernelSize = rangeKernelSize;
        sortKernelSize = rangeKernelSize/2;
        if (rangeKernelSize > length)
            rangeKernelSize = length;
        unsigned int maxLocalBuffer = (unsigned int) ((context.getDevice().getInfo<CL_DEVICE_LOCAL_MEM_SIZE>()/sizeof(typename TRAIT::DataType))/2);
        if (sortKernelSize > maxLocalBuffer)
            sortKernelSize = maxLocalBuffer;
        unsigned int targetBucketSize = sortKernelSize/2;
        unsigned int numBuckets = length/targetBucketSize;
        if (numBuckets < 1)
            numBuckets = 1;
        if (positionsKernelSize > numBuckets)
            positionsKernelSize = numBuckets;

        // Create workspace arrays.

        dataRange = new OpenCLArray<typename TRAIT::KeyType>(context, 2, "sortDataRange");
        bucketOffset = new OpenCLArray<cl_uint>(context, numBuckets, "bucketOffset");
        bucketOfElement = new OpenCLArray<cl_uint>(context, length, "bucketOfElement");
        offsetInBucket = new OpenCLArray<cl_uint>(context, length, "offsetInBucket");
        buckets = new OpenCLArray<typename TRAIT::DataType>(context, length, "buckets");
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
    void sort(OpenCLArray<typename TRAIT::DataType>& data) {

        if (data.getSize() != bucketOfElement->getSize())
            throw OpenMMException("OpenCLSort called with different data size");
        if (data.getSize() == 0)
            return;

        // Compute the range of data values.

        computeRangeKernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
        computeRangeKernel.setArg<cl_uint>(1, data.getSize());
        computeRangeKernel.setArg<cl::Buffer>(2, dataRange->getDeviceBuffer());
        computeRangeKernel.setArg(3, rangeKernelSize*sizeof(typename TRAIT::KeyType), NULL);
        context.executeKernel(computeRangeKernel, rangeKernelSize, rangeKernelSize);

        // Assign array elements to buckets.

        unsigned int numBuckets = bucketOffset->getSize();
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
        sortBucketsKernel.setArg(4, sortKernelSize*sizeof(typename TRAIT::DataType), NULL);
        context.executeKernel(sortBucketsKernel, ((data.getSize()+sortKernelSize-1)/sortKernelSize)*sortKernelSize, sortKernelSize);
    }
private:
    OpenCLContext& context;
    OpenCLArray<typename TRAIT::KeyType>* dataRange;
    OpenCLArray<cl_uint>* bucketOfElement;
    OpenCLArray<cl_uint>* offsetInBucket;
    OpenCLArray<cl_uint>* bucketOffset;
    OpenCLArray<typename TRAIT::DataType>* buckets;
    cl::Kernel computeRangeKernel, assignElementsKernel, computeBucketPositionsKernel, copyToBucketsKernel, sortBucketsKernel;
    unsigned int rangeKernelSize, positionsKernelSize, sortKernelSize;
};

} // namespace OpenMM

#endif // __OPENMM_OPENCLSORT_H__
