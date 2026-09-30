/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform.
 * Source: platforms/opencl/include/OpenCLSort.h
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

#ifndef OPENMM_METALSORT_H_
#define OPENMM_METALSORT_H_

#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeKernel.h"
#include <memory>
#include "MetalContext.h"
#include "openmm/common/ComputeSort.h"
#include "openmm/common/windowsExportCommon.h"

namespace OpenMM {

/**
 * This class sorts arrays of values.  It supports any type of values, not just scalars,
 * so long as an appropriate sorting key can be defined by which to sort them.
 *
 * The sorting behavior is specified by a "trait" class that defines the type of data to
 * sort and the key for sorting it.  Here is an example of a trait class for
 * sorting floats:
 *
 * class FloatTrait : public ComputeSortImpl::SortTrait {
 *     int getDataSize() const {return 4;}
 *     int getKeySize() const {return 4;}
 *     const char* getDataType() const {return "float";}
 *     const char* getKeyType() const {return "float";}
 *     const char* getMinKey() const {return "-MAXFLOAT";}
 *     const char* getMaxKey() const {return "MAXFLOAT";}
 *     const char* getMaxValue() const {return "MAXFLOAT";}
 *     const char* getSortKey() const {return "value";}
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

class OPENMM_EXPORT_COMMON MetalSort : public ComputeSortImpl {
public:
    /**
     * Create an MetalSort object for sorting data of a particular type.
     *
     * @param context    the context in which to perform calculations
     * @param trait      a SortTrait defining the type of data to sort.  It should have been allocated
     *                   on the heap with the "new" operator.  This object takes over ownership of it,
     *                   and deletes it when the MetalSort is deleted.
     * @param length     the length of the arrays this object will be used to sort
     * @param uniform    whether the input data is expected to follow a uniform or nonuniform
     *                   distribution.  This argument is used only as a hint.  It allows parts
     *                   of the algorithm to be tuned for faster performance on the expected
     *                   distribution.
     */
    MetalSort(MetalContext& context, ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform=true);
    /** @brief Release the owned sort trait, kernels, and workspace arrays. */
    ~MetalSort();
    /**
     * Sort an array.
     */
    void sort(ArrayInterface& data);
private:
    MetalContext& context;
    std::unique_ptr<ComputeSortImpl::SortTrait> trait;
    ComputeArray dataRange;
    ComputeArray bucketOfElement;
    ComputeArray offsetInBucket;
    ComputeArray bucketOffset;
    ComputeArray buckets;
    ComputeKernel shortListKernel, computeRangeKernel, assignElementsKernel, computeBucketPositionsKernel, copyToBucketsKernel, sortBucketsKernel;
    unsigned int dataLength, rangeKernelSize, positionsKernelSize, sortKernelSize;
    bool isShortList, uniform;
};

} // namespace OpenMM

#endif // OPENMM_METALSORT_H_
