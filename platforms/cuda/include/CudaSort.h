#ifndef __OPENMM_CUDASORT_H__
#define __OPENMM_CUDASORT_H__

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

#include "CudaArray.h"
#include "openmm/common/ComputeSort.h"
#include "openmm/common/windowsExportCommon.h"
#include "CudaContext.h"

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
 *     const char* getMinKey() const {return "-3.40282e+38f";}
 *     const char* getMaxKey() const {return "3.40282e+38f";}
 *     const char* getMaxValue() const {return "3.40282e+38f";}
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
    
class OPENMM_EXPORT_COMMON CudaSort : public ComputeSortImpl {
public:
    /**
     * Create a CudaSort object for sorting data of a particular type.
     *
     * @param context    the context in which to perform calculations
     * @param trait      a SortTrait defining the type of data to sort.  It should have been allocated
     *                   on the heap with the "new" operator.  This object takes over ownership of it,
     *                   and deletes it when the CudaSort is deleted.
     * @param length     the length of the arrays this object will be used to sort
     * @param uniform    whether the input data is expected to follow a uniform or nonuniform
     *                   distribution.  This argument is used only as a hint.  It allows parts
     *                   of the algorithm to be tuned for faster performance on the expected
     *                   distribution.
     */
    CudaSort(CudaContext& context, ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform=true);
    /**
     * Experimental conditional long-list sort.  Every stage reads the same immutable
     * device int flag.  If zero, all stages leave input and scratch arrays untouched.
     * The flag must belong to this context and outlive the sorter.  It must be set
     * on the current stream before sort() and not modified until sort() completes.
     * Short lists are deliberately unsupported (their host copy would be unguarded).
     */
    CudaSort(CudaContext& context, ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform, CudaArray* executionFlag);
    /**
     * Experimental known-range integer sort.  A nonnegative maximum promises
     * every key is in [0, knownIntegerMaximum].  Only unguarded uniform int-key
     * sorting is supported; -1 retains the original range calculation.
     * The fixed range changes bucket/tie ordering, not key ordering.
     */
    CudaSort(CudaContext& context, ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform,
             CudaArray* executionFlag, int knownIntegerMaximum);
    ~CudaSort();
    /**
     * Sort an array.
     */
    void sort(ArrayInterface& data);
    /**
     * Experimental coarse permutation: partition into the same buckets without
     * sorting within a bucket.  The result is NOT guaranteed to be key-sorted.
     * Short lists and conditional sorts retain the complete sorting path.
     * Only callers that require permutation coverage, such as PME, may use it.
     */
    void bucketize(CudaArray& data);
    /**
     * Generate int2/value.y keys while assigning buckets. The generator receives
     * caller arguments followed by numBuckets, dataRange, bucketOffset,
     * bucketOfElement and offsetInBucket. It must initialize every data element.
     * Returns false without enqueuing work if the fixed-range long-list contract
     * is unavailable. Launch errors propagate; they must not trigger fallback.
     * sortWithinBuckets=false has the same permutation contract as bucketize().
     */
    bool sortWithGeneratedKeys(CudaArray& data, CUfunction generator, void* const* arguments,
            int numArguments, bool sortWithinBuckets);
    /**
     * Experimental physical-index permutation for the bounded PME coarse path.
     * The generator initializes bucketOfElement/offsetInBucket for physical
     * indices [0, length); it need not write data. The final scatter overwrites
     * data with int2(physicalIndex, 0), so keys and original values are discarded.
     * Only consumers of .x permutation coverage may use the result. Before any
     * later generic sort, the caller must regenerate every complete data/key pair.
     * Workspace arguments match sortWithGeneratedKeys. Unsupported contracts
     * return false before any enqueue; errors after enqueue propagate.
     */
    bool bucketizeGeneratedPhysicalIndices(CudaArray& data, CUfunction generator,
            void* const* arguments, int numArguments);
    /**
     * Generate bounded int2 keys with a ComputeKernel.  Arguments [0, 9] are
     * initialized by the caller and [10, 14] must be reserved for workspace.
     * A false result guarantees that no device work was enqueued.
     */
    bool sortWithGeneratedKeys(CudaArray& data, ComputeKernel generator,
            bool sortWithinBuckets, bool directPhysicalIndices);
private:
    void sortImpl(CudaArray& data, bool sortWithinBuckets);
    void sortImpl(CudaArray& data, bool sortWithinBuckets, CUfunction generator,
            void* const* arguments, int numArguments);
    void sortImpl(CudaArray& data, bool sortWithinBuckets, CUfunction generator,
            void* const* arguments, int numArguments, bool directPhysicalIndices, ComputeKernel commonGenerator=ComputeKernel());
    void sortConditional(CudaArray& data);
    CudaArray* executionFlag;
    CudaContext& context;
    ComputeSortImpl::SortTrait* trait;
    CudaArray dataRange;
    CudaArray bucketOfElement;
    CudaArray offsetInBucket;
    CudaArray bucketOffset;
    CudaArray buckets;
    CUfunction shortListKernel, shortList2Kernel, computeRangeKernel, assignElementsKernel, computeBucketPositionsKernel, copyToBucketsKernel, sortBucketsKernel;
    CUfunction scatterPhysicalIndicesKernel;
    unsigned int dataLength, rangeKernelSize, positionsKernelSize, sortKernelSize;
    bool isShortList, uniform, hasFixedRange;
};

} // namespace OpenMM

#endif // __OPENMM_CUDASORT_H__
