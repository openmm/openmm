#ifndef OPENMM_COMPUTESORT_H_
#define OPENMM_COMPUTESORT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "openmm/common/ArrayInterface.h"
#include "openmm/common/windowsExportCommon.h"
#include <memory>

namespace OpenMM {

/**
 * This abstract class represents an algorithm for sorting arrays.  It is created
 * by calling createEvent() on a ComputeContext, which returns an instance of a
 * platform-specific subclass.
 *
 * Instead of referring to this class directly, it is best to use a ComputeSort, which is
 * a typedef for a shared_ptr to a ComputeSortImpl.  This allows you to treat it as having
 * value semantics, and frees you from having to manage memory.  
 */

class OPENMM_EXPORT_COMMON ComputeSortImpl {
public:
    class SortTrait;
    virtual ~ComputeSortImpl() {
    }
    /**
     * Sort an array.
     */
    virtual void sort(ArrayInterface& data) = 0;
};

typedef std::shared_ptr<ComputeSortImpl> ComputeSort;

/**
 * A subclass of SortTrait defines the type of value to sort, and the key for sorting them.
 */
class ComputeSortImpl::SortTrait {
public:
    virtual ~SortTrait() {
    }
    /**
     * Get the size of each data value in bytes.
     */
    virtual int getDataSize() const = 0;
    /**
     * Get the size of each key value in bytes.
     */
    virtual int getKeySize() const = 0;
    /**
     * Get the data type of the values to sort.
     */
    virtual const char* getDataType() const = 0;
    /**
     * Get the data type of the sorting key.
     */
    virtual const char* getKeyType() const = 0;
    /**
     * Get the minimum value a key can take.
     */
    virtual const char* getMinKey() const = 0;
    /**
     * Get the maximum value a key can take.
     */
    virtual const char* getMaxKey() const = 0;
    /**
     * Get a value whose key is guaranteed to equal getMaxKey().
     */
    virtual const char* getMaxValue() const = 0;
    /**
     * Get the source code to select the key from the data value.
     */
    virtual const char* getSortKey() const = 0;
};

} // namespace OpenMM

#endif /*OPENMM_COMPUTESORT_H_*/
