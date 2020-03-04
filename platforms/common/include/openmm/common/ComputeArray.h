#ifndef OPENMM_COMPUTEARRAY_H_
#define OPENMM_COMPUTEARRAY_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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

namespace OpenMM {

/**
 * This is an implementation of ArrayInterface that acts as a wrapper around a platform-specific
 * array implementation (typically CudaArray or OpenCLArray).  This class can be used in code that
 * is not platform-specific, and an appropriate implementation array is created automatically
 * based on the ComputeContext.
 */

class OPENMM_EXPORT_COMMON ComputeArray : public ArrayInterface {
public:
    /**
     * Create an uninitialized ComputeArray object.  It cannot be used until initialize() is called on it.
     */
    ComputeArray();
    /**
     * Release all resources allocated by this object.
     */
    ~ComputeArray();
    /**
     * Get the internal array this object is wrapping.
     */
    ArrayInterface& getArray();
    /**
     * Initialize this array.
     *
     * @param context           the context for which to create the array
     * @param size              the number of elements in the array
     * @param elementSize       the size of each element in bytes
     * @param name              the name of the array
     */
    void initialize(ComputeContext& context, int size, int elementSize, const std::string& name);
    /**
     * Initialize this object.  The template argument is the data type of each array element.
     *
     * @param context           the context for which to create the array
     * @param size              the number of elements in the array
     * @param name              the name of the array
     */
    template <class T>
    void initialize(ComputeContext& context, int size, const std::string& name) {
        initialize(context, size, sizeof(T), name);
    }
    /**
     * Recreate the internal storage to have a different size.
     */
    void resize(int size);
    /**
     * Get whether this array has been initialized.
     */
    bool isInitialized() const;
    /**
     * Get the number of elements in the array.
     */
    int getSize() const;
    /**
     * Get the size of each element in bytes.
     */
    int getElementSize() const;
    /**
     * Get the name of the array.
     */
    const std::string& getName() const;
    /**
     * Get the context this array belongs to.
     */
    ComputeContext& getContext();
    /**
     * Copy the values in a vector to the Buffer.
     */
    template <class T>
    void upload(const std::vector<T>& data, bool convert=false) {
        ArrayInterface::upload(data, convert);
    }
    /**
     * Copy the values in the Buffer to a vector.
     */
    template <class T>
    void download(std::vector<T>& data) const {
        ArrayInterface::download(data);
    }
    /**
     * Copy the values from host memory to the array.
     * 
     * @param data     the data to copy
     * @param blocking if true, this call will block until the transfer is complete.  Subclasses often
     *                 have restrictions on non-blocking copies, such as that the source data must be
     *                 in page-locked memory.
     */
    void upload(const void* data, bool blocking=true);
    /**
     * Copy the values in the array to host memory.
     * 
     * @param data     the destination to copy the value to
     * @param blocking if true, this call will block until the transfer is complete.  Subclasses often
     *                 have restrictions on non-blocking copies, such as that the destination must be
     *                 in page-locked memory.
     */
    void download(void* data, bool blocking=true) const;
    /**
     * Copy the values in this array to a second array.
     * 
     * @param dest     the destination array to copy to
     */
    void copyTo(ArrayInterface& dest) const;
private:
    ArrayInterface* impl;
};

} // namespace OpenMM

#endif /*OPENMM_COMPUTEARRAY_H_*/
