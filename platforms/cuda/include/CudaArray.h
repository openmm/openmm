#ifndef OPENMM_CUDAARRAY_H_
#define OPENMM_CUDAARRAY_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2012 Stanford University and the Authors.      *
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

#include "openmm/OpenMMException.h"
#include "windowsExportCuda.h"
#include <cuda.h>
#include <iostream>
#include <sstream>
#include <vector>

namespace OpenMM {

class CudaContext;

/**
 * This class encapsulates a block of CUDA device memory.  It provides a simplified API
 * for working with it and for copying data to and from device memory.
 */

class OPENMM_EXPORT_CUDA CudaArray {
public:
    /**
     * Create a CudaArray object.  The object is allocated on the heap with the "new" operator.
     * The template argument is the data type of each array element.
     *
     * @param context           the context for which to create the array
     * @param size              the number of elements in the array
     * @param name              the name of the array
     */
    template <class T>
    static CudaArray* create(CudaContext& context, int size, const std::string& name) {
        return new CudaArray(context, size, sizeof(T), name);
    }
    /**
     * Create a CudaArray object.
     *
     * @param context           the context for which to create the array
     * @param size              the number of elements in the array
     * @param elementSize       the size of each element in bytes
     * @param name              the name of the array
     */
    CudaArray(CudaContext& context, int size, int elementSize, const std::string& name);
    ~CudaArray();
    /**
     * Get the number of elements in the array.
     */
    int getSize() const {
        return size;
    }
    /**
     * Get the size of each element in bytes.
     */
    int getElementSize() const {
        return elementSize;
    }
    /**
     * Get the name of the array.
     */
    const std::string& getName() const {
        return name;
    }
    /**
     * Get a pointer to the device memory.
     */
    CUdeviceptr& getDevicePointer() {
        return pointer;
    }
    /**
     * Copy the values in a vector to the device memory.
     */
    template <class T>
    void upload(const std::vector<T>& data) {
        if (sizeof(T) != elementSize || data.size() != size)
            throw OpenMMException("Error uploading array "+name+": The specified vector does not match the size of the array");
        upload(&data[0], true);
    }
    /**
     * Copy the values in the Buffer to a vector.
     */
    template <class T>
    void download(std::vector<T>& data) const {
        if (sizeof(T) != elementSize)
            throw OpenMMException("Error downloading array "+name+": The specified vector has the wrong element size");
        if (data.size() != size)
            data.resize(size);
        download(&data[0], true);
    }
    /**
     * Copy the values in an array to the device memory.
     * 
     * @param data     the data to copy
     * @param blocking if true, this call will block until the transfer is complete.  If false,
     *                 the source array  must be in page-locked memory.
     */
    void upload(const void* data, bool blocking = true);
    /**
     * Copy the values in the device memory to an array.
     * 
     * @param data     the array to copy the memory to
     * @param blocking if true, this call will block until the transfer is complete.  If false,
     *                 the destination array must be in page-locked memory.
     */
    void download(void* data, bool blocking = true) const;
    /**
     * Copy the values in the device memory to a second array.
     * 
     * @param dest     the destination array to copy to
     */
    void copyTo(CudaArray& dest) const;
private:
    CudaContext& context;
    CUdeviceptr pointer;
    int size, elementSize;
    bool ownsMemory;
    std::string name;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAARRAY_H_*/
