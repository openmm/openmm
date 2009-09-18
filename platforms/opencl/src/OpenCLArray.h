#ifndef OPENMM_OPENCLARRAY_H_
#define OPENMM_OPENCLARRAY_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

#include "OpenCLContext.h"
#include "openmm/OpenMMException.h"
#include <vector>

namespace OpenMM {

/**
 * This class encapsulates an OpenCL Buffer.  It provides a simplified API for working with it,
 * an optionally includes a buffer in host memory for copying data to and from the OpenCL Buffer.
 */

template <class T>
class OpenCLArray {
public:
    /**
     * Create an OpenCLArray object.
     *
     * @param context           the context for which to create the array
     * @param size              the number of elements in the array
     * @param name              the name of the array
     * @param createHostBuffer  specifies whether to create a buffer in host memory for copying data to and from
     *                          the OpenCL Buffer
     */
    OpenCLArray(OpenCLContext& context, int size, const std::string& name, bool createHostBuffer = false) :
            context(context), size(size), name(name), local(createHostBuffer ? size : 0), ownsBuffer(true) {
        buffer = new cl::Buffer(context.getContext(), CL_MEM_READ_WRITE, size*sizeof(T));
    }
    /**
     * Create an OpenCLArray object the uses a preexisting Buffer.
     *
     * @param context           the context for which to create the array
     * @param buffer            the OpenCL Buffer this object encapsulates
     * @param size              the number of elements in the array
     * @param name              the name of the array
     * @param createHostBuffer  specifies whether to create a buffer in host memory for copying data to and from
     *                          the OpenCL Buffer
     */
    OpenCLArray(OpenCLContext& context, cl::Buffer* buffer, int size, const std::string& name, bool createHostBuffer = false) :
            context(context), buffer(buffer), size(size), name(name), local(createHostBuffer ? size : 0), ownsBuffer(false) {
    }
    ~OpenCLArray() {
        if (ownsBuffer)
            delete buffer;
    }
    const T& operator[](int index) const {
        return local[index];
    }
    T& operator[](int index) {
        return local[index];
    }
    /**
     * Get the size of the array.
     */
    int getSize() {
        return size;
    }
    /**
     * Get the OpenCL Buffer object.
     */
    cl::Buffer& getDeviceBuffer() {
        return *buffer;
    }
    /**
     * Get a pointer to the host buffer.
     */
    T* getHostBuffer() {
        return &local[0];
    }
    /**
     * Get an element of the host buffer.
     */
    const T& get(int index) const {
        return local[index];
    }
    /**
     * Set an element of the host buffer.
     */
    void set(int index, const T& value) {
        local[index] = value;
    }
    /**
     * Copy the values in a vector to the Buffer.
     */
    void upload(std::vector<T>& data) {
        context.getQueue().enqueueWriteBuffer(*buffer, true, 0, size*sizeof(T), &data[0]);
    }
    /**
     * Copy the values in the Buffer to a vector.
     */
    void download(std::vector<T>& data) const {
        context.getQueue().enqueueReadBuffer(*buffer, true, 0, size*sizeof(T), &data[0]);
    }
    /**
     * Copy the values in the host buffer to the OpenCL Buffer.
     */
    void upload() {
        if (local.size() == 0)
            throw OpenMMException(name+": Called upload() on an OpenCLArray with no host buffer");
        upload(local);
    }
    /**
     * Copy the values in the Buffer to the host buffer.
     */
    void download() {
        if (local.size() == 0)
            throw OpenMMException(name+": Called download() on an OpenCLArray with no host buffer");
        download(local);
    }
private:
    OpenCLContext& context;
    cl::Buffer* buffer;
    std::vector<T> local;
    int size;
    bool ownsBuffer;
    std::string name;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLARRAY_H_*/
