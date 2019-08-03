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
 * Portions copyright (c) 2009-2018 Stanford University and the Authors.      *
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

#define __CL_ENABLE_EXCEPTIONS
#define CL_USE_DEPRECATED_OPENCL_1_1_APIS
#include "openmm/OpenMMException.h"
#include "windowsExportOpenCL.h"
#include <cl.hpp>
#include <iostream>
#include <sstream>
#include <vector>

namespace OpenMM {

class OpenCLContext;

/**
 * This class encapsulates an OpenCL Buffer.  It provides a simplified API for working with it,
 * and for copying data to and from the OpenCL Buffer.
 */

class OPENMM_EXPORT_OPENCL OpenCLArray {
public:
    /**
     * Create an OpenCLArray object.  The object is allocated on the heap with the "new" operator.
     * The template argument is the data type of each array element.
     *
     * @param context           the context for which to create the array
     * @param size              the number of elements in the array
     * @param name              the name of the array
     * @param flags             the set of flags to specify when creating the OpenCL Buffer
     */
    template <class T>
    static OpenCLArray* create(OpenCLContext& context, int size, const std::string& name, cl_int flags = CL_MEM_READ_WRITE) {
        return new OpenCLArray(context, size, sizeof(T), name, flags);
    }
    /**
     * Create an OpenCLArray object that uses a preexisting Buffer.  The object is allocated on the heap with the "new" operator.
     * The template argument is the data type of each array element.
     *
     * @param context           the context for which to create the array
     * @param buffer            the OpenCL Buffer this object encapsulates
     * @param size              the number of elements in the array
     * @param name              the name of the array
     */
    template <class T>
    static OpenCLArray* create(OpenCLContext& context, cl::Buffer* buffer, int size, const std::string& name) {
        return new OpenCLArray(context, buffer, size, sizeof(T), name);
    }
    /**
     * Create an uninitialized OpenCLArray object.  It does not point to any OpenCL Buffer,
     * and cannot be used until initialize() is called on it.
     */
    OpenCLArray();
    /**
     * Create an OpenCLArray object.
     *
     * @param context           the context for which to create the array
     * @param size              the number of elements in the array
     * @param elementSize       the size of each element in bytes
     * @param name              the name of the array
     * @param flags             the set of flags to specify when creating the OpenCL Buffer
     */
    OpenCLArray(OpenCLContext& context, int size, int elementSize, const std::string& name, cl_int flags = CL_MEM_READ_WRITE);
    /**
     * Create an OpenCLArray object that uses a preexisting Buffer.
     *
     * @param context           the context for which to create the array
     * @param buffer            the OpenCL Buffer this object encapsulates
     * @param size              the number of elements in the array
     * @param elementSize       the size of each element in bytes
     * @param name              the name of the array
     */
    OpenCLArray(OpenCLContext& context, cl::Buffer* buffer, int size, int elementSize, const std::string& name);
    ~OpenCLArray();
    /**
     * Initialize this object.
     *
     * @param context           the context for which to create the array
     * @param size              the number of elements in the array
     * @param elementSize       the size of each element in bytes
     * @param name              the name of the array
     * @param flags             the set of flags to specify when creating the OpenCL Buffer
     */
    void initialize(OpenCLContext& context, int size, int elementSize, const std::string& name, cl_int flags = CL_MEM_READ_WRITE);
    /**
     * Initialize this object to use a preexisting Buffer.
     *
     * @param context           the context for which to create the array
     * @param buffer            the OpenCL Buffer this object encapsulates
     * @param size              the number of elements in the array
     * @param elementSize       the size of each element in bytes
     * @param name              the name of the array
     */
    void initialize(OpenCLContext& context, cl::Buffer* buffer, int size, int elementSize, const std::string& name);
    /**
     * Initialize this object.  The template argument is the data type of each array element.
     *
     * @param context           the context for which to create the array
     * @param size              the number of elements in the array
     * @param name              the name of the array
     * @param flags             the set of flags to specify when creating the OpenCL Buffer
     */
    template <class T>
    void initialize(OpenCLContext& context, int size, const std::string& name, cl_int flags = CL_MEM_READ_WRITE) {
        initialize(context, size, sizeof(T), name, flags);
    }
    /**
     * Initialize this object to use a preexisting Buffer.  The template argument
     * is the data type of each array element.
     *
     * @param context           the context for which to create the array
     * @param buffer            the OpenCL Buffer this object encapsulates
     * @param size              the number of elements in the array
     * @param name              the name of the array
     */
    template <class T>
    void initialize(OpenCLContext& context, cl::Buffer* buffer, int size, const std::string& name) {
        initialize(context, buffer, size, sizeof(T), name);
    }
    /**
     * Recreate the internal storage to have a different size.
     */
    void resize(int size);
    /**
     * Get whether this array has been initialized.
     */
    bool isInitialized() const {
        return (buffer != NULL);
    }
    /**
     * Get the size of the array.
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
     * Get the OpenCL Buffer object.
     */
    cl::Buffer& getDeviceBuffer() {
        return *buffer;
    }
    /**
     * Copy the values in a vector to the Buffer.
     */
    template <class T>
    void upload(const std::vector<T>& data, bool blocking = true, bool convert = false) {
        if (convert && data.size() == size && sizeof(T) != elementSize) {
            if (sizeof(T) == 2*elementSize) {
                // Convert values from double to single precision.
                const double* d = reinterpret_cast<const double*>(&data[0]);
                std::vector<float> v(elementSize*size/sizeof(float));
                for (int i = 0; i < v.size(); i++)
                    v[i] = (float) d[i];
                upload(&v[0], blocking);
                return;
            }
            if (2*sizeof(T) == elementSize) {
                // Convert values from single to double precision.
                const float* d = reinterpret_cast<const float*>(&data[0]);
                std::vector<double> v(elementSize*size/sizeof(double));
                for (int i = 0; i < v.size(); i++)
                    v[i] = (double) d[i];
                upload(&v[0], blocking);
                return;
            }
        }
        if (sizeof(T) != elementSize || data.size() != size)
            throw OpenMMException("Error uploading array "+name+": The specified vector does not match the size of the array");
        upload(&data[0], blocking);
    }
    /**
     * Copy the values in the Buffer to a vector.
     */
    template <class T>
    void download(std::vector<T>& data, bool blocking = true) const {
        if (sizeof(T) != elementSize)
            throw OpenMMException("Error downloading array "+name+": The specified vector has the wrong element size");
        if (data.size() != size)
            data.resize(size);
        download(&data[0], blocking);
    }
    /**
     * Copy the values in an array to the Buffer.
     * 
     * @param data     the data to copy
     * @param blocking if true, this call will block until the transfer is complete.
     */
    void upload(const void* data, bool blocking = true);
    /**
     * Copy the values in the Buffer to an array.
     * 
     * @param data     the array to copy the memory to
     * @param blocking if true, this call will block until the transfer is complete.
     */
    void download(void* data, bool blocking = true) const;
    /**
     * Copy the values in the Buffer to a second OpenCLArray.
     * 
     * @param dest     the destination array to copy to
     */
    void copyTo(OpenCLArray& dest) const;
private:
    OpenCLContext* context;
    cl::Buffer* buffer;
    int size, elementSize;
    cl_int flags;
    bool ownsBuffer;
    std::string name;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLARRAY_H_*/
