#ifndef OPENMM_ARRAYINTERFACE_H_
#define OPENMM_ARRAYINTERFACE_H_

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

#include "openmm/OpenMMException.h"
#include "openmm/common/windowsExportCommon.h"
#include <vector>

namespace OpenMM {

class ComputeContext;

/**
 * This abstract class defines the interface for arrays stored on a computing device.
 */

class OPENMM_EXPORT_COMMON ArrayInterface {
public:
    virtual ~ArrayInterface() {
    }
    /**
     * Initialize this array.
     *
     * @param context           the context for which to create the array
     * @param size              the number of elements in the array
     * @param elementSize       the size of each element in bytes
     * @param name              the name of the array
     */
    virtual void initialize(ComputeContext& context, int size, int elementSize, const std::string& name) = 0;
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
    virtual void resize(int size) = 0;
    /**
     * Get whether this array has been initialized.
     */
    virtual bool isInitialized() const = 0;
    /**
     * Get the number of elements in the array.
     */
    virtual int getSize() const = 0;
    /**
     * Get the size of each element in bytes.
     */
    virtual int getElementSize() const = 0;
    /**
     * Get the name of the array.
     */
    virtual const std::string& getName() const = 0;
    /**
     * Get the context this array belongs to.
     */
    virtual ComputeContext& getContext() = 0;
    /**
     * Copy the values in a vector to the device memory.
     * 
     * @param data      the data in host memory to copy
     * @param convert   if true, automatic conversions between single and double
     *                  precision will be performed as necessary
     */
    template <class T>
    void upload(const std::vector<T>& data, bool convert=false) {
        if (convert && data.size() == getSize() && sizeof(T) != getElementSize()) {
            if (sizeof(T) == 2*getElementSize()) {
                // Convert values from double to single precision.
                const double* d = reinterpret_cast<const double*>(&data[0]);
                std::vector<float> v(getElementSize()*getSize()/sizeof(float));
                for (int i = 0; i < v.size(); i++)
                    v[i] = (float) d[i];
                upload(&v[0], true);
                return;
            }
            if (2*sizeof(T) == getElementSize()) {
                // Convert values from single to double precision.
                const float* d = reinterpret_cast<const float*>(&data[0]);
                std::vector<double> v(getElementSize()*getSize()/sizeof(double));
                for (int i = 0; i < v.size(); i++)
                    v[i] = (double) d[i];
                upload(&v[0], true);
                return;
            }
        }
        if (sizeof(T) != getElementSize() || data.size() != getSize())
            throw OpenMMException("Error uploading array "+getName()+": The specified vector does not match the size of the array");
        upload(&data[0], true);
    }
    /**
     * Copy the values in the array to a vector.
     */
    template <class T>
    void download(std::vector<T>& data) const {
        if (sizeof(T) != getElementSize())
            throw OpenMMException("Error downloading array "+getName()+": The specified vector has the wrong element size");
        if (data.size() != getSize())
            data.resize(getSize());
        download(&data[0], true);
    }
    /**
     * Copy the values from host memory to the array.
     * 
     * @param data     the data to copy
     * @param blocking if true, this call will block until the transfer is complete.  Subclasses often
     *                 have restrictions on non-blocking copies, such as that the source data must be
     *                 in page-locked memory.
     */
    virtual void upload(const void* data, bool blocking=true) = 0;
    /**
     * Copy the values in the array to host memory.
     * 
     * @param data     the destination to copy the value to
     * @param blocking if true, this call will block until the transfer is complete.  Subclasses often
     *                 have restrictions on non-blocking copies, such as that the destination must be
     *                 in page-locked memory.
     */
    virtual void download(void* data, bool blocking=true) const = 0;
    /**
     * Copy the values in this array to a second array.
     * 
     * @param dest     the destination array to copy to
     */
    virtual void copyTo(ArrayInterface& dest) const = 0;
};

} // namespace OpenMM

#endif /*OPENMM_ARRAYINTERFACE_H_*/
