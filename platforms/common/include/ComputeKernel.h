#ifndef OPENMM_COMPUTEKERNEL_H_
#define OPENMM_COMPUTEKERNEL_H_

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

#include "ArrayInterface.h"
#include <memory>
#include <string>

namespace OpenMM {

/**
 * This abstract class represents a kernel that can be executed on a computing device.
 * Call createKernel() on a ComputeProgramImpl to create an instance of a platform-specific
 * subclass.  Then call addArg() to specify the values to pass for all of the kernel's arguments.
 * Finally, call execute() to execute the kernel.  If you need to modify the values of kernel
 * arguments between invocations, use setArg() to change the value of an argument.
 * 
 * Instead of referring to this class directly, it is best to use ComputeKernel, which is
 * a typedef for a shared_ptr to a ComputeKernelImpl.  This allows you to treat it as having
 * value semantics, and frees you from having to manage memory.  
 */

class ComputeKernelImpl {
public:
    virtual ~ComputeKernelImpl() {
    }
    /**
     * Get the name of this kernel.
     */
    virtual const std::string& getName() const = 0;
    /**
     * Add an argument to pass the kernel when it is invoked.
     * 
     * @param value     the value to pass to the kernel
     */
    template <class T>
    void addArg(T& value) {
        ArrayInterface* array = dynamic_cast<ArrayInterface*>(&value);
        if (array != NULL)
            addArrayArg(*array);
        else
            addPrimitiveArg(&value, sizeof(value));
    }
    /**
     * Set the value of an argument to pass the kernel when it is invoked.
     * 
     * @param index     the index of the argument to set
     * @param value     the value to pass to the kernel
     */
    template <class T>
    void setArg(int index, T& value) {
        ArrayInterface* array = dynamic_cast<ArrayInterface*>(&value);
        if (array != NULL)
            setArrayArg(index, *array);
        else
            setPrimitiveArg(index, &value, sizeof(value));
    }
    /**
     * Execute this kernel.
     *
     * @param threads      the maximum number of threads that should be used.  Depending on the
     *                     computing device, it may choose to use fewer threads than this number.
     * @param blockSize    the number of threads in each thread block.  If this is omitted, a
     *                     default size that is appropriate for the computing device is used.
     */
    void execute(int workUnits, int blockSize=-1);
protected:
    /**
     * Add an argument to pass the kernel when it is invoked, where the value is a
     * subclass of ArrayInterface.
     * 
     * @param value     the value to pass to the kernel
     */
    virtual void addArrayArg(ArrayInterface& value) = 0;
    /**
     * Add an argument to pass the kernel when it is invoked, where the value is a primitive type.
     * 
     * @param value    a pointer to the argument value
     * @param size     the size of the value in bytes
     */
    virtual void addPrimitiveArg(void* value, int size) = 0;
    /**
     * Add an argument to pass the kernel when it is invoked, where the value is a
     * subclass of ArrayInterface.
     * 
     * @param index     the index of the argument to set
     * @param value     the value to pass to the kernel
     */
    virtual void setArrayArg(int index, ArrayInterface& value) = 0;
    /**
     * Add an argument to pass the kernel when it is invoked, where the value is a primitive type.
     * 
     * @param index     the index of the argument to set
     * @param value    a pointer to the argument value
     * @param size     the size of the value in bytes
     */
    virtual void setPrimitiveArg(int index, void* value, int size) = 0;
};

typedef std::shared_ptr<ComputeKernelImpl> ComputeKernel;

} // namespace OpenMM

#endif /*OPENMM_COMPUTEKERNEL_H_*/
