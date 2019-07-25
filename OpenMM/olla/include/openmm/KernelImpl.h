#ifndef OPENMM_KERNELIMPL_H_
#define OPENMM_KERNELIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "Platform.h"
#include <string>
#include <cassert>
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

/**
 * A KernelImpl defines the internal implementation of a Kernel object.  A subclass will typically
 * declare an abstract execute() method which defines the API for executing the kernel.  Other classes
 * will in turn subclass it and provide concrete implementations of the execute() method.
 */

class OPENMM_EXPORT KernelImpl {
public:
    /**
     * Create a KernelImpl.
     * 
     * @param name      the name of the kernel to create
     * @param platform  the Platform that created this kernel
     */
    KernelImpl(std::string name, const Platform& platform);
    virtual ~KernelImpl() {
        assert(referenceCount == 0);
    }
    /**
     * Get the name of this kernel.
     */
    std::string getName() const;
    /**
     * Get the Platform that created this KernelImpl.
     */
    const Platform& getPlatform();
private:
    friend class Kernel;
    std::string name;
    const Platform* platform;
    int referenceCount;
};

} // namespace OpenMM

#endif /*OPENMM_KERNELIMPL_H_*/
