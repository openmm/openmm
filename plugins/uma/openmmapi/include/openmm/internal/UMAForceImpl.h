#ifndef OPENMM_UMAFORCEIMPL_H_
#define OPENMM_UMAFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                             *
 * Portions copyright (c) 2025 Stanford University and the Authors.            *
 * Authors: Muhammad Hasyim                                                    *
 * Contributors:                                                               *
 *                                                                             *
 * Permission is hereby granted, free of charge, to any person obtaining a     *
 * copy of this software and associated documentation files (the "Software"),  *
 * to deal in the Software without restriction, including without limitation   *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,    *
 * and/or sell copies of the Software, and to permit persons to whom the       *
 * Software is furnished to do so, subject to the following conditions:        *
 *                                                                             *
 * The above copyright notice and this permission notice shall be included in  *
 * all copies or substantial portions of the Software.                         *
 *                                                                             *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL     *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,     *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR       *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE   *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                      *
 * -------------------------------------------------------------------------- */

#include "openmm/UMAForce.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

class System;

/**
 * This is the internal implementation of UMAForce.
 */
class OPENMM_EXPORT_UMA UMAForceImpl : public ForceImpl {
public:
    UMAForceImpl(const UMAForce& owner);
    ~UMAForceImpl();
    
    void initialize(ContextImpl& context);
    
    const UMAForce& getOwner() const {
        return owner;
    }
    
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
        // This force doesn't update the state directly
    }
    
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    
    std::map<std::string, double> getDefaultParameters() {
        return std::map<std::string, double>(); // No adjustable parameters
    }
    
    std::vector<std::string> getKernelNames();
    
    void updateParametersInContext(ContextImpl& context);
    
    /**
     * Get the list of atoms to which this force is applied, or an empty set if it is applied to all atoms.
     */
    std::set<int> getAtomIndices() const;
    
private:
    const UMAForce& owner;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_UMAFORCEIMPL_H_*/
