#ifndef OPENMM_CUSTOMCPPFORCEIMPL_H_
#define OPENMM_CUSTOMCPPFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2023 Stanford University and the Authors.           *
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

#include "ForceImpl.h"
#include "openmm/Kernel.h"
#include "openmm/Vec3.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This is the internal implementation of CustomBondForce.
 */

class CustomCPPForceImpl : public ForceImpl {
public:
    CustomCPPForceImpl(const Force& owner);
    void initialize(ContextImpl& context);
    virtual double computeForce(ContextImpl& context, const std::vector<Vec3>& positions, std::vector<Vec3>& forces) = 0;
    /**
     * Override this if the force updates the context state directly.
     */
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
    }
    /**
     * Override this if the force defines global parameters.
     */
    std::map<std::string, double> getDefaultParameters() {
        return {};
    }
    /**
     * Override this if the force defines bonds between particles.
     */
    std::vector<std::pair<int, int> > getBondedParticles() const {
        return {};
    }
    /**
     * Override this if the force defines an updateParametersInContext() method.
     */
    void updateParametersInContext(ContextImpl& context) {
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::vector<std::string> getKernelNames();
private:
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMCPPFORCEIMPL_H_*/
