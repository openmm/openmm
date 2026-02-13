#ifndef OPENMM_MULTIMODECAVITYFORCEIMPL_H_
#define OPENMM_MULTIMODECAVITYFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Muhammad Hasyim                                                   *
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
#include "openmm/MultiModeCavityForce.h"
#include "openmm/Kernel.h"
#include <string>
#include <vector>

namespace OpenMM {

class System;

/**
 * This is the internal implementation of MultiModeCavityForce.
 */
class OPENMM_EXPORT MultiModeCavityForceImpl : public ForceImpl {
public:
    MultiModeCavityForceImpl(const MultiModeCavityForce& owner);
    void initialize(ContextImpl& context);
    const MultiModeCavityForce& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid);
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(ContextImpl& context);
    /**
     * Get the total harmonic energy component from the last force calculation.
     */
    double getHarmonicEnergy() const {
        return harmonicEnergy;
    }
    /**
     * Get the total coupling energy component from the last force calculation.
     */
    double getCouplingEnergy() const {
        return couplingEnergy;
    }
    /**
     * Get the dipole self-energy component from the last force calculation.
     */
    double getDipoleSelfEnergy() const {
        return dipoleSelfEnergy;
    }
private:
    const MultiModeCavityForce& owner;
    Kernel kernel;
    double harmonicEnergy;
    double couplingEnergy;
    double dipoleSelfEnergy;
};

} // namespace OpenMM

#endif /*OPENMM_MULTIMODECAVITYFORCEIMPL_H_*/
