#ifndef OPENMM_CAVITYPARTICLEDISPLACERIMPL_H_
#define OPENMM_CAVITYPARTICLEDISPLACERIMPL_H_

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
#include "openmm/CavityParticleDisplacer.h"
#include "openmm/Kernel.h"
#include <string>
#include <vector>

namespace OpenMM {

class System;

/**
 * This is the internal implementation of CavityParticleDisplacer.
 */
class OPENMM_EXPORT CavityParticleDisplacerImpl : public ForceImpl {
public:
    CavityParticleDisplacerImpl(const CavityParticleDisplacer& owner);
    void initialize(ContextImpl& context);
    const CavityParticleDisplacer& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid);
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
        // This force doesn't apply forces - it only displaces positions
        return 0.0;
    }
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
    /**
     * Displace the cavity particle to its equilibrium position.
     * 
     * @param context        the context to modify
     * @param lambdaCoupling the coupling value to use
     */
    void displaceToEquilibrium(ContextImpl& context, double lambdaCoupling);
private:
    const CavityParticleDisplacer& owner;
    Kernel kernel;
    int stepCount;
    bool hasTriggered;
    double lastCoupling;
};

} // namespace OpenMM

#endif /*OPENMM_CAVITYPARTICLEDISPLACERIMPL_H_*/
