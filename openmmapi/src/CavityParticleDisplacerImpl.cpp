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

#include "openmm/internal/CavityParticleDisplacerImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/kernels.h"
#include <vector>

using namespace OpenMM;
using std::vector;

CavityParticleDisplacerImpl::CavityParticleDisplacerImpl(const CavityParticleDisplacer& owner) : 
        owner(owner), stepCount(0), hasTriggered(false), lastCoupling(0.0) {
}

void CavityParticleDisplacerImpl::initialize(ContextImpl& context) {
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    
    // Validate cavity particle index
    int cavityIdx = owner.getCavityParticleIndex();
    if (cavityIdx < 0 || cavityIdx >= numParticles)
        throw OpenMMException("CavityParticleDisplacer: cavity particle index out of range");
    
    // Create the kernel
    kernel = context.getPlatform().createKernel(ApplyCavityDisplacementKernel::Name(), context);
    kernel.getAs<ApplyCavityDisplacementKernel>().initialize(context.getSystem(), owner);
}

void CavityParticleDisplacerImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
    // Check if we should trigger the displacement
    int switchStep = owner.getSwitchOnStep();
    double switchLambda = owner.getSwitchOnLambda();
    
    // Trigger displacement when:
    // 1. We've reached the switch-on step
    // 2. The coupling is non-zero
    // 3. We haven't already triggered for this transition
    if (stepCount == switchStep && switchLambda > 0.0 && !hasTriggered) {
        displaceToEquilibrium(context, switchLambda);
        hasTriggered = true;
        forcesInvalid = true;  // Forces need to be recalculated after displacement
    }
    
    stepCount++;
}

void CavityParticleDisplacerImpl::displaceToEquilibrium(ContextImpl& context, double lambdaCoupling) {
    kernel.getAs<ApplyCavityDisplacementKernel>().execute(context, lambdaCoupling);
}

std::map<std::string, double> CavityParticleDisplacerImpl::getDefaultParameters() {
    return {};
}

std::vector<std::string> CavityParticleDisplacerImpl::getKernelNames() {
    return {ApplyCavityDisplacementKernel::Name()};
}
