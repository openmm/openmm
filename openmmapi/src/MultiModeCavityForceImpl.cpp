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

#include "openmm/internal/MultiModeCavityForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/kernels.h"
#include <sstream>
#include <vector>

using namespace OpenMM;
using std::vector;

MultiModeCavityForceImpl::MultiModeCavityForceImpl(const MultiModeCavityForce& owner) :
        owner(owner), harmonicEnergy(0), couplingEnergy(0), dipoleSelfEnergy(0) {
}

void MultiModeCavityForceImpl::initialize(ContextImpl& context) {
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    
    // Validate that the correct number of cavity particles have been added
    const auto& indices = owner.getCavityParticleIndices();
    if ((int)indices.size() != owner.getNumModes()) {
        std::stringstream msg;
        msg << "MultiModeCavityForce: expected " << owner.getNumModes()
            << " cavity particles but only " << indices.size() << " were added";
        throw OpenMMException(msg.str());
    }
    
    // Validate all cavity particle indices
    for (int i = 0; i < (int)indices.size(); i++) {
        if (indices[i] < 0 || indices[i] >= numParticles) {
            std::stringstream msg;
            msg << "MultiModeCavityForce: cavity particle index " << indices[i]
                << " for mode " << (i+1) << " is out of range [0, " << numParticles << ")";
            throw OpenMMException(msg.str());
        }
    }
    
    // Create the kernel
    kernel = context.getPlatform().createKernel(CalcMultiModeCavityForceKernel::Name(), context);
    kernel.getAs<CalcMultiModeCavityForceKernel>().initialize(context.getSystem(), owner);
}

void MultiModeCavityForceImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
    // Nothing to do here
}

double MultiModeCavityForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups & (1 << owner.getForceGroup())) == 0)
        return 0.0;
    
    // Execute the kernel and get energy components
    double totalEnergy = kernel.getAs<CalcMultiModeCavityForceKernel>().execute(context, includeForces, includeEnergy);
    
    // Retrieve energy components from the kernel
    harmonicEnergy = kernel.getAs<CalcMultiModeCavityForceKernel>().getHarmonicEnergy();
    couplingEnergy = kernel.getAs<CalcMultiModeCavityForceKernel>().getCouplingEnergy();
    dipoleSelfEnergy = kernel.getAs<CalcMultiModeCavityForceKernel>().getDipoleSelfEnergy();
    
    return totalEnergy;
}

std::map<std::string, double> MultiModeCavityForceImpl::getDefaultParameters() {
    return {};
}

std::vector<std::string> MultiModeCavityForceImpl::getKernelNames() {
    return {CalcMultiModeCavityForceKernel::Name()};
}

void MultiModeCavityForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcMultiModeCavityForceKernel>().copyParametersToContext(context, owner);
}
