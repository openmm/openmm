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

#include "openmm/internal/CavityForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/kernels.h"
#include <vector>

using namespace OpenMM;
using std::vector;

CavityForceImpl::CavityForceImpl(const CavityForce& owner) : 
        owner(owner), harmonicEnergy(0), couplingEnergy(0), dipoleSelfEnergy(0) {
}

void CavityForceImpl::initialize(ContextImpl& context) {
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    
    // Validate cavity particle index
    int cavityIdx = owner.getCavityParticleIndex();
    if (cavityIdx < 0 || cavityIdx >= numParticles)
        throw OpenMMException("CavityForce: cavity particle index out of range");
    
    // Create the kernel
    kernel = context.getPlatform().createKernel(CalcCavityForceKernel::Name(), context);
    kernel.getAs<CalcCavityForceKernel>().initialize(context.getSystem(), owner);
}

void CavityForceImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
    // Nothing to do here - the cavity force doesn't modify state
}

double CavityForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups & (1 << owner.getForceGroup())) == 0)
        return 0.0;
    
    // Execute the kernel and get energy components
    double totalEnergy = kernel.getAs<CalcCavityForceKernel>().execute(context, includeForces, includeEnergy);
    
    // Retrieve energy components from the kernel
    harmonicEnergy = kernel.getAs<CalcCavityForceKernel>().getHarmonicEnergy();
    couplingEnergy = kernel.getAs<CalcCavityForceKernel>().getCouplingEnergy();
    dipoleSelfEnergy = kernel.getAs<CalcCavityForceKernel>().getDipoleSelfEnergy();
    
    return totalEnergy;
}

std::map<std::string, double> CavityForceImpl::getDefaultParameters() {
    return {
        {CavityForce::OmegaC(), owner.getOmegac()},
        {CavityForce::LambdaCoupling(), owner.getLambdaCoupling()},
        {CavityForce::PhotonMass(), owner.getPhotonMass()}
    };
}

std::vector<std::string> CavityForceImpl::getKernelNames() {
    return {CalcCavityForceKernel::Name()};
}

void CavityForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcCavityForceKernel>().copyParametersToContext(context, owner);
}
