/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Pretti                                        *
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

#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/LCPOForceImpl.h"
#include "openmm/kernels.h"

using namespace OpenMM;
using namespace std;

LCPOForceImpl::LCPOForceImpl(const LCPOForce& owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
}

void LCPOForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcLCPOForceKernel::Name(), context);

    const System& system = context.getSystem();
    if (owner.getNumParticles() != system.getNumParticles()) {
        throw OpenMMException("LCPOForce must have exactly as many particles as the System it belongs to.");
    }

    for (int i = 0; i < owner.getNumParticles(); i++) {
        double radius, p1, p2, p3, p4;
        owner.getParticleParameters(i, radius, p1, p2, p3, p4);
        if (radius < 0.0) {
            throw OpenMMException("LCPOForce: radius for a particle cannot be negative");
        }
    }

    kernel.getAs<CalcLCPOForceKernel>().initialize(context.getSystem(), owner);
}

double LCPOForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups & (1 << forceGroup)) != 0) {
        return kernel.getAs<CalcLCPOForceKernel>().execute(context, includeForces, includeEnergy);
    }
    return 0.0;
}

std::vector<std::string> LCPOForceImpl::getKernelNames() {
    return {CalcLCPOForceKernel::Name()};
}

void LCPOForceImpl::updateParametersInContext(ContextImpl& context, int firstParticle, int lastParticle) {
    kernel.getAs<CalcLCPOForceKernel>().copyParametersToContext(context, owner, firstParticle, lastParticle);
    context.systemChanged();
}
