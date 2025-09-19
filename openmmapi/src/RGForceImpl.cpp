/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/RGForceImpl.h"
#include "openmm/kernels.h"
#include <set>
#include <sstream>

using namespace OpenMM;
using namespace std;

RGForceImpl::RGForceImpl(const RGForce& owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
}

RGForceImpl::~RGForceImpl() {
}

void RGForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcRGForceKernel::Name(), context);

    // Check for errors in the specification of particles.
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    set<int> particles;
    for (int i : owner.getParticles()) {
        if (i < 0 || i >= numParticles) {
            stringstream msg;
            msg << "RGForce: Illegal particle index for Rg: ";
            msg << i;
            throw OpenMMException(msg.str());
        }
        if (particles.find(i) != particles.end()) {
            stringstream msg;
            msg << "RGForce: Duplicated particle index for Rg: ";
            msg << i;
            throw OpenMMException(msg.str());
        }
        particles.insert(i);
    }
    kernel.getAs<CalcRGForceKernel>().initialize(context.getSystem(), owner);
}

double RGForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<forceGroup)) != 0)
        return kernel.getAs<CalcRGForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> RGForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcRGForceKernel::Name());
    return names;
}
