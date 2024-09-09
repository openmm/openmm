/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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
#include "openmm/internal/PeriodicTorsionForceImpl.h"
#include "openmm/kernels.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

PeriodicTorsionForceImpl::PeriodicTorsionForceImpl(const PeriodicTorsionForce& owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
}

PeriodicTorsionForceImpl::~PeriodicTorsionForceImpl() {
}

void PeriodicTorsionForceImpl::initialize(ContextImpl& context) {
    const System& system = context.getSystem();
    for (int i = 0; i < owner.getNumTorsions(); i++) {
        int particle[4], periodicity;
        double phase, k;
        owner.getTorsionParameters(i, particle[0], particle[1], particle[2], particle[3], periodicity, phase, k);
        for (int j = 0; j < 4; j++) {
            if (particle[j] < 0 || particle[j] >= system.getNumParticles()) {
                stringstream msg;
                msg << "PeriodicTorsionForce: Illegal particle index for a torsion: ";
                msg << particle[j];
                throw OpenMMException(msg.str());
            }
        }
        if (periodicity < 1)
            throw OpenMMException("PeriodicTorsionForce: periodicity must be positive");
    }
    kernel = context.getPlatform().createKernel(CalcPeriodicTorsionForceKernel::Name(), context);
    kernel.getAs<CalcPeriodicTorsionForceKernel>().initialize(context.getSystem(), owner);
}

double PeriodicTorsionForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<forceGroup)) != 0)
        return kernel.getAs<CalcPeriodicTorsionForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> PeriodicTorsionForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcPeriodicTorsionForceKernel::Name());
    return names;
}

void PeriodicTorsionForceImpl::updateParametersInContext(ContextImpl& context, int firstTorsion, int lastTorsion) {
    kernel.getAs<CalcPeriodicTorsionForceKernel>().copyParametersToContext(context, owner, firstTorsion, lastTorsion);
    context.systemChanged();
}
