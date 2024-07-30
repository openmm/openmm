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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/HarmonicAngleForceImpl.h"
#include "openmm/kernels.h"
#include <cmath>
#include <sstream>

using namespace OpenMM;
using namespace std;

HarmonicAngleForceImpl::HarmonicAngleForceImpl(const HarmonicAngleForce& owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
}

HarmonicAngleForceImpl::~HarmonicAngleForceImpl() {
}

void HarmonicAngleForceImpl::initialize(ContextImpl& context) {
    const System& system = context.getSystem();
    for (int i = 0; i < owner.getNumAngles(); i++) {
        int particle[3];
        double angle, k;
        owner.getAngleParameters(i, particle[0], particle[1], particle[2], angle, k);
        for (int j = 0; j < 3; j++) {
            if (particle[j] < 0 || particle[j] >= system.getNumParticles()) {
                stringstream msg;
                msg << "HarmonicAngleForce: Illegal particle index for an angle: ";
                msg << particle[j];
                throw OpenMMException(msg.str());
            }
        }
        if (angle < 0 || angle > M_PI*1.000001)
            throw OpenMMException("HarmonicAngleForce: angle must be between 0 and pi");
    }
    kernel = context.getPlatform().createKernel(CalcHarmonicAngleForceKernel::Name(), context);
    kernel.getAs<CalcHarmonicAngleForceKernel>().initialize(context.getSystem(), owner);
}

double HarmonicAngleForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<forceGroup)) != 0)
        return kernel.getAs<CalcHarmonicAngleForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> HarmonicAngleForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcHarmonicAngleForceKernel::Name());
    return names;
}

void HarmonicAngleForceImpl::updateParametersInContext(ContextImpl& context, int firstAngle, int lastAngle) {
    kernel.getAs<CalcHarmonicAngleForceKernel>().copyParametersToContext(context, owner, firstAngle, lastAngle);
    context.systemChanged();
}
