/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
 * Authors:                                                                   *
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
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaCTForceImpl.h"
#include "openmm/amoebaKernels.h"
#include <map>
#include <cmath>

using namespace OpenMM;
using namespace std;

using std::pair;
using std::vector;
using std::set;

AmoebaCTForceImpl::AmoebaCTForceImpl(const AmoebaCTForce& owner)
    : owner(owner) {}

AmoebaCTForceImpl::~AmoebaCTForceImpl() {}

void AmoebaCTForceImpl::initialize(ContextImpl& context) {
    const System& system = context.getSystem();

    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("AmoebaCTForce must have exactly as many particles as the System it belongs to.");

    // check that cutoff < 0.5*boxSize

    if (owner.getNonbondedMethod() == AmoebaCTForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5 * boxVectors[0][0] || cutoff > 0.5 * boxVectors[1][1] || cutoff > 0.5 * boxVectors[2][2])
            throw OpenMMException("AmoebaCTForce: The cutoff distance cannot be greater than half the periodic box size.");
    }

    kernel = context.getPlatform().createKernel(CalcAmoebaCTForceKernel::Name(), context);
    kernel.getAs<CalcAmoebaCTForceKernel>().initialize(context.getSystem(), owner);
}

double AmoebaCTForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups & (1 << owner.getForceGroup())) != 0)
        return kernel.getAs<CalcAmoebaCTForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> AmoebaCTForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcAmoebaCTForceKernel::Name());
    return names;
}

void AmoebaCTForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcAmoebaCTForceKernel>().copyParametersToContext(context, owner);
}
