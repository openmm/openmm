/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
#include "openmm/internal/NonbondedSoftcoreForceImpl.h"
#include "openmm/freeEnergyKernels.h"
#include <sstream>

using namespace OpenMM;
using std::pair;
using std::vector;
using std::set;
using std::stringstream;

NonbondedSoftcoreForceImpl::NonbondedSoftcoreForceImpl(NonbondedSoftcoreForce& owner) : owner(owner) {
}

NonbondedSoftcoreForceImpl::~NonbondedSoftcoreForceImpl() {
}

void NonbondedSoftcoreForceImpl::initialize(ContextImpl& context) {

    kernel = context.getPlatform().createKernel(CalcNonbondedSoftcoreForceKernel::Name(), context);

    // Check for errors in the specification of exceptions.

    System& system = context.getSystem();
    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("NonbondedSoftcoreForce must have exactly as many particles as the System it belongs to.");
    vector<set<int> > exceptions(owner.getNumParticles());
    for (int i = 0; i < owner.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon, softcoreLJLambda;
        owner.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda);
        if (particle1 < 0 || particle1 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "NonbondedSoftcoreForce: Illegal particle index for an exception: ";
            msg << particle1;
            throw OpenMMException(msg.str());
        }
        if (particle2 < 0 || particle2 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "NonbondedSoftcoreForce: Illegal particle index for an exception: ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        if (exceptions[particle1].count(particle2) > 0 || exceptions[particle2].count(particle1) > 0) {
            stringstream msg;
            msg << "NonbondedSoftcoreForce: Multiple exceptions are specified for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exceptions[particle1].insert(particle2);
        exceptions[particle2].insert(particle1);
    }
    dynamic_cast<CalcNonbondedSoftcoreForceKernel&>(kernel.getImpl()).initialize(context.getSystem(), owner);
}

double NonbondedSoftcoreForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return dynamic_cast<CalcNonbondedSoftcoreForceKernel&>(kernel.getImpl()).execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> NonbondedSoftcoreForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcNonbondedSoftcoreForceKernel::Name());
    return names;
}
