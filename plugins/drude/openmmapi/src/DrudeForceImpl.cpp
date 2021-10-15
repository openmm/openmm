/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2021 Stanford University and the Authors.      *
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
#include "openmm/internal/DrudeForceImpl.h"
#include "openmm/DrudeKernels.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>

using namespace OpenMM;
using namespace std;

DrudeForceImpl::DrudeForceImpl(const DrudeForce& owner) : owner(owner) {
}

DrudeForceImpl::~DrudeForceImpl() {
}

void DrudeForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcDrudeForceKernel::Name(), context);
    const System& system = context.getSystem();

    // Check for errors in the specification of particles.

    set<int> usedParticles;
    for (int i = 0; i < owner.getNumParticles(); i++) {
        int particle[5];
        double charge, k, k2, k3;
        owner.getParticleParameters(i, particle[0], particle[1], particle[2], particle[3], particle[4], charge, k, k2, k3);
        for (int i = 0; i < 2; i++) {
            if (particle[i] < 0 || particle[i] >= system.getNumParticles()) {
                stringstream msg;
                msg << "DrudeForce: Illegal particle index: ";
                msg << particle[i];
                throw OpenMMException(msg.str());
            }
            if (usedParticles.find(particle[i]) != usedParticles.end()) {
                stringstream msg;
                msg << "DrudeForce: Particle index is used by two different Drude particles: ";
                msg << particle[i];
                throw OpenMMException(msg.str());
            }
            usedParticles.insert(particle[i]);
        }
        for (int i = 2; i < 5; i++) {
            if (particle[i] < -1 || particle[i] >= system.getNumParticles()) {
                stringstream msg;
                msg << "DrudeForce: Illegal particle index: ";
                msg << particle[i];
                throw OpenMMException(msg.str());
            }
        }
    }

    // Check for errors in the specification of screened pairs.

    vector<set<int> > screenedPairs(owner.getNumParticles());
    for (int i = 0; i < owner.getNumScreenedPairs(); i++) {
        int particle[2];
        double thole;
        owner.getScreenedPairParameters(i, particle[0], particle[1], thole);
        for (int i = 0; i < 2; i++) {
            if (particle[i] < 0 || particle[i] >= owner.getNumParticles()) {
                stringstream msg;
                msg << "DrudeForce: Illegal particle index for a screened pair: ";
                msg << particle[i];
                throw OpenMMException(msg.str());
            }
        }
        if (screenedPairs[particle[0]].count(particle[1]) > 0 || screenedPairs[particle[1]].count(particle[0]) > 0) {
            stringstream msg;
            msg << "DrudeForce: Multiple screened pairs are specified for particles ";
            msg << particle[0];
            msg << " and ";
            msg << particle[1];
            throw OpenMMException(msg.str());
        }
        screenedPairs[particle[0]].insert(particle[1]);
        screenedPairs[particle[1]].insert(particle[0]);
    }
    kernel.getAs<CalcDrudeForceKernel>().initialize(context.getSystem(), owner);
}

double DrudeForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcDrudeForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> DrudeForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcDrudeForceKernel::Name());
    return names;
}

void DrudeForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcDrudeForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}

vector<pair<int, int> > DrudeForceImpl::getBondedParticles() const {
    int numParticles = owner.getNumParticles();
    vector<pair<int, int> > bonds(numParticles);
    for (int i = 0; i < numParticles; i++) {
        int p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        owner.getParticleParameters(i, bonds[i].first, bonds[i].second, p2, p3, p4, charge, polarizability, aniso12, aniso34);
    }
    return bonds;
}
