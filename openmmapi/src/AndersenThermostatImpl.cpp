/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2010 Stanford University and the Authors.      *
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

#include "openmm/internal/AndersenThermostatImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Integrator.h"
#include "openmm/System.h"
#include "openmm/kernels.h"
#include <vector>

using namespace OpenMM;
using std::vector;

AndersenThermostatImpl::AndersenThermostatImpl(const AndersenThermostat& owner) : owner(owner) {
}

void AndersenThermostatImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(ApplyAndersenThermostatKernel::Name(), context);
    kernel.getAs<ApplyAndersenThermostatKernel>().initialize(context.getSystem(), owner);
}

void AndersenThermostatImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
    kernel.getAs<ApplyAndersenThermostatKernel>().execute(context);
}

std::map<std::string, double> AndersenThermostatImpl::getDefaultParameters() {
    std::map<std::string, double> parameters;
    parameters[AndersenThermostat::Temperature()] = getOwner().getDefaultTemperature();
    parameters[AndersenThermostat::CollisionFrequency()] = getOwner().getDefaultCollisionFrequency();
    return parameters;
}

std::vector<std::string> AndersenThermostatImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(ApplyAndersenThermostatKernel::Name());
    return names;
}

vector<vector<int> > AndersenThermostatImpl::calcParticleGroups(const System& system) {
    // First make a list of every other particle to which each particle is connected by a constraint.

    int numParticles = system.getNumParticles();
    vector<vector<int> > particleConstraints(numParticles);
    for (int i = 0; i < system.getNumConstraints(); i++) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        particleConstraints[particle1].push_back(particle2);
        particleConstraints[particle2].push_back(particle1);
    }

    // Now tag particles by which molecule they belong to.

    vector<int> particleGroup(numParticles, -1);
    int numGroups = 0;
    for (int i = 0; i < numParticles; i++)
        if (particleGroup[i] == -1)
            tagParticlesInGroup(i, numGroups++, particleGroup, particleConstraints);
    vector<vector<int> > particleIndices(numGroups);
    for (int i = 0; i < numParticles; i++)
        particleIndices[particleGroup[i]].push_back(i);
    return particleIndices;
}

void AndersenThermostatImpl::tagParticlesInGroup(int particle, int group, vector<int>& particleGroup, vector<vector<int> >& particleConstraints) {
    // Recursively tag particles as belonging to a particular group.

    particleGroup[particle] = group;
    for (int constrained : particleConstraints[particle])
        if (particleGroup[constrained] == -1)
            tagParticlesInGroup(constrained, group, particleGroup, particleConstraints);
}
