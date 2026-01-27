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
#include "openmm/internal/PythonForceImpl.h"
#include "openmm/kernels.h"
#include <set>
#include <sstream>

using namespace OpenMM;
using namespace std;

PythonForceImpl::PythonForceImpl(const PythonForce& owner) : owner(owner), computation(owner.getComputation()),
        defaultParameters(owner.getGlobalParameters()), usePeriodic(owner.usesPeriodicBoundaryConditions()) {
    forceGroup = owner.getForceGroup();
}

PythonForceImpl::~PythonForceImpl() {
}

void PythonForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcPythonForceKernel::Name(), context);
    kernel.getAs<CalcPythonForceKernel>().initialize(context.getSystem(), owner);
}

double PythonForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<forceGroup)) != 0)
        return kernel.getAs<CalcPythonForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

double PythonForceImpl::calcForcesAndEnergyBatched(ContextImpl& context,
                                                     const std::vector<std::vector<Vec3>>& allBeadPositions,
                                                     std::vector<std::vector<Vec3>>& allBeadForces)
{
    if (allBeadPositions.empty())
        return 0.0;
    
    int numBeads = allBeadPositions.size();
    int numParticles = allBeadPositions[0].size();
    
    // Build vector of State objects for all beads
    // Use vector of pointers to avoid State copy issues
    std::vector<State> stateStorage;
    std::vector<State> states;
    stateStorage.reserve(numBeads);
    states.reserve(numBeads);
    
    Vec3 boxVectors[3];
    bool hasBox = false;
    if (usePeriodic) {
        context.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        hasBox = true;
    }
    
    for (const auto& beadPos : allBeadPositions) {
        State::StateBuilder builder(context.getTime(), context.getStepCount());
        builder.setPositions(beadPos);
        builder.setParameters(context.getParameters());
        
        if (hasBox) {
            builder.setPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        }
        
        stateStorage.push_back(builder.getState());
    }
    
    // Now copy references (should be safe since stateStorage won't be modified)
    for (auto& state : stateStorage) {
        states.push_back(state);
    }
    
    // Allocate flat forces array
    std::vector<double> forcesFlat(3 * numParticles * numBeads);
    double totalEnergy = 0.0;
    
    // Call batched computation
    computation.computeBatch(states, totalEnergy, forcesFlat.data(), true);
    
    // Unpack forces from flat array to vector of vectors
    for (int b = 0; b < numBeads; b++) {
        allBeadForces[b].resize(numParticles);
        for (int p = 0; p < numParticles; p++) {
            int idx = (b * numParticles + p) * 3;
            allBeadForces[b][p] = Vec3(
                forcesFlat[idx],
                forcesFlat[idx+1],
                forcesFlat[idx+2]
            );
        }
    }
    
    return totalEnergy;
}

vector<string> PythonForceImpl::getKernelNames() {
    return {CalcCustomCPPForceKernel::Name()};
}

map<string, double> PythonForceImpl::getDefaultParameters() {
    return defaultParameters;
}
