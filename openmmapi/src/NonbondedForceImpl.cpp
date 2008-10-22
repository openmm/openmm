/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "internal/OpenMMContextImpl.h"
#include "internal/NonbondedForceImpl.h"
#include "kernels.h"

using namespace OpenMM;
using std::pair;
using std::vector;
using std::set;

NonbondedForceImpl::NonbondedForceImpl(NonbondedForce& owner) : owner(owner) {
}

NonbondedForceImpl::~NonbondedForceImpl() {
}

void NonbondedForceImpl::initialize(OpenMMContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcNonbondedForceKernel::Name(), context);
    
    // See if the system contains a HarmonicBondForce.  If so, use it to identify exclusions.
    
    System& system = context.getSystem();
    vector<set<int> > exclusions(owner.getNumParticles());
    for (int i = 0; i < system.getNumForces(); i++) {
        if (dynamic_cast<HarmonicBondForce*>(&system.getForce(i)) != NULL) {
            const HarmonicBondForce& force = dynamic_cast<const HarmonicBondForce&>(system.getForce(i));
            vector<vector<int> > bondIndices(force.getNumBonds());
            set<pair<int, int> > bonded14set;
            for (int i = 0; i < force.getNumBonds(); ++i) {
                int particle1, particle2;
                double length, k;
                force.getBondParameters(i, particle1, particle2, length, k);
                bondIndices[i].push_back(particle1);
                bondIndices[i].push_back(particle2);
            }
            findExclusions(bondIndices, exclusions, bonded14set);
        }
    }
    dynamic_cast<CalcNonbondedForceKernel&>(kernel.getImpl()).initialize(context.getSystem(), owner, exclusions);
}

void NonbondedForceImpl::calcForces(OpenMMContextImpl& context, Stream& forces) {
    dynamic_cast<CalcNonbondedForceKernel&>(kernel.getImpl()).executeForces(context);
}

double NonbondedForceImpl::calcEnergy(OpenMMContextImpl& context) {
    return dynamic_cast<CalcNonbondedForceKernel&>(kernel.getImpl()).executeEnergy(context);
}

std::vector<std::string> NonbondedForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcNonbondedForceKernel::Name());
    return names;
}

void NonbondedForceImpl::findExclusions(const vector<vector<int> >& bondIndices, vector<set<int> >& exclusions, set<pair<int, int> >& bonded14Indices) const {
    vector<set<int> > bonded12(exclusions.size());
    for (int i = 0; i < (int) bondIndices.size(); ++i) {
        bonded12[bondIndices[i][0]].insert(bondIndices[i][1]);
        bonded12[bondIndices[i][1]].insert(bondIndices[i][0]);
    }
    for (int i = 0; i < (int) exclusions.size(); ++i)
        addExclusionsToSet(bonded12, exclusions[i], i, i, 2);
    for (int i = 0; i < (int) exclusions.size(); ++i) {
        set<int> bonded13;
        addExclusionsToSet(bonded12, bonded13, i, i, 1);
        for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter)
            if (*iter < i && bonded13.find(*iter) == bonded13.end())
                bonded14Indices.insert(pair<int, int> (*iter, i));
    }
}

void NonbondedForceImpl::addExclusionsToSet(const vector<set<int> >& bonded12, set<int>& exclusions, int baseParticle, int fromParticle, int currentLevel) const {
    for (set<int>::const_iterator iter = bonded12[fromParticle].begin(); iter != bonded12[fromParticle].end(); ++iter) {
        if (*iter != baseParticle)
            exclusions.insert(*iter);
        if (currentLevel > 0)
            addExclusionsToSet(bonded12, exclusions, baseParticle, *iter, currentLevel-1);
    }
}

