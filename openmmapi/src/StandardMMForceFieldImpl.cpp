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
#include "internal/StandardMMForceFieldImpl.h"
#include "kernels.h"

using namespace OpenMM;
using std::pair;
using std::vector;
using std::set;

StandardMMForceFieldImpl::StandardMMForceFieldImpl(StandardMMForceField& owner) : owner(owner) {
}

StandardMMForceFieldImpl::~StandardMMForceFieldImpl() {
}

void StandardMMForceFieldImpl::initialize(OpenMMContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcStandardMMForceFieldKernel::Name(), context);
    vector<set<int> > exclusions(owner.getNumAtoms());
    vector<vector<int> > bondIndices(owner.getNumBonds());
    set<pair<int, int> > bonded14set;
    for (int i = 0; i < owner.getNumBonds(); ++i) {
        int atom1, atom2;
        double length, k;
        owner.getBondParameters(i, atom1, atom2, length, k);
        bondIndices[i].push_back(atom1);
        bondIndices[i].push_back(atom2);
    }
    findExclusions(bondIndices, exclusions, bonded14set);
    dynamic_cast<CalcStandardMMForceFieldKernel&>(kernel.getImpl()).initialize(context.getSystem(), owner, exclusions);
}

void StandardMMForceFieldImpl::calcForces(OpenMMContextImpl& context, Stream& forces) {
    dynamic_cast<CalcStandardMMForceFieldKernel&>(kernel.getImpl()).executeForces(context);
}

double StandardMMForceFieldImpl::calcEnergy(OpenMMContextImpl& context) {
    return dynamic_cast<CalcStandardMMForceFieldKernel&>(kernel.getImpl()).executeEnergy(context);
}

std::vector<std::string> StandardMMForceFieldImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcStandardMMForceFieldKernel::Name());
    return names;
}

void StandardMMForceFieldImpl::findExclusions(const vector<vector<int> >& bondIndices, vector<set<int> >& exclusions, set<pair<int, int> >& bonded14Indices) const {
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

void StandardMMForceFieldImpl::addExclusionsToSet(const vector<set<int> >& bonded12, set<int>& exclusions, int baseAtom, int fromAtom, int currentLevel) const {
    for (set<int>::const_iterator iter = bonded12[fromAtom].begin(); iter != bonded12[fromAtom].end(); ++iter) {
        if (*iter != baseAtom)
            exclusions.insert(*iter);
        if (currentLevel > 0)
            addExclusionsToSet(bonded12, exclusions, baseAtom, *iter, currentLevel-1);
    }
}

