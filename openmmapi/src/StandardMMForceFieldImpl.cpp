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

StandardMMForceFieldImpl::StandardMMForceFieldImpl(StandardMMForceField& owner, OpenMMContextImpl& context) : owner(owner) {
    kernel = context.getPlatform().createKernel(CalcStandardMMForceFieldKernel::Name());
    vector<vector<int> > bondIndices(owner.getNumBonds());
    vector<vector<double> > bondParameters(owner.getNumBonds());
    vector<vector<int> > angleIndices(owner.getNumAngles());
    vector<vector<double> > angleParameters(owner.getNumAngles());
    vector<vector<int> > periodicTorsionIndices(owner.getNumPeriodicTorsions());
    vector<vector<double> > periodicTorsionParameters(owner.getNumPeriodicTorsions());
    vector<vector<int> > rbTorsionIndices(owner.getNumRBTorsions());
    vector<vector<double> > rbTorsionParameters(owner.getNumRBTorsions());
    vector<vector<int> > bonded14Indices;
    vector<set<int> > exclusions(owner.getNumAtoms());
    vector<vector<double> > nonbondedParameters(owner.getNumAtoms());
    for (int i = 0; i < owner.getNumBonds(); ++i) {
        int atom1, atom2;
        double length, k;
        owner.getBondParameters(i, atom1, atom2, length, k);
        bondIndices[i].push_back(atom1);
        bondIndices[i].push_back(atom2);
        bondParameters[i].push_back(length);
        bondParameters[i].push_back(k);
    }
    for (int i = 0; i < owner.getNumAngles(); ++i) {
        int atom1, atom2, atom3;
        double angle, k;
        owner.getAngleParameters(i, atom1, atom2, atom3, angle, k);
        angleIndices[i].push_back(atom1);
        angleIndices[i].push_back(atom2);
        angleIndices[i].push_back(atom3);
        angleParameters[i].push_back(angle);
        angleParameters[i].push_back(k);
    }
    for (int i = 0; i < owner.getNumPeriodicTorsions(); ++i) {
        int atom1, atom2, atom3, atom4, periodicity;
        double phase, k;
        owner.getPeriodicTorsionParameters(i, atom1, atom2, atom3, atom4, periodicity, phase, k);
        periodicTorsionIndices[i].push_back(atom1);
        periodicTorsionIndices[i].push_back(atom2);
        periodicTorsionIndices[i].push_back(atom3);
        periodicTorsionIndices[i].push_back(atom4);
        periodicTorsionParameters[i].push_back(k);
        periodicTorsionParameters[i].push_back(phase);
        periodicTorsionParameters[i].push_back(periodicity);
    }
    for (int i = 0; i < owner.getNumRBTorsions(); ++i) {
        int atom1, atom2, atom3, atom4;
        double c0, c1, c2, c3, c4, c5;
        owner.getRBTorsionParameters(i, atom1, atom2, atom3, atom4, c0, c1, c2, c3, c4, c5);
        rbTorsionIndices[i].push_back(atom1);
        rbTorsionIndices[i].push_back(atom2);
        rbTorsionIndices[i].push_back(atom3);
        rbTorsionIndices[i].push_back(atom4);
        rbTorsionParameters[i].push_back(c0);
        rbTorsionParameters[i].push_back(c1);
        rbTorsionParameters[i].push_back(c2);
        rbTorsionParameters[i].push_back(c3);
        rbTorsionParameters[i].push_back(c4);
        rbTorsionParameters[i].push_back(c5);
    }
    for (int i = 0; i < owner.getNumAtoms(); ++i) {
        double charge, radius, depth;
        owner.getAtomParameters(i, charge, radius, depth);
        nonbondedParameters[i].push_back(charge);
        nonbondedParameters[i].push_back(radius);
        nonbondedParameters[i].push_back(depth);
    }
    set<pair<int, int> > bonded14set;
    findExclusions(bondIndices, exclusions, bonded14set);
    bonded14Indices.resize(bonded14set.size());
    int index = 0;
    for (set<pair<int, int> >::const_iterator iter = bonded14set.begin(); iter != bonded14set.end(); ++iter) {
        bonded14Indices[index].push_back(iter->first);
        bonded14Indices[index++].push_back(iter->second);
    }
    dynamic_cast<CalcStandardMMForceFieldKernel&>(kernel.getImpl()).initialize(bondIndices, bondParameters, angleIndices, angleParameters,
            periodicTorsionIndices, periodicTorsionParameters, rbTorsionIndices, rbTorsionParameters, bonded14Indices, 0.5, 1.0/1.2, exclusions, nonbondedParameters);
}

StandardMMForceFieldImpl::~StandardMMForceFieldImpl() {
}

void StandardMMForceFieldImpl::calcForces(OpenMMContextImpl& context, Stream& forces) {
    dynamic_cast<CalcStandardMMForceFieldKernel&>(kernel.getImpl()).executeForces(context.getPositions(), forces);
}

double StandardMMForceFieldImpl::calcEnergy(OpenMMContextImpl& context) {
    return dynamic_cast<CalcStandardMMForceFieldKernel&>(kernel.getImpl()).executeEnergy(context.getPositions());
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

