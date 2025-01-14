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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/PeriodicTorsionForceImpl.h"

using namespace OpenMM;
using namespace std;

PeriodicTorsionForce::PeriodicTorsionForce() : usePeriodic(false), numContexts(0) {
}

int PeriodicTorsionForce::addTorsion(int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    periodicTorsions.push_back(PeriodicTorsionInfo(particle1, particle2, particle3, particle4, periodicity, phase, k));
    return periodicTorsions.size()-1;
}

void PeriodicTorsionForce::getTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, int& periodicity, double& phase, double& k) const {
    ASSERT_VALID_INDEX(index, periodicTorsions);
    particle1 = periodicTorsions[index].particle1;
    particle2 = periodicTorsions[index].particle2;
    particle3 = periodicTorsions[index].particle3;
    particle4 = periodicTorsions[index].particle4;
    periodicity = periodicTorsions[index].periodicity;
    phase = periodicTorsions[index].phase;
    k = periodicTorsions[index].k;
}

void PeriodicTorsionForce::setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    ASSERT_VALID_INDEX(index, periodicTorsions);
    periodicTorsions[index].particle1 = particle1;
    periodicTorsions[index].particle2 = particle2;
    periodicTorsions[index].particle3 = particle3;
    periodicTorsions[index].particle4 = particle4;
    periodicTorsions[index].periodicity = periodicity;
    periodicTorsions[index].phase = phase;
    periodicTorsions[index].k = k;
    if (numContexts > 0) {
        firstChangedTorsion = min(index, firstChangedTorsion);
        lastChangedTorsion = max(index, lastChangedTorsion);
    }
}

ForceImpl* PeriodicTorsionForce::createImpl() const {
    if (numContexts == 0) {
        // Begin tracking changes to torsions.
        firstChangedTorsion = periodicTorsions.size();
        lastChangedTorsion = -1;
    }
    numContexts++;
    return new PeriodicTorsionForceImpl(*this);
}

void PeriodicTorsionForce::updateParametersInContext(Context& context) {
    dynamic_cast<PeriodicTorsionForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context), firstChangedTorsion, lastChangedTorsion);
    if (numContexts == 1) {
        // We just updated the only existing context for this force, so we can reset
        // the tracking of changed torsions.
        firstChangedTorsion = periodicTorsions.size();
        lastChangedTorsion = -1;
    }
}

void PeriodicTorsionForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool PeriodicTorsionForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
