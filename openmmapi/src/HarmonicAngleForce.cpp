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
#include "openmm/HarmonicAngleForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/HarmonicAngleForceImpl.h"

using namespace OpenMM;
using namespace std;

HarmonicAngleForce::HarmonicAngleForce() : usePeriodic(false), numContexts(0) {
}

int HarmonicAngleForce::addAngle(int particle1, int particle2, int particle3, double angle, double k) {
    angles.push_back(AngleInfo(particle1, particle2, particle3, angle, k));
    return angles.size()-1;
}

void HarmonicAngleForce::getAngleParameters(int index, int& particle1, int& particle2, int& particle3, double& angle, double& k) const {
    ASSERT_VALID_INDEX(index, angles);
    particle1 = angles[index].particle1;
    particle2 = angles[index].particle2;
    particle3 = angles[index].particle3;
    angle = angles[index].angle;
    k = angles[index].k;
}

void HarmonicAngleForce::setAngleParameters(int index, int particle1, int particle2, int particle3, double angle, double k) {
    ASSERT_VALID_INDEX(index, angles);
    angles[index].particle1 = particle1;
    angles[index].particle2 = particle2;
    angles[index].particle3 = particle3;
    angles[index].angle = angle;
    angles[index].k = k;
    if (numContexts > 0) {
        firstChangedAngle = min(index, firstChangedAngle);
        lastChangedAngle = max(index, lastChangedAngle);
    }
}

ForceImpl* HarmonicAngleForce::createImpl() const {
    if (numContexts == 0) {
        // Begin tracking changes to angles.
        firstChangedAngle = angles.size();
        lastChangedAngle = -1;
    }
    numContexts++;
    return new HarmonicAngleForceImpl(*this);
}

void HarmonicAngleForce::updateParametersInContext(Context& context) {
    dynamic_cast<HarmonicAngleForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context), firstChangedAngle, lastChangedAngle);
    if (numContexts == 1) {
        // We just updated the only existing context for this force, so we can reset
        // the tracking of changed angles.
        firstChangedAngle = angles.size();
        lastChangedAngle = -1;
    }
}

void HarmonicAngleForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool HarmonicAngleForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
