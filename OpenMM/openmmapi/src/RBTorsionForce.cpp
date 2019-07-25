/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
#include "openmm/RBTorsionForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/RBTorsionForceImpl.h"

using namespace OpenMM;

RBTorsionForce::RBTorsionForce() : usePeriodic(false) {
}

int RBTorsionForce::addTorsion(int particle1, int particle2, int particle3, int particle4, double c0, double c1, double c2, double c3, double c4, double c5) {
    rbTorsions.push_back(RBTorsionInfo(particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5));
    return rbTorsions.size()-1;
}

void RBTorsionForce::getTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, double& c0, double& c1, double& c2, double& c3, double& c4, double& c5) const {
    ASSERT_VALID_INDEX(index, rbTorsions);
    particle1 = rbTorsions[index].particle1;
    particle2 = rbTorsions[index].particle2;
    particle3 = rbTorsions[index].particle3;
    particle4 = rbTorsions[index].particle4;
    c0 = rbTorsions[index].c[0];
    c1 = rbTorsions[index].c[1];
    c2 = rbTorsions[index].c[2];
    c3 = rbTorsions[index].c[3];
    c4 = rbTorsions[index].c[4];
    c5 = rbTorsions[index].c[5];
}

void RBTorsionForce::setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, double c0, double c1, double c2, double c3, double c4, double c5) {
    ASSERT_VALID_INDEX(index, rbTorsions);
    rbTorsions[index].particle1 = particle1;
    rbTorsions[index].particle2 = particle2;
    rbTorsions[index].particle3 = particle3;
    rbTorsions[index].particle4 = particle4;
    rbTorsions[index].c[0] = c0;
    rbTorsions[index].c[1] = c1;
    rbTorsions[index].c[2] = c2;
    rbTorsions[index].c[3] = c3;
    rbTorsions[index].c[4] = c4;
    rbTorsions[index].c[5] = c5;
}

ForceImpl* RBTorsionForce::createImpl() const {
    return new RBTorsionForceImpl(*this);
}

void RBTorsionForce::updateParametersInContext(Context& context) {
    dynamic_cast<RBTorsionForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void RBTorsionForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool RBTorsionForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
