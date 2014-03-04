/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "openmm/internal/AssertionUtilities.h"

using namespace OpenMM;

System::System() {
    periodicBoxVectors[0] = Vec3(2, 0, 0);
    periodicBoxVectors[1] = Vec3(0, 2, 0);
    periodicBoxVectors[2] = Vec3(0, 0, 2);
}

System::~System() {
    for (int i = 0; i < (int) forces.size(); ++i)
        delete forces[i];
    for (int i = 0; i < (int) virtualSites.size(); ++i)
        delete virtualSites[i];
}

double System::getParticleMass(int index) const {
    ASSERT_VALID_INDEX(index, masses);
    return masses[index];
}

void System::setParticleMass(int index, double mass) {
    ASSERT_VALID_INDEX(index, masses);
    masses[index] = mass;
}


void System::setVirtualSite(int index, VirtualSite* virtualSite) {
    if (index >= (int) virtualSites.size())
        virtualSites.resize(getNumParticles(), NULL);
    virtualSites[index] = virtualSite;
}

const VirtualSite& System::getVirtualSite(int index) const {
    if (index >= (int) virtualSites.size() || virtualSites[index] == NULL)
        throw OpenMMException("This particle is not a virtual site");
    return *virtualSites[index];
}

int System::addConstraint(int particle1, int particle2, double distance) {
    constraints.push_back(ConstraintInfo(particle1, particle2, distance));
    return constraints.size()-1;
}

void System::getConstraintParameters(int index, int& particle1, int& particle2, double& distance) const {
    ASSERT_VALID_INDEX(index, constraints);
    particle1 = constraints[index].particle1;
    particle2 = constraints[index].particle2;
    distance = constraints[index].distance;
}

void System::setConstraintParameters(int index, int particle1, int particle2, double distance) {
    ASSERT_VALID_INDEX(index, constraints);
    constraints[index].particle1 = particle1;
    constraints[index].particle2 = particle2;
    constraints[index].distance = distance;
}

const Force& System::getForce(int index) const {
    ASSERT_VALID_INDEX(index, forces);
    return *forces[index];
}

Force& System::getForce(int index) {
    ASSERT_VALID_INDEX(index, forces);
    return *forces[index];
}

void System::getDefaultPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const {
    a = periodicBoxVectors[0];
    b = periodicBoxVectors[1];
    c = periodicBoxVectors[2];
}

void System::setDefaultPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) {
    if (a[1] != 0.0 || a[2] != 0.0)
        throw OpenMMException("First periodic box vector must be parallel to x.");
    if (b[0] != 0.0 || b[2] != 0.0)
        throw OpenMMException("Second periodic box vector must be parallel to y.");
    if (c[0] != 0.0 || c[1] != 0.0)
        throw OpenMMException("Third periodic box vector must be parallel to z.");
    periodicBoxVectors[0] = a;
    periodicBoxVectors[1] = b;
    periodicBoxVectors[2] = c;
}
