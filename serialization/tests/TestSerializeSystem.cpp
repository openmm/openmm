/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2017 Stanford University and the Authors.      *
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void compareSystems(System& system, System& system2) {
    ASSERT_EQUAL(system.getNumParticles(), system2.getNumParticles());
    for (int i = 0; i < system.getNumParticles(); i++)
        ASSERT_EQUAL(system.getParticleMass(i), system2.getParticleMass(i));
    ASSERT_EQUAL(system.getNumConstraints(), system2.getNumConstraints());
    for (int i = 0; i < system.getNumConstraints(); i++) {
        int p1, p2, p3, p4;
        double d1, d2;
        system.getConstraintParameters(i, p1, p2, d1);
        system2.getConstraintParameters(i, p3, p4, d2);
        ASSERT_EQUAL(p1, p3);
        ASSERT_EQUAL(p2, p4);
        ASSERT_EQUAL(d1, d2);
    }
    Vec3 a, b, c;
    Vec3 a2, b2, c2;
    system.getDefaultPeriodicBoxVectors(a, b, c);
    system2.getDefaultPeriodicBoxVectors(a2, b2, c2);
    ASSERT_EQUAL_VEC(a, a2, 0);
    ASSERT_EQUAL_VEC(b, b2, 0);
    ASSERT_EQUAL_VEC(c, c2, 0);
    for (int i = 0; i < system.getNumParticles(); i++)
        ASSERT_EQUAL(system.isVirtualSite(i), system2.isVirtualSite(i));
    const TwoParticleAverageSite& site5 = dynamic_cast<const TwoParticleAverageSite&>(system2.getVirtualSite(5));
    ASSERT_EQUAL(0, site5.getParticle(0));
    ASSERT_EQUAL(1, site5.getParticle(1));
    ASSERT_EQUAL(0.3, site5.getWeight(0));
    ASSERT_EQUAL(0.7, site5.getWeight(1));
    const ThreeParticleAverageSite& site6 = dynamic_cast<const ThreeParticleAverageSite&>(system2.getVirtualSite(6));
    ASSERT_EQUAL(2, site6.getParticle(0));
    ASSERT_EQUAL(4, site6.getParticle(1));
    ASSERT_EQUAL(3, site6.getParticle(2));
    ASSERT_EQUAL(0.5, site6.getWeight(0));
    ASSERT_EQUAL(0.2, site6.getWeight(1));
    ASSERT_EQUAL(0.3, site6.getWeight(2));
    const OutOfPlaneSite& site7 = dynamic_cast<const OutOfPlaneSite&>(system2.getVirtualSite(7));
    ASSERT_EQUAL(0, site7.getParticle(0));
    ASSERT_EQUAL(3, site7.getParticle(1));
    ASSERT_EQUAL(1, site7.getParticle(2));
    ASSERT_EQUAL(0.1, site7.getWeight12());
    ASSERT_EQUAL(0.2, site7.getWeight13());
    ASSERT_EQUAL(0.5, site7.getWeightCross());
    const LocalCoordinatesSite& site8 = dynamic_cast<const LocalCoordinatesSite&>(system2.getVirtualSite(8));
    ASSERT_EQUAL(4, site8.getNumParticles());
    ASSERT_EQUAL(4, site8.getParticle(0));
    ASSERT_EQUAL(3, site8.getParticle(1));
    ASSERT_EQUAL(2, site8.getParticle(2));
    ASSERT_EQUAL(1, site8.getParticle(3));
    ASSERT_EQUAL(Vec3(-0.5, 1.0, 1.5), site8.getLocalPosition());
    vector<double> wo, wx, wy;
    site8.getOriginWeights(wo);
    site8.getXWeights(wx);
    site8.getYWeights(wy);
    vector<double> woExpected = {0.1, 0.2, 0.3, 0.4};
    vector<double> wxExpected = {-1.0, 0.4, 0.4, 0.2};
    vector<double> wyExpected = {0.3, 0.7, 0.0, -1.0};
    ASSERT_EQUAL_CONTAINERS(woExpected, wo);
    ASSERT_EQUAL_CONTAINERS(wxExpected, wx);
    ASSERT_EQUAL_CONTAINERS(wyExpected, wy);
    ASSERT_EQUAL(system.getNumForces(), system2.getNumForces());
    for (int i = 0; i < system.getNumForces(); i++)
        ASSERT(typeid(system.getForce(i)) == typeid(system2.getForce(i)))
}

void testSerialization() {
    // Create a System.

    System system;
    for (int i = 0; i < 5; i++)
        system.addParticle(0.1*i+1);
    for (int i = 0; i < 5; i++)
        system.addParticle(0.0);
    system.addConstraint(0, 1, 3.0);
    system.addConstraint(1, 2, 2.5);
    system.addConstraint(4, 1, 1.001);
    system.setDefaultPeriodicBoxVectors(Vec3(5, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 1.5));
    system.setVirtualSite(5, new TwoParticleAverageSite(0, 1, 0.3, 0.7));
    system.setVirtualSite(6, new ThreeParticleAverageSite(2, 4, 3, 0.5, 0.2, 0.3));
    system.setVirtualSite(7, new OutOfPlaneSite(0, 3, 1, 0.1, 0.2, 0.5));
    system.setVirtualSite(8, new LocalCoordinatesSite({4, 3, 2, 1}, {0.1, 0.2, 0.3, 0.4}, {-1.0, 0.4, 0.4, 0.2}, {0.3, 0.7, 0.0, -1.0}, Vec3(-0.5, 1.0, 1.5)));
    system.addForce(new HarmonicBondForce());

    // Serialize and then deserialize it, then make sure the systems are identical.

    stringstream buffer;
    XmlSerializer::serialize<System>(&system, "System", buffer);
    System* copy = XmlSerializer::deserialize<System>(buffer);
    compareSystems(system, *copy);
    delete copy;

    // Now do the same thing but by calling clone().

    copy = XmlSerializer::clone(system);
    compareSystems(system, *copy);
    delete copy;
}

int main() {
    try {
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}


