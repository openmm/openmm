/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "../../../tests/AssertionUtilities.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/System.h"
#include "openmm/serialization/HarmonicBondForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/SerializationProxy.h"
#include "openmm/serialization/SystemProxy.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a System.

    System system;
    for (int i = 0; i < 5; i++)
        system.addParticle(0.1*i+1);
    system.addConstraint(0, 1, 3.0);
    system.addConstraint(1, 2, 2.5);
    system.addConstraint(4, 1, 1.001);
    system.setDefaultPeriodicBoxVectors(Vec3(5, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 1.5));
    system.addForce(new HarmonicBondForce());

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<System>(&system, "System", buffer);
    System* copy = XmlSerializer::deserialize<System>(buffer);

    // Compare the two systems to see if they are identical.

    System& system2 = *copy;
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
    ASSERT_EQUAL(system.getNumForces(), system2.getNumForces());
    for (int i = 0; i < system.getNumForces(); i++)
        ASSERT(typeid(system.getForce(i)) == typeid(system2.getForce(i)))
}

int main() {
    try {
        SerializationProxy::registerProxy(typeid(System), new SystemProxy());
        SerializationProxy::registerProxy(typeid(HarmonicBondForce), new HarmonicBondForceProxy());
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}


