/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2010-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Pretti                                        *
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
#include "openmm/LCPOForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    LCPOForce force;
    force.setForceGroup(3);
    force.setName("custom name");
    force.setUsesPeriodicBoundaryConditions(true);
    force.addParticle(1.0, 2.0, 3.0, 4.0, 5.0);
    force.addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
    force.addParticle(6.0, 7.0, 8.0, 9.0, 10.0);
    force.addParticle(1.0, 2.0, 3.0, 4.0, 5.0);
    force.addParticle(0.0, 0.0, 0.0, 0.0, 0.0);
    force.addParticle(6.0, 7.0, 8.0, 9.0, 10.0);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<LCPOForce>(&force, "Force", buffer);
    LCPOForce* copy = XmlSerializer::deserialize<LCPOForce>(buffer);

    // Compare the two forces to see if they are identical.

    LCPOForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.getName(), force2.getName());
    ASSERT_EQUAL(force.usesPeriodicBoundaryConditions(), force2.usesPeriodicBoundaryConditions());
    ASSERT_EQUAL(force.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); i++) {
        double radius1, radius2, p11, p12, p21, p22, p31, p32, p41, p42;
        force.getParticleParameters(i, radius1, p11, p21, p31, p41);
        force2.getParticleParameters(i, radius2, p12, p22, p32, p42);
        ASSERT_EQUAL(radius1, radius2);
        ASSERT_EQUAL(p11, p12);
        ASSERT_EQUAL(p21, p22);
        ASSERT_EQUAL(p31, p32);
        ASSERT_EQUAL(p41, p42);
    }
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
