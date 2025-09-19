/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.            *
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
#include "openmm/OrientationRestraintForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    vector<Vec3> refPos;
    for (int i = 0; i < 10; i++)
        refPos.push_back(Vec3(i/5.0, i*1.2, i*i/3.5));
    vector<int> particles;
    for (int i = 0; i < 5; i++)
        particles.push_back(i*i);
    OrientationRestraintForce force(2.3, refPos, particles);
    force.setForceGroup(3);
    force.setName("custom name");

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<OrientationRestraintForce>(&force, "Force", buffer);
    OrientationRestraintForce* copy = XmlSerializer::deserialize<OrientationRestraintForce>(buffer);

    // Compare the two forces to see if they are identical.

    OrientationRestraintForce& force2 = *copy;
    ASSERT_EQUAL(force.getK(), force2.getK());
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.getName(), force2.getName());
    ASSERT_EQUAL(force.getReferencePositions().size(), force2.getReferencePositions().size());
    for (int i = 0; i < force.getReferencePositions().size(); i++)
        ASSERT_EQUAL_VEC(force.getReferencePositions()[i], force2.getReferencePositions()[i], 0.0);
    ASSERT_EQUAL(force.getParticles().size(), force2.getParticles().size());
    for (int i = 0; i < force.getParticles().size(); i++)
        ASSERT_EQUAL(force.getParticles()[i], force2.getParticles()[i]);
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
