/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
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
#include "openmm/CMAPTorsionForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    CMAPTorsionForce force;
    force.setForceGroup(3);
    vector<double> map1(9);
    for (int i = 0; i < 9; i++)
        map1[i] = 0.1*i;
    force.addMap(3, map1);
    vector<double> map2(16);
    for (int i = 0; i < 16; i++)
        map2[i] = 0.2*i;
    force.addMap(4, map2);
    force.addTorsion(0, 0, 1, 2, 3, 2, 3, 4, 5);
    force.addTorsion(0, 0, 2, 3, 4, 5, 6, 7, 8);
    force.addTorsion(1, 2, 3, 4, 7, 1, 2, 3, 4);
    force.addTorsion(1, 5, 1, 2, 3, 2, 3, 4, 8);
    force.setUsesPeriodicBoundaryConditions(true);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<CMAPTorsionForce>(&force, "Force", buffer);
    CMAPTorsionForce* copy = XmlSerializer::deserialize<CMAPTorsionForce>(buffer);

    // Compare the two forces to see if they are identical.

    CMAPTorsionForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.usesPeriodicBoundaryConditions(), force2.usesPeriodicBoundaryConditions());
    ASSERT_EQUAL(force.getNumMaps(), force2.getNumMaps());
    for (int i = 0; i < force.getNumMaps(); i++) {
        int size1, size2;
        vector<double> energy1, energy2;
        force.getMapParameters(i, size1, energy1);
        force2.getMapParameters(i, size2, energy2);
        ASSERT_EQUAL(size1, size2);
        ASSERT_EQUAL(energy1.size(), energy2.size());
        for (int j = 0; j < (int) energy1.size(); j++)
            ASSERT_EQUAL(energy1[j], energy2[j]);
    }
    ASSERT_EQUAL(force.getNumTorsions(), force2.getNumTorsions());
    for (int i = 0; i < force.getNumTorsions(); i++) {
        int map1, map2;
        int a11, a21, a31, a41, a12, a22, a32, a42;
        int b11, b21, b31, b41, b12, b22, b32, b42;
        force.getTorsionParameters(i, map1, a11, a21, a31, a41, b11, b21, b31, b41);
        force2.getTorsionParameters(i, map2, a12, a22, a32, a42, b12, b22, b32, b42);
        ASSERT_EQUAL(map1, map2);
        ASSERT_EQUAL(a11, a12);
        ASSERT_EQUAL(a21, a22);
        ASSERT_EQUAL(a31, a32);
        ASSERT_EQUAL(a41, a42);
        ASSERT_EQUAL(b11, b12);
        ASSERT_EQUAL(b21, b22);
        ASSERT_EQUAL(b31, b32);
        ASSERT_EQUAL(b41, b42);
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
