/* -------------------------------------------------------------------------- *
 *                                OpenMMDrude                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/DrudeForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

extern "C" void registerDrudeSerializationProxies();

void testSerialization() {
    // Create a Force.

    DrudeForce force1;
    force1.addParticle(0, 1, 2, 3, 4, 0.5, 1.0, 1.5, 2.0);
    force1.addParticle(2, 3, 7, 8, 9, 0.1, 1e-4, 1.0, 0.9);
    force1.addParticle(5, 6, -1, -1, -1, 0.2, 0.1, 1.0, 1.0);
    force1.addScreenedPair(0, 1, 2.6);
    force1.addScreenedPair(1, 2, 2.2);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<DrudeForce>(&force1, "Force", buffer);
    DrudeForce* copy = XmlSerializer::deserialize<DrudeForce>(buffer);

    // Compare the two forces to see if they are identical.  
    DrudeForce& force2 = *copy;
    ASSERT_EQUAL(force1.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < (int) force1.getNumParticles(); i++) {
        int a1, a2, a3, a4, a5, b1, b2, b3, b4, b5;
        double charge1, charge2;
        double polar1, polar2;
        double aa12, ba12, aa34, ba34;
        force1.getParticleParameters(i, a1, a2, a3, a4, a5, charge1, polar1, aa12, aa34);
        force2.getParticleParameters(i, b1, b2, b3, b4, b5, charge2, polar2, ba12, ba34);
        ASSERT_EQUAL(a1, b1);
        ASSERT_EQUAL(a2, b2);
        ASSERT_EQUAL(a3, b3);
        ASSERT_EQUAL(a4, b4);
        ASSERT_EQUAL(a5, b5);
        ASSERT_EQUAL(charge1, charge2);
        ASSERT_EQUAL(polar1, polar2);
        ASSERT_EQUAL(aa12, ba12);
        ASSERT_EQUAL(aa34, ba34);
    }
    ASSERT_EQUAL(force1.getNumScreenedPairs(), force2.getNumScreenedPairs());
    for (int i = 0; i < (int) force1.getNumScreenedPairs(); i++) {
        int a1, a2, b1, b2;
        double thole1, thole2;
        force1.getScreenedPairParameters(i, a1, a2, thole1);
        force1.getScreenedPairParameters(i, b1, b2, thole2);
        ASSERT_EQUAL(a1, b1);
        ASSERT_EQUAL(a2, b2);
        ASSERT_EQUAL(thole1, thole2);
    }
}

int main() {
    try {
        registerDrudeSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

