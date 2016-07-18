/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
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

#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/AmoebaStretchBendForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaSerializationProxies();

void testSerialization() {
    // Create a Force.

    AmoebaStretchBendForce force1;
    force1.setForceGroup(3);
    force1.addStretchBend(0, 1, 3, 1.0, 1.2, 150.1, 83.2, 100.);
    force1.addStretchBend(2, 4, 4, 1.1, 2.2, 180.1, 89.2, 100.);
    force1.addStretchBend(5, 0, 1, 3.1, 8.2, 140.1, 98.2, 100.);
    force1.setUsesPeriodicBoundaryConditions(true);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaStretchBendForce>(&force1, "Force", buffer);
    AmoebaStretchBendForce* copy = XmlSerializer::deserialize<AmoebaStretchBendForce>(buffer);

    // Compare the two forces to see if they are identical.  
    AmoebaStretchBendForce& force2 = *copy;
    ASSERT_EQUAL(force1.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force1.usesPeriodicBoundaryConditions(), force2.usesPeriodicBoundaryConditions());
    ASSERT_EQUAL(force1.getNumStretchBends(), force2.getNumStretchBends());
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumStretchBends()); ii++) {
        int p11, p12, p13;
        int p21, p22, p23;
        double dAB1, dAB2;
        double dCB1, dCB2;
        double angle1, angle2;
        double k11, k12, k21, k22;

        force1.getStretchBendParameters(ii, p11, p12, p13, dAB1, dCB1, angle1, k11, k12);
        force2.getStretchBendParameters(ii, p21, p22, p23, dAB2, dCB2, angle2, k21, k22);

        ASSERT_EQUAL(p11, p21);
        ASSERT_EQUAL(p12, p22);
        ASSERT_EQUAL(p13, p23);
        ASSERT_EQUAL(dAB1, dAB2);
        ASSERT_EQUAL(dCB1, dCB2);
        ASSERT_EQUAL(angle1, angle2);
        ASSERT_EQUAL(k11, k21);
        ASSERT_EQUAL(k12, k22);
    }
}

int main() {
    try {
        registerAmoebaSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

