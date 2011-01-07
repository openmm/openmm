/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
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
#include "openmm/AmoebaOutOfPlaneBendForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    AmoebaOutOfPlaneBendForce force;
    force.setAmoebaGlobalOutOfPlaneBendCubic( 12.3 );
    force.setAmoebaGlobalOutOfPlaneBendQuartic( 98.7 );
    force.setAmoebaGlobalOutOfPlaneBendPentic( 91.7 );
    force.setAmoebaGlobalOutOfPlaneBendSextic( 93.7 );
    force.addOutOfPlaneBend(0, 1, 3, 4, 2.0);
    force.addOutOfPlaneBend(0, 2, 3, 5, 2.1);
    force.addOutOfPlaneBend(2, 3, 5, 6, 2.2);
    force.addOutOfPlaneBend(5, 1, 8, 8, 2.3);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaOutOfPlaneBendForce>(&force, "Force", buffer);
    AmoebaOutOfPlaneBendForce* copy = XmlSerializer::deserialize<AmoebaOutOfPlaneBendForce>(buffer);

    // Compare the two forces to see if they are identical.  
    AmoebaOutOfPlaneBendForce& force2 = *copy;
    ASSERT_EQUAL(force.getAmoebaGlobalOutOfPlaneBendCubic(), force2.getAmoebaGlobalOutOfPlaneBendCubic());
    ASSERT_EQUAL(force.getAmoebaGlobalOutOfPlaneBendQuartic(), force2.getAmoebaGlobalOutOfPlaneBendQuartic());
    ASSERT_EQUAL(force.getAmoebaGlobalOutOfPlaneBendPentic(), force2.getAmoebaGlobalOutOfPlaneBendPentic());
    ASSERT_EQUAL(force.getAmoebaGlobalOutOfPlaneBendSextic(), force2.getAmoebaGlobalOutOfPlaneBendSextic());
    ASSERT_EQUAL(force.getNumOutOfPlaneBends(), force2.getNumOutOfPlaneBends());
    for (int i = 0; i < force.getNumOutOfPlaneBends(); i++) {
        int a1, a2, a3, a4, b1, b2, b3, b4;
        double ka, kb;
        force.getOutOfPlaneBendParameters(i, a1, a2, a3, a4, ka);
        force2.getOutOfPlaneBendParameters(i, b1, b2, b3, b4, kb);
        ASSERT_EQUAL(a1, b1);
        ASSERT_EQUAL(a2, b2);
        ASSERT_EQUAL(a3, b3);
        ASSERT_EQUAL(a4, b4);
        ASSERT_EQUAL(ka, kb);
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

