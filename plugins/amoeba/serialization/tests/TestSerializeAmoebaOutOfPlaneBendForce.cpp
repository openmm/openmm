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
#include "openmm/AmoebaOutOfPlaneBendForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaSerializationProxies();

void testSerialization() {
    // Create a Force.

    AmoebaOutOfPlaneBendForce force1;

    force1.setForceGroup(3);
    force1.setAmoebaGlobalOutOfPlaneBendCubic(12.3);
    force1.setAmoebaGlobalOutOfPlaneBendQuartic(98.7);
    force1.setAmoebaGlobalOutOfPlaneBendPentic(91.7);
    force1.setAmoebaGlobalOutOfPlaneBendSextic(93.7);

    force1.addOutOfPlaneBend(0, 1, 3, 4, 2.0);
    force1.addOutOfPlaneBend(0, 2, 3, 5, 2.1);
    force1.addOutOfPlaneBend(2, 3, 5, 6, 2.2);
    force1.addOutOfPlaneBend(5, 1, 8, 8, 2.3);
    force1.setUsesPeriodicBoundaryConditions(true);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaOutOfPlaneBendForce>(&force1, "Force", buffer);
    AmoebaOutOfPlaneBendForce* copy = XmlSerializer::deserialize<AmoebaOutOfPlaneBendForce>(buffer);

    // Compare the two forces to see if they are identical.  
    AmoebaOutOfPlaneBendForce& force2 = *copy;

    ASSERT_EQUAL(force1.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force1.usesPeriodicBoundaryConditions(), force2.usesPeriodicBoundaryConditions());
    ASSERT_EQUAL(force1.getAmoebaGlobalOutOfPlaneBendCubic(), force2.getAmoebaGlobalOutOfPlaneBendCubic());
    ASSERT_EQUAL(force1.getAmoebaGlobalOutOfPlaneBendQuartic(), force2.getAmoebaGlobalOutOfPlaneBendQuartic());
    ASSERT_EQUAL(force1.getAmoebaGlobalOutOfPlaneBendPentic(), force2.getAmoebaGlobalOutOfPlaneBendPentic());
    ASSERT_EQUAL(force1.getAmoebaGlobalOutOfPlaneBendSextic(), force2.getAmoebaGlobalOutOfPlaneBendSextic());
    ASSERT_EQUAL(force1.getNumOutOfPlaneBends(), force2.getNumOutOfPlaneBends());

    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumOutOfPlaneBends()); ii++) {
        int a1, a2, a3, a4, b1, b2, b3, b4;
        double ka, kb;
        force1.getOutOfPlaneBendParameters(ii, a1, a2, a3, a4, ka);
        force2.getOutOfPlaneBendParameters(ii, b1, b2, b3, b4, kb);
        ASSERT_EQUAL(a1, b1);
        ASSERT_EQUAL(a2, b2);
        ASSERT_EQUAL(a3, b3);
        ASSERT_EQUAL(a4, b4);
        ASSERT_EQUAL(ka, kb);
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

