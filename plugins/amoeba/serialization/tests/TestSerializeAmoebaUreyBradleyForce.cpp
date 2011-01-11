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
#include "openmm/AmoebaUreyBradleyForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    AmoebaUreyBradleyForce force1;
    force1.setAmoebaGlobalUreyBradleyCubic( 12.3 );
    force1.setAmoebaGlobalUreyBradleyQuartic( 98.7 );
    force1.addUreyBradley(0, 1, 1.0, 2.0);
    force1.addUreyBradley(0, 2, 2.0, 2.1);
    force1.addUreyBradley(2, 3, 3.0, 2.2);
    force1.addUreyBradley(5, 1, 4.0, 2.3);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaUreyBradleyForce>(&force1, "Force", buffer);
    AmoebaUreyBradleyForce* copy = XmlSerializer::deserialize<AmoebaUreyBradleyForce>(buffer);

    // Compare the two forces to see if they are identical.  
    AmoebaUreyBradleyForce& force2 = *copy;
    ASSERT_EQUAL(force1.getAmoebaGlobalUreyBradleyCubic(), force2.getAmoebaGlobalUreyBradleyCubic());
    ASSERT_EQUAL(force1.getAmoebaGlobalUreyBradleyQuartic(), force2.getAmoebaGlobalUreyBradleyQuartic());
    ASSERT_EQUAL(force1.getNumInteractions(), force2.getNumInteractions());
    for (unsigned int ii = 0; ii < force1.getNumInteractions(); ii++) {
        int a1, a2, b1, b2;
        double da, db, ka, kb;

        force1.getUreyBradleyParameters(ii, a1, a2, da, ka);
        force2.getUreyBradleyParameters(ii, b1, b2, db, kb);

        ASSERT_EQUAL(a1, b1);
        ASSERT_EQUAL(a2, b2);
        ASSERT_EQUAL(da, db);
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

