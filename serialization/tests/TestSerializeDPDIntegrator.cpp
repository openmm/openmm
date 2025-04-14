/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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
#include "openmm/DPDIntegrator.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerializeDPDIntegrator() {
    // Create an integrator.

    DPDIntegrator integrator(120.0, 1.1, 1.2, 0.05);
    integrator.setRandomNumberSeed(4);
    integrator.setIntegrationForceGroups(3);
    integrator.setParticleType(1, 3);
    integrator.setParticleType(3, 8);
    integrator.addTypePair(0, 3, 1.5, 2.1);
    integrator.addTypePair(3, 8, 0.5, 0.9);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<DPDIntegrator>(&integrator, "Integrator", buffer);
    DPDIntegrator* copy = XmlSerializer::deserialize<DPDIntegrator>(buffer);

    // Compare the two integrators to see if they are identical.

    DPDIntegrator& integrator2 = *copy;
    ASSERT_EQUAL(integrator.getTemperature(), integrator2.getTemperature());
    ASSERT_EQUAL(integrator.getDefaultFriction(), integrator2.getDefaultFriction());
    ASSERT_EQUAL(integrator.getDefaultCutoff(), integrator2.getDefaultCutoff());
    ASSERT_EQUAL(integrator.getStepSize(), integrator2.getStepSize());
    ASSERT_EQUAL(integrator.getRandomNumberSeed(), integrator2.getRandomNumberSeed());
    ASSERT_EQUAL(integrator.getIntegrationForceGroups(), integrator2.getIntegrationForceGroups());
    for (int i = 0; i < 10; i++)
        ASSERT_EQUAL(integrator.getParticleType(i), integrator2.getParticleType(i));
    ASSERT_EQUAL_CONTAINERS(integrator.getParticleTypes(), integrator2.getParticleTypes());
    ASSERT_EQUAL(integrator.getNumTypePairs(), integrator2.getNumTypePairs());
    for (int i = 0; i < integrator.getNumTypePairs(); i++) {
        int a1, a2, b1, b2;
        double f1, f2, c1, c2;
        integrator.getTypePairParameters(i, a1, b1, f1, c1);
        integrator.getTypePairParameters(i, a2, b2, f2, c2);
        ASSERT_EQUAL(a1, a2);
        ASSERT_EQUAL(b1, b2);
        ASSERT_EQUAL(f1, f2);
        ASSERT_EQUAL(c1, c2);
    }
}

int main() {
    try {
        testSerializeDPDIntegrator();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
