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
#include "openmm/MonteCarloAnisotropicBarostat.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    MonteCarloAnisotropicBarostat force(Vec3(15.1, 18.2, 19.3), 250.0, true, false, true, 14);
    force.setForceGroup(3);
    force.setRandomNumberSeed(3);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<MonteCarloAnisotropicBarostat>(&force, "Force", buffer);
    MonteCarloAnisotropicBarostat* copy = XmlSerializer::deserialize<MonteCarloAnisotropicBarostat>(buffer);

    // Compare the two forces to see if they are identical.

    MonteCarloAnisotropicBarostat& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL_VEC(force.getDefaultPressure(), force2.getDefaultPressure(), 0.0);
    ASSERT_EQUAL(force.getDefaultTemperature(), force2.getDefaultTemperature());
    ASSERT_EQUAL(force.getScaleX(), force2.getScaleX());
    ASSERT_EQUAL(force.getScaleY(), force2.getScaleY());
    ASSERT_EQUAL(force.getScaleZ(), force2.getScaleZ());
    ASSERT_EQUAL(force.getFrequency(), force2.getFrequency());
    ASSERT_EQUAL(force.getRandomNumberSeed(), force2.getRandomNumberSeed());
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
