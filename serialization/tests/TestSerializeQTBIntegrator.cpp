/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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
#include "openmm/QTBIntegrator.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>

#include <fstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    QTBIntegrator integ(310.0, 10.0, 0.002);
    integ.setSegmentLength(0.5);
    integ.setCutoffFrequency(600.0);
    integ.setDefaultAdaptationRate(0.02);
    integ.setRandomNumberSeed(10);
    integ.setIntegrationForceGroups(3);
    for (int i = 0; i < 5; i++)
        integ.setParticleType(i, i%3);
    for (int i = 0; i < 3; i++)
        integ.setTypeAdaptationRate(i, 0.1*i);
    stringstream ss;
    XmlSerializer::serialize<Integrator>(&integ, "QTBIntegrator", ss);
    QTBIntegrator *copy = dynamic_cast<QTBIntegrator*>(XmlSerializer::deserialize<Integrator>(ss));
    QTBIntegrator& integ2 = *copy;
    ASSERT_EQUAL(integ.getConstraintTolerance(), integ2.getConstraintTolerance());
    ASSERT_EQUAL(integ.getStepSize(), integ2.getStepSize());
    ASSERT_EQUAL(integ.getTemperature(), integ2.getTemperature());
    ASSERT_EQUAL(integ.getFriction(), integ2.getFriction());
    ASSERT_EQUAL(integ.getSegmentLength(), integ2.getSegmentLength());
    ASSERT_EQUAL(integ.getCutoffFrequency(), integ2.getCutoffFrequency());
    ASSERT_EQUAL(integ.getDefaultAdaptationRate(), integ2.getDefaultAdaptationRate());
    ASSERT_EQUAL(integ.getRandomNumberSeed(), integ2.getRandomNumberSeed());
    ASSERT_EQUAL(integ.getIntegrationForceGroups(), integ2.getIntegrationForceGroups());
    ASSERT_EQUAL_CONTAINERS(integ.getParticleTypes(), integ2.getParticleTypes());
    ASSERT_EQUAL_CONTAINERS(integ.getTypeAdaptationRates(), integ2.getTypeAdaptationRates());
    delete copy;
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
