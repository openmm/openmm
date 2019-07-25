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
#include "openmm/GayBerneForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    GayBerneForce force;
    force.setForceGroup(3);
    force.setNonbondedMethod(GayBerneForce::CutoffPeriodic);
    force.setSwitchingDistance(1.5);
    force.setUseSwitchingFunction(true);
    force.setCutoffDistance(2.0);
    force.addParticle(0.1, 0.01, -1, -1, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0);
    force.addParticle(0.2, 0.02, -1, -1, 0.7, 0.7, 0.7, 1.2, 1.2, 1.2);
    force.addParticle(0.3, 0.03, 1, 0, 0.8, 0.9, 1.0, 0.5, 0.6, 0.7);
    force.addException(0, 1, 0.5, 0.1);
    force.addException(1, 2, 0.4, 0.2);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<GayBerneForce>(&force, "Force", buffer);
    GayBerneForce* copy = XmlSerializer::deserialize<GayBerneForce>(buffer);

    // Compare the two forces to see if they are identical.

    GayBerneForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.getNonbondedMethod(), force2.getNonbondedMethod());
    ASSERT_EQUAL(force.getSwitchingDistance(), force2.getSwitchingDistance());
    ASSERT_EQUAL(force.getUseSwitchingFunction(), force2.getUseSwitchingFunction());
    ASSERT_EQUAL(force.getCutoffDistance(), force2.getCutoffDistance());
    ASSERT_EQUAL(force.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); i++) {
        double sigma1, epsilon1, sx1, sy1, sz1, ex1, ey1, ez1;
        double sigma2, epsilon2, sx2, sy2, sz2, ex2, ey2, ez2;
        int xparticle1, yparticle1, xparticle2, yparticle2;
        force.getParticleParameters(i, sigma1, epsilon1, xparticle1, yparticle1, sx1, sy1, sz1, ex1, ey1, ez1);
        force2.getParticleParameters(i, sigma2, epsilon2, xparticle2, yparticle2, sx2, sy2, sz2, ex2, ey2, ez2);
        ASSERT_EQUAL(sigma1, sigma2);
        ASSERT_EQUAL(epsilon1, epsilon2);
        ASSERT_EQUAL(xparticle1, xparticle1);
        ASSERT_EQUAL(xparticle2, xparticle2);
        ASSERT_EQUAL(sx1, sx2);
        ASSERT_EQUAL(sy1, sy2);
        ASSERT_EQUAL(sz1, sz2);
        ASSERT_EQUAL(ex1, ex2);
        ASSERT_EQUAL(ey1, ey2);
        ASSERT_EQUAL(ez1, ez2);
    }
    ASSERT_EQUAL(force.getNumExceptions(), force2.getNumExceptions());
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int a1, a2, b1, b2;
        double sigma1, epsilon1;
        double sigma2, epsilon2;
        force.getExceptionParameters(i, a1, b1, sigma1, epsilon1);
        force2.getExceptionParameters(i, a2, b2, sigma2, epsilon2);
        ASSERT_EQUAL(a1, a2);
        ASSERT_EQUAL(b1, b2);
        ASSERT_EQUAL(sigma1, sigma2);
        ASSERT_EQUAL(epsilon1, epsilon2);
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
