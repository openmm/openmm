/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/GBVISoftcoreForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    GBVISoftcoreForce force;
    force.setNonbondedMethod(GBVISoftcoreForce::CutoffPeriodic);
    force.setBornRadiusScalingMethod(GBVISoftcoreForce::QuinticSpline);
    force.setCutoffDistance(2.0);
    force.setSoluteDielectric(5.1);
    force.setSolventDielectric(50.0);
    force.setQuinticLowerLimitFactor(0.1);
    force.setQuinticUpperBornRadiusLimit(5.1);
    //force.setTanhParameters(5.1, 4.2, 1.3);
    
    force.addParticle(1, 0.1, 0.01, 0.1);
    force.addParticle(0.5, 0.2, 0.02, 0.2);
    force.addParticle(-0.5, 0.3, 0.03, 0.3);
    force.addBond(0, 1, 2.0);
    force.addBond(3, 5, 1.2);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<GBVISoftcoreForce>(&force, "Force", buffer);
    GBVISoftcoreForce* copy = XmlSerializer::deserialize<GBVISoftcoreForce>(buffer);

    // Compare the two forces to see if they are identical.

    GBVISoftcoreForce& force2 = *copy;
    ASSERT_EQUAL(force.getNonbondedMethod(), force2.getNonbondedMethod());
    ASSERT_EQUAL(force.getCutoffDistance(), force2.getCutoffDistance());
    ASSERT_EQUAL(force.getSoluteDielectric(), force2.getSoluteDielectric());
    ASSERT_EQUAL(force.getSolventDielectric(), force2.getSolventDielectric());
    ASSERT_EQUAL(force.getQuinticLowerLimitFactor(), force2.getQuinticLowerLimitFactor());
    ASSERT_EQUAL(force.getQuinticUpperBornRadiusLimit(), force2.getQuinticUpperBornRadiusLimit());
    ASSERT_EQUAL(force.getBornRadiusScalingMethod(), force2.getBornRadiusScalingMethod());
/*
    double alpha1, beta1, gamma1;
    double alpha2, beta2, gamma2;
    force.getTanhParameters( alpha1, beta1, gamma1 );
    force2.getTanhParameters( alpha2, beta2, gamma2 );
    ASSERT_EQUAL(alpha1, alpha2);
    ASSERT_EQUAL(beta1, beta2);
    ASSERT_EQUAL(gamma1, gamma2);
*/
    ASSERT_EQUAL(force.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge1, radius1, scale1, bornRadiusScaleFactor1;
        double charge2, radius2, scale2, bornRadiusScaleFactor2;
        force.getParticleParameters(i, charge1, radius1, scale1, bornRadiusScaleFactor1);
        force2.getParticleParameters(i, charge2, radius2, scale2, bornRadiusScaleFactor2);
        ASSERT_EQUAL(charge1, charge2);
        ASSERT_EQUAL(radius1, radius2);
        ASSERT_EQUAL(scale1, scale2);
        ASSERT_EQUAL(bornRadiusScaleFactor1, bornRadiusScaleFactor2);
    }
    ASSERT_EQUAL(force.getNumBonds(), force2.getNumBonds());
    for (int i = 0; i < force.getNumBonds(); i++) {
        int a1, a2, b1, b2;
        double da, db;
        force.getBondParameters(i, a1, a2, da);
        force2.getBondParameters(i, b1, b2, db);
        ASSERT_EQUAL(a1, b1);
        ASSERT_EQUAL(a2, b2);
        ASSERT_EQUAL(da, db);
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
