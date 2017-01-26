/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2014 Stanford University and the Authors.      *
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
#include "openmm/NonbondedForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    NonbondedForce force;
    force.setForceGroup(3);
    force.setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    force.setSwitchingDistance(1.5);
    force.setUseSwitchingFunction(true);
    force.setCutoffDistance(2.0);
    force.setEwaldErrorTolerance(1e-3);
    force.setReactionFieldDielectric(50.0);
    force.setUseDispersionCorrection(false);
    double alpha = 0.5;
    int nx = 3, ny = 5, nz = 7;
    force.setPMEParameters(alpha, nx, ny, nz);
    double dalpha = 0.8;
    int dnx = 4, dny = 6, dnz = 7;
    force.setLJPMEParameters(dalpha, dnx, dny, dnz);
    force.addParticle(1, 0.1, 0.01);
    force.addParticle(0.5, 0.2, 0.02);
    force.addParticle(-0.5, 0.3, 0.03);
    force.addException(0, 1, 2, 0.5, 0.1);
    force.addException(1, 2, 0.2, 0.4, 0.2);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<NonbondedForce>(&force, "Force", buffer);
    NonbondedForce* copy = XmlSerializer::deserialize<NonbondedForce>(buffer);

    // Compare the two forces to see if they are identical.

    NonbondedForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.getNonbondedMethod(), force2.getNonbondedMethod());
    ASSERT_EQUAL(force.getSwitchingDistance(), force2.getSwitchingDistance());
    ASSERT_EQUAL(force.getUseSwitchingFunction(), force2.getUseSwitchingFunction());
    ASSERT_EQUAL(force.getCutoffDistance(), force2.getCutoffDistance());
    ASSERT_EQUAL(force.getEwaldErrorTolerance(), force2.getEwaldErrorTolerance());
    ASSERT_EQUAL(force.getReactionFieldDielectric(), force2.getReactionFieldDielectric());
    ASSERT_EQUAL(force.getUseDispersionCorrection(), force2.getUseDispersionCorrection());
    ASSERT_EQUAL(force.getNumParticles(), force2.getNumParticles());
    double alpha2;
    int nx2, ny2, nz2;
    force2.getPMEParameters(alpha2, nx2, ny2, nz2);
    ASSERT_EQUAL(alpha, alpha2);
    ASSERT_EQUAL(nx, nx2);
    ASSERT_EQUAL(ny, ny2);
    ASSERT_EQUAL(nz, nz2);    
    double dalpha2;
    int dnx2, dny2, dnz2;
    force2.getLJPMEParameters(dalpha2, dnx2, dny2, dnz2);
    ASSERT_EQUAL(dalpha, dalpha2);
    ASSERT_EQUAL(dnx, dnx2);
    ASSERT_EQUAL(dny, dny2);
    ASSERT_EQUAL(dnz, dnz2);    
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge1, sigma1, epsilon1;
        double charge2, sigma2, epsilon2;
        force.getParticleParameters(i, charge1, sigma1, epsilon1);
        force2.getParticleParameters(i, charge2, sigma2, epsilon2);
        ASSERT_EQUAL(charge1, charge2);
        ASSERT_EQUAL(sigma1, sigma2);
        ASSERT_EQUAL(epsilon1, epsilon2);
    }
    ASSERT_EQUAL(force.getNumExceptions(), force2.getNumExceptions());
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int a1, a2, b1, b2;
        double charge1, sigma1, epsilon1;
        double charge2, sigma2, epsilon2;
        force.getExceptionParameters(i, a1, b1, charge1, sigma1, epsilon1);
        force2.getExceptionParameters(i, a2, b2, charge2, sigma2, epsilon2);
        ASSERT_EQUAL(a1, a2);
        ASSERT_EQUAL(b1, b2);
        ASSERT_EQUAL(charge1, charge2);
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
