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
#include "openmm/GBSAOBCForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    GBSAOBCForce force;
    force.setForceGroup(3);
    force.setNonbondedMethod(GBSAOBCForce::CutoffPeriodic);
    force.setCutoffDistance(2.0);
    force.setSoluteDielectric(5.1);
    force.setSolventDielectric(50.0);
    force.setSurfaceAreaEnergy(1.7);
    force.addParticle(1, 0.1, 0.01);
    force.addParticle(0.5, 0.2, 0.02);
    force.addParticle(-0.5, 0.3, 0.03);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<GBSAOBCForce>(&force, "Force", buffer);
    GBSAOBCForce* copy = XmlSerializer::deserialize<GBSAOBCForce>(buffer);

    // Compare the two forces to see if they are identical.

    GBSAOBCForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.getNonbondedMethod(), force2.getNonbondedMethod());
    ASSERT_EQUAL(force.getCutoffDistance(), force2.getCutoffDistance());
    ASSERT_EQUAL(force.getSoluteDielectric(), force2.getSoluteDielectric());
    ASSERT_EQUAL(force.getSolventDielectric(), force2.getSolventDielectric());
    ASSERT_EQUAL(force.getSurfaceAreaEnergy(), force2.getSurfaceAreaEnergy());
    ASSERT_EQUAL(force.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge1, radius1, scale1;
        double charge2, radius2, scale2;
        force.getParticleParameters(i, charge1, radius1, scale1);
        force2.getParticleParameters(i, charge2, radius2, scale2);
        ASSERT_EQUAL(charge1, charge2);
        ASSERT_EQUAL(radius1, radius2);
        ASSERT_EQUAL(scale1, scale2);
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
