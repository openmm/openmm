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
#include "openmm/AmoebaGeneralizedKirkwoodForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaSerializationProxies();

void testSerialization() {
    // Create a Force.

    AmoebaGeneralizedKirkwoodForce force1;
    force1.setForceGroup(3);
    force1.setSolventDielectric(  80.0);
    force1.setSoluteDielectric(   1.0);
    //force1.setDielectricOffset(   0.09);
    force1.setProbeRadius(        1.40);
    force1.setSurfaceAreaFactor(  0.888);
    force1.setIncludeCavityTerm(  1);

    force1.addParticle(1.0, 2.0, 0.9);
    force1.addParticle(-1.1,2.1, 0.8);
    force1.addParticle(0.1, 2.2, 0.7);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaGeneralizedKirkwoodForce>(&force1, "Force", buffer);
    AmoebaGeneralizedKirkwoodForce* copy = XmlSerializer::deserialize<AmoebaGeneralizedKirkwoodForce>(buffer);

    // Compare the two forces to see if they are identical.  
    AmoebaGeneralizedKirkwoodForce& force2 = *copy;
    ASSERT_EQUAL(force1.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force1.getSolventDielectric(),    force2.getSolventDielectric());
    ASSERT_EQUAL(force1.getSoluteDielectric(),     force2.getSoluteDielectric());
    //ASSERT_EQUAL(force1.getDielectricOffset(),     force2.getDielectricOffset());
    ASSERT_EQUAL(force1.getProbeRadius(),          force2.getProbeRadius());
    ASSERT_EQUAL(force1.getSurfaceAreaFactor(),    force2.getSurfaceAreaFactor());
    ASSERT_EQUAL(force1.getIncludeCavityTerm(),    force2.getIncludeCavityTerm());

    ASSERT_EQUAL(force1.getNumParticles(), force2.getNumParticles());
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumParticles()); ii++) {

        double radius1, charge1, scaleFactor1;
        double radius2, charge2, scaleFactor2;

        force1.getParticleParameters(ii, charge1, radius1, scaleFactor1);
        force2.getParticleParameters(ii, charge2, radius2, scaleFactor2);

        ASSERT_EQUAL(charge1,      charge2);
        ASSERT_EQUAL(radius1,      radius2);
        ASSERT_EQUAL(scaleFactor1, scaleFactor2);
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

