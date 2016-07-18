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
#include "openmm/AmoebaWcaDispersionForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaSerializationProxies();

void testSerialization() {
    // Create a Force.

    AmoebaWcaDispersionForce force1;
    force1.setForceGroup(3);
    force1.setEpso(   1.0);
    force1.setEpsh(   1.1);
    force1.setRmino(  1.2);
    force1.setRminh(  1.3);
    force1.setAwater( 1.4);
    force1.setShctd(  1.5);
    force1.setDispoff(1.6);
    force1.setSlevy(  1.7);

    force1.addParticle(1.0, 2.0);
    force1.addParticle(1.1, 2.1);
    force1.addParticle(1.2, 2.2);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaWcaDispersionForce>(&force1, "Force", buffer);
    AmoebaWcaDispersionForce* copy = XmlSerializer::deserialize<AmoebaWcaDispersionForce>(buffer);

    // Compare the two forces to see if they are identical.  

    AmoebaWcaDispersionForce& force2 = *copy;

    ASSERT_EQUAL(force1.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force1.getEpso(),    force2.getEpso());
    ASSERT_EQUAL(force1.getEpsh(),    force2.getEpsh());
    ASSERT_EQUAL(force1.getRmino(),   force2.getRmino());
    ASSERT_EQUAL(force1.getRminh(),   force2.getRminh());
    ASSERT_EQUAL(force1.getAwater(),  force2.getAwater());
    ASSERT_EQUAL(force1.getShctd(),   force2.getShctd());
    ASSERT_EQUAL(force1.getDispoff(), force2.getDispoff());
    ASSERT_EQUAL(force1.getSlevy(),   force2.getSlevy());

    ASSERT_EQUAL(force1.getNumParticles(), force2.getNumParticles());
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumParticles()); ii++) {

        double radius1, epsilon1;
        double radius2, epsilon2;

        force1.getParticleParameters(ii, radius1, epsilon1);
        force2.getParticleParameters(ii, radius2, epsilon2);

        ASSERT_EQUAL(radius1,  radius2);
        ASSERT_EQUAL(epsilon1, epsilon2);
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

