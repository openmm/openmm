/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2010-2023 Stanford University and the Authors.      *
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
#include "openmm/ATMForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    ATMForce force("2*x*(u0+u1+y)");
    force.setForceGroup(3);
    force.setName("custom name");
    force.addGlobalParameter("x", 1.3);
    force.addGlobalParameter("y", 2.221);
    force.addEnergyParameterDerivative("y");
    HarmonicBondForce* v1 = new HarmonicBondForce();
    v1->addBond(2, 3, 5.2, 1.1);
    force.addForce(v1);
    HarmonicAngleForce* v2 = new HarmonicAngleForce();
    v2->addAngle(3, 11, 15, 0.4, 0.2);
    force.addForce(v2);
    force.addParticle(new ATMForce::FixedDisplacement(Vec3(1, 2, 3)));
    force.addParticle(new ATMForce::ParticleOffsetDisplacement(0, 1));

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<ATMForce>(&force, "Force", buffer);
    ATMForce* copy = XmlSerializer::deserialize<ATMForce>(buffer);

    // Compare the two forces to see if they are identical.

    ATMForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.getName(), force2.getName());
    ASSERT_EQUAL(force.getEnergyFunction(), force2.getEnergyFunction());
    ASSERT_EQUAL(force.getNumGlobalParameters(), force2.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        ASSERT_EQUAL(force.getGlobalParameterName(i), force2.getGlobalParameterName(i));
        ASSERT_EQUAL(force.getGlobalParameterDefaultValue(i), force2.getGlobalParameterDefaultValue(i));
    }
    ASSERT_EQUAL(force.getNumEnergyParameterDerivatives(), force2.getNumEnergyParameterDerivatives());
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
        ASSERT_EQUAL(force.getEnergyParameterDerivativeName(i), force2.getEnergyParameterDerivativeName(i));
    ASSERT_EQUAL(force.getNumForces(), force2.getNumForces());
    for (int i = 0; i < force.getNumForces(); i++) {
        stringstream buffer1, buffer2;
        XmlSerializer::serialize<Force>(&force.getForce(i), "Force", buffer1);
        XmlSerializer::serialize<Force>(&force2.getForce(i), "Force", buffer2);
        ASSERT_EQUAL(buffer1.str(), buffer2.str());
    }
    ASSERT_EQUAL(force.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); i++) {
        const ATMForce::CoordinateTransformation& transformation  = force.getParticleTransformation(i);
        const ATMForce::CoordinateTransformation& transformation2 = force2.getParticleTransformation(i);
        if (dynamic_cast<const ATMForce::FixedDisplacement*>(&transformation) != nullptr) {
            const ATMForce::FixedDisplacement* fd  = dynamic_cast<const ATMForce::FixedDisplacement*>(&transformation);
            const ATMForce::FixedDisplacement* fd2 = dynamic_cast<const ATMForce::FixedDisplacement*>(&transformation2);
            const Vec3 d1a = fd->getFixedDisplacement1();
            const Vec3 d0a = fd->getFixedDisplacement0();
            const Vec3 d1b = fd2->getFixedDisplacement1();
            const Vec3 d0b = fd2->getFixedDisplacement0();
            ASSERT_EQUAL_VEC(d1a, d1b, 0.0);
            ASSERT_EQUAL_VEC(d0a, d0b, 0.0);
        }
        else if (dynamic_cast<const ATMForce::ParticleOffsetDisplacement*>(&transformation) != nullptr) {
            const ATMForce::ParticleOffsetDisplacement* vd  = dynamic_cast<const ATMForce::ParticleOffsetDisplacement*>(&transformation);
            const ATMForce::ParticleOffsetDisplacement* vd2 = dynamic_cast<const ATMForce::ParticleOffsetDisplacement*>(&transformation2);
            int j1a = vd->getDestinationParticle1();
            int i1a = vd->getOriginParticle1();
            int j0a = vd->getDestinationParticle0();
            int i0a = vd->getOriginParticle0();
            int j1b = vd2->getDestinationParticle1();
            int i1b = vd2->getOriginParticle1();
            int j0b = vd2->getDestinationParticle0();
            int i0b = vd2->getOriginParticle0();
            ASSERT_EQUAL(j1a, j1b);
            ASSERT_EQUAL(i1a, i1b);
            ASSERT_EQUAL(j0a, j0b);
            ASSERT_EQUAL(i0a, i0b);
        }
        else {
            throwException(__FILE__, __LINE__, "Unknown CoordinateTransformation type");
        }

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
