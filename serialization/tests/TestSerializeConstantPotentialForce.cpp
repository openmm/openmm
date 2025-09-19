/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2010-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Pretti                                        *
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
#include "openmm/ConstantPotentialForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    ConstantPotentialForce force;
    force.setForceGroup(3);
    force.setName("custom name");
    force.setCutoffDistance(2.0);
    force.setEwaldErrorTolerance(1e-3);
    force.setExceptionsUsePeriodicBoundaryConditions(true);
    double alpha = 0.5;
    int nx = 3, ny = 5, nz = 7;
    force.setPMEParameters(alpha, nx, ny, nz);
    force.setConstantPotentialMethod(ConstantPotentialForce::Matrix);
    force.setUsePreconditioner(false);
    force.setCGErrorTolerance(2e-8);
    force.setUseChargeConstraint(true);
    force.setChargeConstraintTarget(0.1);
    Vec3 externalField(1.0, 2.0, 3.0);
    force.setExternalField(externalField);
    force.addParticle(1);
    force.addParticle(0.0);
    force.addParticle(0.0);
    force.addParticle(0.0);
    force.addParticle(0.0);
    force.addParticle(0.0);
    force.addParticle(0.5);
    force.addParticle(-0.5);
    force.addException(0, 6, 2);
    force.addException(6, 7, 0.2);
    std::set<int> electrode1Particles{1, 2, 4};
    std::set<int> electrode2Particles{3, 5};
    force.addElectrode(electrode1Particles, 1.0, 0.5, 2.0);
    force.addElectrode(electrode2Particles, -0.5, 0.75, 1.5);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<ConstantPotentialForce>(&force, "Force", buffer);
    ConstantPotentialForce* copy = XmlSerializer::deserialize<ConstantPotentialForce>(buffer);

    // Compare the two forces to see if they are identical.

    ConstantPotentialForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.getName(), force2.getName());
    ASSERT_EQUAL(force.getCutoffDistance(), force2.getCutoffDistance());
    ASSERT_EQUAL(force.getEwaldErrorTolerance(), force2.getEwaldErrorTolerance());
    ASSERT_EQUAL(force.getExceptionsUsePeriodicBoundaryConditions(), force2.getExceptionsUsePeriodicBoundaryConditions());
    double alpha2;
    int nx2, ny2, nz2;
    force2.getPMEParameters(alpha2, nx2, ny2, nz2);
    ASSERT_EQUAL(alpha, alpha2);
    ASSERT_EQUAL(nx, nx2);
    ASSERT_EQUAL(ny, ny2);
    ASSERT_EQUAL(nz, nz2);
    ASSERT_EQUAL(force.getConstantPotentialMethod(), force2.getConstantPotentialMethod());
    ASSERT_EQUAL(force.getUsePreconditioner(), force2.getUsePreconditioner());
    ASSERT_EQUAL(force.getCGErrorTolerance(), force2.getCGErrorTolerance());
    ASSERT_EQUAL(force.getUseChargeConstraint(), force2.getUseChargeConstraint());
    ASSERT_EQUAL(force.getChargeConstraintTarget(), force2.getChargeConstraintTarget());
    Vec3 externalField2;
    force2.getExternalField(externalField2);
    ASSERT_EQUAL(externalField, externalField2);
    ASSERT_EQUAL(force.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge1, charge2;
        force.getParticleParameters(i, charge1);
        force2.getParticleParameters(i, charge2);
        ASSERT_EQUAL(charge1, charge2);
    }
    ASSERT_EQUAL(force.getNumExceptions(), force2.getNumExceptions());
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int a1, a2, b1, b2;
        double charge1, charge2;
        force.getExceptionParameters(i, a1, b1, charge1);
        force2.getExceptionParameters(i, a2, b2, charge2);
        ASSERT_EQUAL(a1, a2);
        ASSERT_EQUAL(b1, b2);
        ASSERT_EQUAL(charge1, charge2);
    }
    ASSERT_EQUAL(force.getNumElectrodes(), force2.getNumElectrodes());
    for (int i = 0; i < force.getNumElectrodes(); i++) {
        std::set<int> electrodeParticles1, electrodeParticles2;
        double potential1, potential2;
        double gaussianWidth1, gaussianWidth2;
        double thomasFermiScale1, thomasFermiScale2;
        force.getElectrodeParameters(i, electrodeParticles1, potential1, gaussianWidth1, thomasFermiScale1);
        force2.getElectrodeParameters(i, electrodeParticles2, potential2, gaussianWidth2, thomasFermiScale2);
        ASSERT_EQUAL_CONTAINERS(electrodeParticles1, electrodeParticles2);
        ASSERT_EQUAL(potential1, potential2);
        ASSERT_EQUAL(gaussianWidth1, gaussianWidth2);
        ASSERT_EQUAL(thomasFermiScale1, thomasFermiScale2);
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
