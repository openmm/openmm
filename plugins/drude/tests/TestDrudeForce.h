/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2024 Stanford University and the Authors.      *
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
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/DrudeForce.h"
#include "SimTKOpenMMUtilities.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void validateForce(System& system, vector<Vec3>& positions, double expectedEnergy) {
    // Given a System containing a Drude force, check that its energy has the expected value.

    VerletIntegrator integ(1.0);
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-5);

    // Try moving each particle along each axis, and see if the energy changes by the correct amount.
    double offset = 1e-3;
    double TOL = 1e-5;
    try {
        if (platform.getPropertyValue(context, "Precision") != "double") {
            offset = 1e-2;
            TOL = 5e-4;
        }
    } catch(OpenMMException) {
        // The defaults above are for double precision, which is assumed in this case
    }
    for (int i = 0; i < system.getNumParticles(); i++)
        for (int j = 0; j < 3; j++) {
            vector<Vec3> offsetPos = positions;
            offsetPos[i][j] = positions[i][j]-offset;
            context.setPositions(offsetPos);
            double e1 = context.getState(State::Energy | State::Forces).getPotentialEnergy();
            offsetPos[i][j] = positions[i][j]+offset;
            context.setPositions(offsetPos);
            double e2 = context.getState(State::Energy | State::Forces).getPotentialEnergy();
            ASSERT_EQUAL_TOL(state.getForces()[i][j], (e1-e2)/(2*offset), TOL);
        }
}

void testSingleParticle() {
    const double k = ONE_4PI_EPS0*1.5;
    const double charge = 0.1;
    const double alpha = ONE_4PI_EPS0*charge*charge/k;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    DrudeForce* drude = new DrudeForce();
    drude->addParticle(1, 0, -1, -1, -1, charge, alpha, 1, 1);
    system.addForce(drude);
    ASSERT(!drude->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    validateForce(system, positions, 0.5*k*3*3);
}

void testAnisotropicParticle() {
    const double k = ONE_4PI_EPS0*1.5;
    const double charge = 0.1;
    const double alpha = ONE_4PI_EPS0*charge*charge/k;
    const double a1 = 0.8;
    const double a2 = 1.1;
    const double k1 = k/a1;
    const double k2 = k/a2;
    const double k3 = k/(3-a1-a2);
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    DrudeForce* drude = new DrudeForce();
    drude->addParticle(1, 0, 2, 3, 4, charge, alpha, a1, a2);
    system.addForce(drude);
    vector<Vec3> positions(5);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0.1, -0.5, 0.8);
    positions[2] = Vec3(0, 2, 0);
    positions[3] = Vec3(1, 2, 0);
    positions[4] = Vec3(1, 2, 3);
    validateForce(system, positions, 0.5*k1*0.5*0.5 + 0.5*k2*0.8*0.8 + 0.5*k3*0.1*0.1);
}

double computeScreening(double r, double thole, double alpha1, double alpha2) {
    double u = r*thole/pow(alpha1*alpha2, 1.0/6.0);
    return 1.0-(1.0+u/2)*exp(-u);
}

void testThole(bool periodic) {
    const double k = ONE_4PI_EPS0*1.5;
    const double charge = 0.1;
    const double alpha = ONE_4PI_EPS0*charge*charge/k;
    const double thole = 2.5;
    const double box = 2.0;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.setDefaultPeriodicBoxVectors(Vec3(box, 0, 0), Vec3(0, box, 0), Vec3(0, 0, box));
    DrudeForce* drude = new DrudeForce();
    drude->addParticle(1, 0, -1, -1, -1, charge, alpha, 1, 1);
    drude->addParticle(3, 2, -1, -1, -1, charge, alpha, 1, 1);
    drude->addScreenedPair(0, 1, thole);
    drude->setUsesPeriodicBoundaryConditions(periodic);
    system.addForce(drude);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, -0.5, 0);
    positions[2] = Vec3(1.1, 0, 0);
    positions[3] = Vec3(1.1, 0, 0.3);
    double energySpring1 = 0.5*k*0.5*0.5;
    double energySpring2 = 0.5*k*0.3*0.3;
    double energyDipole = 0.0;
    double q[] = {-charge, charge, -charge, charge};
    for (int i = 0; i < 2; i++)
        for (int j = 2; j < 4; j++) {
            Vec3 delta = positions[i]-positions[j];
            if (periodic) {
                delta[0] -= box*floor(delta[0]/box+0.5);
                delta[1] -= box*floor(delta[1]/box+0.5);
                delta[2] -= box*floor(delta[2]/box+0.5);
            }
            double r = sqrt(delta.dot(delta));
            energyDipole += ONE_4PI_EPS0*q[i]*q[j]*computeScreening(r, thole, alpha, alpha)/r;
        }
    validateForce(system, positions, energySpring1+energySpring2+energyDipole);
}

void testChangingParameters() {
    const double k = ONE_4PI_EPS0*1.5;
    const double charge = 0.1;
    const double alpha = ONE_4PI_EPS0*charge*charge/k;

    // Create the system.

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    DrudeForce* drude = new DrudeForce();
    drude->addParticle(1, 0, -1, -1, -1, charge, alpha, 1, 1);
    system.addForce(drude);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(2, 0, 0);

    // Check the energy.

    VerletIntegrator integ(1.0);
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(0.5*k*3*3, state.getPotentialEnergy(), 1e-5);

    // Modify the parameters.

    const double k2 = ONE_4PI_EPS0*2.2;
    const double charge2 = 0.3;
    const double alpha2 = ONE_4PI_EPS0*charge2*charge2/k2;
    drude->setParticleParameters(0, 1, 0, -1, -1, -1, charge2, alpha2, 1, 1);
    drude->updateParametersInContext(context);
    state = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(0.5*k2*3*3, state.getPotentialEnergy(), 1e-5);
}

void setupKernels(int argc, char* argv[]);
void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        setupKernels(argc, argv);
        testSingleParticle();
        testAnisotropicParticle();
        testThole(false);
        testThole(true);
        testChangingParameters();
        runPlatformTests();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
