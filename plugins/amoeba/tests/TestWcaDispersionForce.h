/* -------------------------------------------------------------------------- *
 *                                   OpenMMAmoeba                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs                                                   *
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
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/AmoebaWcaDispersionForce.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol, testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};


using namespace OpenMM;
const double TOL = 1e-4;

void compareForcesEnergy(std::string &testName, double expectedEnergy, double energy,
                         const std::vector<Vec3> &expectedForces,
                         const std::vector<Vec3> &forces, double tolerance) {

    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        ASSERT_EQUAL_VEC_MOD(expectedForces[ii], forces[ii], tolerance, testName)
    }
    ASSERT_EQUAL_TOL_MOD(expectedEnergy, energy, tolerance, testName)
}

// Test Wca Dispersion
void testWcaDispersionAmmonia() {

    std::string testName = "testWcaDispersionAmmonia";
    int numberOfParticles = 8;

    // Create the system.

    System system;
    AmoebaWcaDispersionForce *amoebaWcaDispersionForce = new AmoebaWcaDispersionForce();

    // Convert kcal/mol to kJ/mol
    double epso = 0.1100e0 * 4.184e0;
    double epsh = 0.0135e0 * 4.184e0;
    // Convert A to nm.
    double rmino = 1.7025e-01;
    double rminh = 1.3275e-01;
    double dispoff = 1.056e-01;
    // Convert water number density from water / A^3 to water / nm^3.
    double awater = 0.033428e03;
    // No units.
    double slevy = 1.0e0;
    double shctd = 0.75e0;

    amoebaWcaDispersionForce->setEpso(epso);
    amoebaWcaDispersionForce->setEpsh(epsh);
    amoebaWcaDispersionForce->setRmino(rmino);
    amoebaWcaDispersionForce->setRminh(rminh);
    amoebaWcaDispersionForce->setDispoff(dispoff);
    amoebaWcaDispersionForce->setAwater(awater);
    amoebaWcaDispersionForce->setSlevy(slevy);
    amoebaWcaDispersionForce->setShctd(shctd);

    // Amoeba 2009
    // vdwtype                 BUFFERED-14-7
    // radiusrule              CUBIC-MEAN
    // radiustype              R-MIN
    // radiussize              DIAMETER
    // epsilonrule             HHG
    // vdw     45  3.710000000  0.105000000
    // vdw     46  2.700000000  0.020000000  0.910

    for (unsigned int ii = 0; ii < 2; ii++) {
        system.addParticle(1.4007000e+01);
        amoebaWcaDispersionForce->addParticle(3.71e-01 / 2.0e0, 0.105e0 * 4.184e0);

        system.addParticle(1.0080000e+00);
        amoebaWcaDispersionForce->addParticle(2.7e-01 / 2.0e0, 0.02e0 * 4.184e0);

        system.addParticle(1.0080000e+00);
        amoebaWcaDispersionForce->addParticle(2.7e-01 / 2.0e0, 0.02e0 * 4.184e0);

        system.addParticle(1.0080000e+00);
        amoebaWcaDispersionForce->addParticle(2.7e-01 / 2.0e0, 0.02e0 * 4.184e0);
    }

    std::vector<Vec3> positions(numberOfParticles);

    positions[0] = Vec3(1.57493055e-01, 0.00085545e-01, 0.12707740e-01);
    positions[1] = Vec3(1.94285651e-01, -0.81095826e-01, 0.60622179e-01);
    positions[2] = Vec3(1.95030610e-01, 0.80808023e-01, 0.60809749e-01);
    positions[3] = Vec3(1.99309092e-01, -0.00015413e-01, -0.79434283e-01);
    positions[4] = Vec3(-1.65212994e-01, 0.02174061e-01, -0.04650966e-01);
    positions[5] = Vec3(-1.97210029e-01, 0.82501943e-01, 0.47898953e-01);
    positions[6] = Vec3(-0.64164714e-01, 0.01492667e-01, 0.03495515e-01);
    positions[7] = Vec3(-1.98070588e-01, -0.79434787e-01, 0.45336787e-01);

    system.addForce(amoebaWcaDispersionForce);
    ASSERT(!amoebaWcaDispersionForce->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());

    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, platform);

    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> &forces = state.getForces();
    double energy = state.getPotentialEnergy();

    // Force Field X / Tinker expected energy -5.55412658 * 4.184 and forces.
    double expectedEnergy = -5.55412658 * 4.184;
    std::vector<Vec3> expectedForces(numberOfParticles);

    expectedForces[0] = Vec3(-0.33370779, 0.00478895, -0.18780836);
    expectedForces[1] = Vec3(1.15078615, 0.07363123, -0.04881103);
    expectedForces[2] = Vec3(1.14509754, -0.08727205, -0.05086123);
    expectedForces[3] = Vec3(0.97122614, -0.00473535, 0.09417216);
    expectedForces[4] = Vec3(-1.64874538, 0.01074006, -0.08477238);
    expectedForces[5] = Vec3(-1.28585908, 0.00494013, 0.05864744);
    expectedForces[6] = Vec3(1.28677671, -0.00915430, 0.16081161);
    expectedForces[7] = Vec3(-1.28557430, 0.00706132, 0.05862179);

    double tolerance = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);

    // Try changing the particle parameters and make sure it's still correct.
    for (int i = 0; i < numberOfParticles; i++) {
        double radius, epsilon;
        amoebaWcaDispersionForce->getParticleParameters(i, radius, epsilon);
        amoebaWcaDispersionForce->setParticleParameters(i, 0.9 * radius, 2.0 * epsilon);
    }
    LangevinIntegrator integrator2(0.0, 0.1, 0.01);
    Context context2(system, integrator2, platform);
    context2.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy);
    State state2 = context2.getState(State::Forces | State::Energy);
    bool exceptionThrown = false;
    try {
        // This should throw an exception.
        compareForcesEnergy(testName, state1.getPotentialEnergy(), state2.getPotentialEnergy(), state1.getForces(),
                            state2.getForces(), tolerance);
    }
    catch (std::exception ex) {
        exceptionThrown = true;
    }
    ASSERT(exceptionThrown);
    amoebaWcaDispersionForce->updateParametersInContext(context);
    state1 = context.getState(State::Forces | State::Energy);
    compareForcesEnergy(testName, state1.getPotentialEnergy(), state2.getPotentialEnergy(), state1.getForces(),
                        state2.getForces(), tolerance);
}

void setupKernels(int argc, char* argv[]);
void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        setupKernels(argc, argv);
        testWcaDispersionAmmonia();
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
