/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

/**
 * This tests the Reference implementation of DrudeLangevinIntegrator.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/DrudeForce.h"
#include "openmm/DrudeLangevinIntegrator.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testSinglePair() {
    const double temperature = 300.0;
    const double temperatureDrude = 10.0;
    const double k = 1.5;
    const double charge = 0.1;
    const double alpha = charge*charge/k;
    const double mass1 = 1.0;
    const double mass2 = 0.1;
    const double totalMass = mass1+mass2;
    const double reducedMass = (mass1*mass2)/(mass1+mass2);
    System system;
    system.addParticle(mass1);
    system.addParticle(mass2);
    DrudeForce* drude = new DrudeForce();
    drude->addParticle(1, 0, -1, -1, -1, charge, alpha, 1, 1);
    system.addForce(drude);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);
    DrudeLangevinIntegrator integ(temperature, 10.0, temperatureDrude, 10.0, 0.003);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    context.setPositions(positions);
    
    // Equilibrate.
    
    integ.step(1000);
    
    // Compute the internal and center of mass temperatures.
    
    double keCM = 0, keInternal = 0;
    int numSteps = 10000;
    for (int i = 0; i < numSteps; i++) {
        integ.step(10);
        State state = context.getState(State::Velocities);
        const vector<Vec3>& vel = state.getVelocities();
        Vec3 velCM = vel[0]*(mass1/totalMass) + vel[1]*(mass2/totalMass);
        keCM += 0.5*totalMass*velCM.dot(velCM);
        Vec3 velInternal = vel[0]-vel[1];
        keInternal += 0.5*reducedMass*velInternal.dot(velInternal);
    }
    ASSERT_USUALLY_EQUAL_TOL(3*0.5*BOLTZ*temperature, keCM/numSteps, 0.1);
    ASSERT_USUALLY_EQUAL_TOL(3*0.5*BOLTZ*temperatureDrude, keInternal/numSteps, 0.01);
}

int main() {
    try {
        testSinglePair();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
