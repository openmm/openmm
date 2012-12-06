/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
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
 * This tests the Reference implementation of RPMDIntegrator.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/RPMDIntegrator.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testFreeParticles() {
    const int numParticles = 100;
    const int numCopies = 30;
    const double temperature = 300.0;
    const double mass = 1.0;
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(mass);
    RPMDIntegrator integ(numCopies, temperature, 10.0, 0.001);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numCopies; i++)
    {
        for (int j = 0; j < numParticles; j++)
            positions[j] = Vec3(0.02*genrand_real2(sfmt), 0.02*genrand_real2(sfmt), 0.02*genrand_real2(sfmt));
        integ.setPositions(i, positions);
    }
    const int numSteps = 1000;
    integ.step(1000);
    vector<double> ke(numCopies, 0.0);
    vector<double> rg(numParticles, 0.0);
    const RealOpenMM hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        vector<State> state(numCopies);
        for (int j = 0; j < numCopies; j++)
            state[j] = integ.getState(j, State::Positions | State::Velocities);
        for (int j = 0; j < numParticles; j++) {
            double rg2 = 0.0;
            for (int k = 0; k < numCopies; k++) {
                Vec3 v = state[k].getVelocities()[j];
                ke[k] += 0.5*mass*v.dot(v);
                for (int m = 0; m < numCopies; m++) {
                    Vec3 delta = state[k].getPositions()[j]-state[m].getPositions()[j];
                    rg2 += delta.dot(delta);
                }
            }
            rg[j] += rg2/(2*numCopies*numCopies);
        }
    }
    double meanKE = 0.0;
    for (int i = 0; i < numCopies; i++)
        meanKE += ke[i];
    meanKE /= numSteps*numCopies;
    double expectedKE = 0.5*numCopies*numParticles*3*BOLTZ*temperature;
    ASSERT_USUALLY_EQUAL_TOL(expectedKE, meanKE, 1e-2);
    double meanRg2 = 0.0;
    for (int i = 0; i < numParticles; i++)
        meanRg2 += rg[i];
    meanRg2 /= numSteps*numParticles;
    double expectedRg = hbar/(2*sqrt(mass*BOLTZ*temperature));
    ASSERT_USUALLY_EQUAL_TOL(expectedRg, sqrt(meanRg2), 1e-3);
}

int main() {
    try {
        testFreeParticles();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
