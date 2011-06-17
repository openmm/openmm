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

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/RPMDIntegrator.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testIntegration() {
    const int numParticles = 5;
    const int numCopies = 50;
    const double temperature = 300.0;
    System system;
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(i+1);
        positions[i] = Vec3(i, 0, 0);
    }
    HarmonicBondForce* bonds = new HarmonicBondForce();
    system.addForce(bonds);
    for (int i = 0; i < numParticles-1; i++)
        bonds->addBond(i, i+1, 1.0, 100.0);
    RPMDIntegrator integ(numCopies, temperature, 1.0, 0.001);
    Context context(system, integ);
    context.setPositions(positions);
    const int numSteps = 5000;
    integ.step(1000);
    vector<double> ke(numCopies, 0.0);
    vector<double> rg(numParticles, 0.0);
    const RealOpenMM hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    const double wn = numCopies*BOLTZ*temperature/hbar;
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        vector<State> state(numCopies);
        for (int j = 0; j < numCopies; j++)
            state[j] = integ.getState(j, State::Positions | State::Velocities | State::Energy);
        for (int j = 0; j < numCopies; j++)
            ke[j] += state[j].getKineticEnergy();
        double totalEnergy = 0.0;
        for (int j = 0; j < numCopies; j++) {
            totalEnergy += state[j].getKineticEnergy()+state[j].getPotentialEnergy();
            for (int k = 0; k < numParticles; k++) {
                Vec3 delta = state[j].getPositions()[k]-state[j].getPositions()[j == 0 ? numCopies-1 : j-1];
                totalEnergy += 0.5*system.getParticleMass(k)*wn*wn*delta.dot(delta);
            }
        }
        for (int j = 0; j < numParticles; j++) {
            double rg2 = 0.0;
            for (int k = 0; k < numCopies; k++)
                for (int m = 0; m < numCopies; m++) {
                    Vec3 delta = state[k].getPositions()[j]-state[m].getPositions()[j];
                    rg2 += delta.dot(delta);
                }
            rg[j] += sqrt(rg2/(2*numCopies*numCopies));
        }
    }
    for (int i = 0; i < numCopies; i++) {
        double value = ke[i]/numSteps;
        double expected = 0.5*numCopies*numParticles*3*BOLTZ*temperature;
        printf("%d: %g %g %g\n", i, value, expected, value/expected);
    }
    for (int i = 0; i < numParticles; i++) {
        double value = rg[i]/numSteps;
        double expected = hbar/(2*sqrt(system.getParticleMass(i)*BOLTZ*temperature));
        printf("%d: %g %g %g\n", i, value, expected, value/expected);
    }
}

int main() {
    try {
        Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
        testIntegration();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
