/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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
#include "openmm/CustomBondForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/QTBIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testHarmonic() {
    // Create a collection of uncoupled harmonic oscillators.

    int numParticles = 10;
    double temperature = 300.0;
    double mass = 1.0;
    System system;
    vector<Vec3> positions;
    vector<double> k;
    CustomExternalForce* force = new CustomExternalForce("0.5*k*(x*x+y*y+z*z)");
    system.addForce(force);
    force->addPerParticleParameter("k");
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(mass);
        positions.push_back(Vec3());
        k.push_back(1000.0*(i+1)*(i+1));
        force->addParticle(i, {k[i]});
    }
    QTBIntegrator integrator(temperature, 10.0, 0.001);
    integrator.setCutoffFrequency(500.0);
    integrator.setDefaultAdaptationRate(1e-4);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temperature);

    // Compute the average energy of each particle over a simulation.

    integrator.step(10000);
    vector<double> energy(numParticles, 0.0);
    int numSteps = 100000;
    for (int i = 0; i < numSteps; i++) {
        integrator.step(10);
        State state = context.getState(State::Positions);
        for (int j = 0; j < numParticles; j++) {
            Vec3 p = state.getPositions()[j];
            energy[j] += 0.5*k[j]*p.dot(p);
        }
    }
    for (int i = 0; i < numParticles; i++)
        energy[i] /= numSteps;

    // Compare to the expected distribution.

    for (int i = 0; i < numParticles; i++) {
        double w = sqrt(k[i]/mass);
        double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
        double kT = BOLTZ*temperature;
        double expected = 1.5*hbar*w*(0.5+1/(exp(hbar*w/kT)-1));
        printf("%g %g %g %g\n", w, expected, energy[i], energy[i]/expected);
        ASSERT_USUALLY_EQUAL_TOL(expected, energy[i], 0.05);
    }
}

void testCoupledHarmonic() {
    // Create a collection of weakly coupled harmonic oscillators.

    int numFrequencies = 8;
    int numReplicas = 4;
    int numParticles = numFrequencies*numReplicas;
    double temperature = 10.0;
    double mass = 1.0;
    System system;
    vector<Vec3> positions(numParticles);
    vector<double> k;
    CustomExternalForce* force = new CustomExternalForce("0.5*k*(x*x+y*y+z*z)");
    system.addForce(force);
    force->addPerParticleParameter("k");
    CustomBondForce* bonds = new CustomBondForce("sin(100*r)");
    system.addForce(bonds);
    QTBIntegrator integrator(temperature, 50.0, 0.001);
    integrator.setCutoffFrequency(500.0);
    integrator.setDefaultAdaptationRate(0.5);
    for (int j = 0; j < numReplicas; j++) {
        int base = system.getNumParticles();
        for (int i = 0; i < numFrequencies; i++) {
            system.addParticle(mass);
            k.push_back(1000.0*(i+1)*(i+1));
            force->addParticle(i+base, {k[i]});
            integrator.setParticleType(i+base, i);
        }
        for (int i = 0; i < numFrequencies; i++)
            for (int k = 0; k < i; k++)
                bonds->addBond(i+base, k+base);
    }
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temperature);

    // Compute the average energy of each particle over a simulation.

    integrator.step(200000);
    integrator.setDefaultAdaptationRate(0.05);
    context.reinitialize(true);
    vector<double> energy(numFrequencies, 0.0);
    int numSteps = 50000;
    for (int i = 0; i < numSteps; i++) {
        integrator.step(10);
        State state = context.getState(State::Positions);
        for (int j = 0; j < numParticles; j++) {
            Vec3 p = state.getPositions()[j];
            energy[j%numFrequencies] += 0.5*k[j]*p.dot(p);
        }
    }
    for (int i = 0; i < numFrequencies; i++)
        energy[i] /= numSteps*numReplicas;

    // Compare to the expected distribution.

    for (int i = 0; i < numFrequencies; i++) {
        double w = sqrt(k[i]/mass);
        double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
        double kT = BOLTZ*temperature;
        double expected = 1.5*hbar*w*(0.5+1/(exp(hbar*w/kT)-1));
        printf("%g %g %g %g\n", w, expected, energy[i], energy[i]/expected);
        ASSERT_USUALLY_EQUAL_TOL(expected, energy[i], 0.15);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testHarmonic();
        testCoupledHarmonic();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
