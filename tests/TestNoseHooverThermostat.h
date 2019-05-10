/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
 * Authors: Andreas Krämer and Andrew C. Simmonett                            *
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
#include "openmm/NoseHooverChain.h"
#include "openmm/VelocityVerletIntegrator.h"
#include "openmm/Context.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace OpenMM;
using namespace std;

void testSingleBond() {
    // Check conservation of system + bath energy for a harmonic oscillator
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    VelocityVerletIntegrator integrator(0.001);
    double temperature = 300; // kelvin
    double collisionFrequency = 0.1; // 1/ps
    int numMTS = 3;
    int numYS = 3;
    int chainLength = 5;
    int numDOF = 6;
    integrator.addNoseHooverChainThermostat(system, temperature, collisionFrequency,
                                            chainLength, numMTS, numYS);
    HarmonicBondForce* forceField = new HarmonicBondForce();
    forceField->addBond(0, 1, 1.5, 1);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);

    integrator.step(10000);

    double mean_temp = 0.0;
    for (int i = 0; i < 10000000; ++i) {
        State state = context.getState(State::Energy);
        double KE = state.getKineticEnergy();
        double PE = state.getPotentialEnergy();
        double time = state.getTime();
        double temperature = 2 * KE / (BOLTZ * numDOF);
        mean_temp = (i*mean_temp + temperature)/(i+1);
        double energy = KE + PE + integrator.computeHeatBathEnergy();
        if(i % 100 == 0)
        std::cout << " Time: " << std::setw(4) << time 
                  << " Temperature: " << std::setw(8) << std::setprecision(3) << temperature
                  << " Mean Temperature: " << std::setw(8) << std::setprecision(3) << mean_temp
                  << " Total Energy: " << std::setw(12) << std::setprecision(8) << energy << std::endl;
        integrator.step(1);
    }
}

void testConstraints() {
    const int numParticles = 8;
    const int numConstraints = 5;
    System system;
    VelocityVerletIntegrator integrator(0.001);
    integrator.setConstraintTolerance(1e-5);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(i%2 == 0 ? 5.0 : 10.0);
        forceField->addParticle((i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
    }
    system.addConstraint(0, 1, 1.0);
    system.addConstraint(1, 2, 1.0);
    system.addConstraint(2, 3, 1.0);
    system.addConstraint(4, 5, 1.0);
    system.addConstraint(6, 7, 1.0);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; ++i) {
        positions[i] = Vec3(i/2, (i+1)/2, 0);
        velocities[i] = Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
    }
    context.setPositions(positions);
    context.setVelocities(velocities);

    // Simulate it and see whether the constraints remain satisfied.

    double initialEnergy = 0.0;
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Positions | State::Energy | State::Velocities | State::Forces);
        for (int j = 0; j < numConstraints; ++j) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(j, particle1, particle2, distance);
            Vec3 p1 = state.getPositions()[particle1];
            Vec3 p2 = state.getPositions()[particle2];
            double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
            ASSERT_EQUAL_TOL(distance, dist, 1e-4);
        }
        double energy = state.getPotentialEnergy()+state.getKineticEnergy();
        if (i == 1)
            initialEnergy = energy;
        else if (i > 1)
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        integrator.step(1);
    }
}
void testNHCPropagation() {
    // test if the velocity scale factor goes to one for a single particle
    // with no forces in the system
    for (int numMTS = 1; numMTS < 5; numMTS++){
        for (int numYS = 1; numYS <=5; numYS+=2){
            for (int chainLength = 2; chainLength < 6; chainLength += 2){
            double temperature = 300; // kelvin
            double collisionFrequency = 0.01; // 1/ps
            VelocityVerletIntegrator integrator(0.001);
            // make system
            System system;
            double mass = 1;
            system.addParticle(mass);
            std::vector<int> particleList {0};
            int chainID = integrator.addNoseHooverChainThermostat(system, temperature, collisionFrequency,
                                                                  chainLength, numMTS, numYS);
            Context context(system, integrator, platform);
            // propagate the chain
            double velocity = 0.1;
            double scale, kineticEnergy;
            for(int i = 0; i < 100; ++i) {
                kineticEnergy = 3 * 0.5 * mass * velocity * velocity;
                scale = integrator.propagateChain(kineticEnergy, chainID);
                velocity *= scale;
            }
            ASSERT_EQUAL_TOL(1, scale,  1e-3);
            }
        }
    }
}

void testPropagateChainConsistentWithPythonReference() {
    VelocityVerletIntegrator integrator(0.001);
    System system;
    double mass = 1;
    system.addParticle(mass);
    double kineticEnergy = 1e6;
    double temperature=300, collisionFrequency=1, chainLength=3, numMTS=3, numYS=3;
    int chainID = integrator.addNoseHooverChainThermostat(system, temperature, collisionFrequency,
                                                                  chainLength, numMTS, numYS);
    Context context(system, integrator, platform);
    double scale = integrator.propagateChain(kineticEnergy, chainID);
    std::cout << std::setw(12) << std::setprecision(10) << scale << std::endl;
    ASSERT_EQUAL_TOL(0.9674732261005896, scale, 1e-5)
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        //testNHCPropagation();
        testPropagateChainConsistentWithPythonReference();
        //testSingleBond();
        //testConstraints();
        //runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
