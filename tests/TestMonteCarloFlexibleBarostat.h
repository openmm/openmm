/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2021 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Lee-Ping Wang                                      *
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
#include "openmm/HarmonicBondForce.h"
#include "openmm/MonteCarloFlexibleBarostat.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/CustomIntegrator.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include "SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testIdealGas() {
    const int numParticles = 64;
    const int frequency = 10;
    const int steps = 1000;
    const double pressure = 3.0;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double temp[] = {300.0, 600.0, 1000.0};
    const double initialVolume = numParticles*BOLTZ*temp[1]/pressureInMD;
    const double initialLength = std::pow(initialVolume, 1.0/3.0);

    // Create a gas of noninteracting particles.

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(initialLength, 0, 0), Vec3(0, 0.5*initialLength, 0), Vec3(0, 0, 2*initialLength));
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(1.0);
        positions[i] = Vec3(initialLength*genrand_real2(sfmt), 0.5*initialLength*genrand_real2(sfmt), 2*initialLength*genrand_real2(sfmt));
    }
    MonteCarloFlexibleBarostat* barostat = new MonteCarloFlexibleBarostat(pressure, temp[0], frequency);
    system.addForce(barostat);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->setUsesPeriodicBoundaryConditions(true);
    system.addForce(bonds); // So it won't complain the system is non-periodic.

    // Test it for three different temperatures.

    for (int i = 0; i < 3; i++) {
        barostat->setDefaultTemperature(temp[i]);
        LangevinIntegrator integrator(temp[i], 0.1, 0.01);
        Context context(system, integrator, platform);
        context.setPositions(positions);

        // Let it equilibrate.

        integrator.step(5000);

        // Now run it for a while and see if the volume is correct.

        double volume = 0.0;
        vector<double> avgPressure(6, 0.0);
        for (int j = 0; j < steps; ++j) {
            Vec3 box[3];
            context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
            volume += box[0][0]*box[1][1]*box[2][2];
            integrator.step(frequency);
            vector<double> pressure;
            barostat->computeCurrentPressure(context, pressure);
            for (int j = 0; j < 6; j++)
                avgPressure[j] += pressure[j];
        }
        volume /= steps;
        double expected = (numParticles+1)*BOLTZ*temp[i]/pressureInMD;
        ASSERT_USUALLY_EQUAL_TOL(expected, volume, 3/std::sqrt((double) steps));
        for (int i = 0; i < 3; i++) {
            avgPressure[i] /= steps;
            ASSERT_USUALLY_EQUAL_TOL(pressure, avgPressure[i], 0.2);
        }
        for (int i = 3; i < 6; i++) {
            avgPressure[i] /= steps;
            ASSERT_USUALLY_EQUAL_TOL(0.0, avgPressure[i], 0.4);
        }
    }
}

void testMoleculeScaling(bool rigid) {
    int numMolecules = 10;
    double initialWidth = 3.0;

    // Create a system of diatomic molecules.
    
    System system;
    Vec3 initialBox[] = {Vec3(initialWidth, 0, 0), Vec3(0, initialWidth, 0), Vec3(0, 0, initialWidth)};
    system.setDefaultPeriodicBoxVectors(initialBox[0], initialBox[1], initialBox[2]);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    system.addForce(bonds);
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    system.addForce(nonbonded);
    MonteCarloFlexibleBarostat* barostat = new MonteCarloFlexibleBarostat(1.0, 300.0, 1, rigid);
    system.addForce(barostat);
    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numMolecules; i++) {
        system.addParticle(1.0);
        system.addParticle(1.0);
        bonds->addBond(2*i, 2*i+1, 0.2, 1000.0);
        nonbonded->addParticle(0.0, 0.1, 1.0);
        nonbonded->addParticle(0.0, 0.1, 1.0);
        Vec3 pos1(initialWidth*genrand_real2(sfmt), initialWidth*genrand_real2(sfmt), initialWidth*genrand_real2(sfmt));
        Vec3 delta(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
        delta /= sqrt(delta.dot(delta));
        positions.push_back(pos1);
        positions.push_back(pos1+delta);
    }

    // Use an integrator that applies the barostat but nothing else.
    
    CustomIntegrator integrator(1.0);
    integrator.addUpdateContextState();

    // Let the barostat make some moves.

    Context context(system, integrator, platform);
    context.setPositions(positions);
    integrator.step(100);
    State state = context.getState(State::Positions);

    // All elements of the box vectors should have changed.

    Vec3 finalBox[3];
    state.getPeriodicBoxVectors(finalBox[0], finalBox[1], finalBox[2]);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j <= i; j++)
            ASSERT(finalBox[i][j] != initialBox[i][j]);

    // See if the molecules were scaled correctly.

    Vec3 boxScale(finalBox[0][0]/initialBox[0][0], finalBox[1][1]/initialBox[1][1], finalBox[2][2]/initialBox[2][2]);
    for (int i = 0; i < numMolecules; i++) {
        Vec3 delta1 = positions[2*i+1]-positions[2*i];
        Vec3 delta2 = state.getPositions()[2*i+1]-state.getPositions()[2*i];
        if (rigid) {
            ASSERT_EQUAL_VEC(delta1, delta2, 1e-5);
        }
        else {
            Vec3 expected(delta1[0]*boxScale[0], delta1[1]*boxScale[1], delta1[2]*boxScale[2]);
            ASSERT_EQUAL_VEC(expected, delta2, 1e-5);
        }
    }
}

void testMolecularGas(bool rigid) {
    const int numMolecules = 64;
    const int frequency = 5;
    const int steps = 5000;
    const double pressure = 3.0;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double temp =300.0;
    const double initialVolume = numMolecules*BOLTZ*temp/pressureInMD;
    const double initialLength = std::pow(initialVolume, 1.0/3.0);

    // Create a gas of noninteracting molecules.

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(initialLength, 0, 0), Vec3(0, 0.5*initialLength, 0), Vec3(0, 0, 2*initialLength));
    MonteCarloFlexibleBarostat* barostat = new MonteCarloFlexibleBarostat(pressure, temp, frequency, rigid);
    system.addForce(barostat);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->setUsesPeriodicBoundaryConditions(true);
    system.addForce(bonds);
    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numMolecules; ++i) {
        system.addParticle(1.0);
        system.addParticle(1.0);
        system.addParticle(1.0);
        Vec3 pos(initialLength*genrand_real2(sfmt), 0.5*initialLength*genrand_real2(sfmt), 2*initialLength*genrand_real2(sfmt));
        bonds->addBond(positions.size(), positions.size()+1, 0.1, 0.0);
        bonds->addBond(positions.size(), positions.size()+2, 0.1, 0.0);
        positions.push_back(pos);
        positions.push_back(pos+Vec3(0.1, 0.0, 0.0));
        positions.push_back(pos+Vec3(0.0, 0.1, 0.0));
    }

    // Simulate it and see if the pressure is correct.

    LangevinIntegrator integrator(temp, 0.1, 0.005);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temp);
    integrator.step(5000);
    vector<double> avgPressure(6, 0.0);
    for (int j = 0; j < steps; ++j) {
        integrator.step(frequency);
        vector<double> pressure;
        barostat->computeCurrentPressure(context, pressure);
        for (int j = 0; j < 6; j++)
            avgPressure[j] += pressure[j];
    }
    for (int i = 0; i < 3; i++) {
        avgPressure[i] /= steps;
        ASSERT_USUALLY_EQUAL_TOL(pressure, avgPressure[i], 0.2);
    }
    for (int i = 3; i < 6; i++) {
        avgPressure[i] /= steps;
        ASSERT_USUALLY_EQUAL_TOL(0.0, avgPressure[i], 0.4);
    }
}

void testRandomSeed() {
    const int numParticles = 8;
    const double temp = 100.0;
    const double pressure = 1.5;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(8, 0, 0), Vec3(0, 8, 0), Vec3(0, 0, 8));
    VerletIntegrator integrator(0.01);
    NonbondedForce* forceField = new NonbondedForce();
    forceField->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(2.0);
        forceField->addParticle((i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(forceField);
    MonteCarloFlexibleBarostat* barostat = new MonteCarloFlexibleBarostat(pressure, temp, 1);
    system.addForce(barostat);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
        velocities[i] = Vec3(0, 0, 0);
    }

    // Try twice with the same random seed.

    barostat->setRandomNumberSeed(5);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(10);
    State state1 = context.getState(State::Positions);
    context.reinitialize();
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(10);
    State state2 = context.getState(State::Positions);

    // Try twice with a different random seed.

    barostat->setRandomNumberSeed(10);
    context.reinitialize();
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(10);
    State state3 = context.getState(State::Positions);
    context.reinitialize();
    context.setPositions(positions);
    context.setVelocities(velocities);
    integrator.step(10);
    State state4 = context.getState(State::Positions);

    // Compare the results.

    for (int i = 0; i < numParticles; i++) {
        for (int j = 0; j < 3; j++) {
            ASSERT(state1.getPositions()[i][j] == state2.getPositions()[i][j]);
            ASSERT(state3.getPositions()[i][j] == state4.getPositions()[i][j]);
            ASSERT(state1.getPositions()[i][j] != state3.getPositions()[i][j]);
        }
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testIdealGas();
        testMoleculeScaling(true);
        testMoleculeScaling(false);
        testMolecularGas(true);
        testMolecularGas(false);
        testRandomSeed();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

