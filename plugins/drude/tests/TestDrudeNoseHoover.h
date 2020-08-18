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
#include "openmm/CMMotionRemover.h"
#include "openmm/DrudeNoseHooverIntegrator.h"
#include "openmm/Context.h"
#include "openmm/State.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/VirtualSite.h"
#include "openmm/NonbondedForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/System.h"
#include "openmm/DrudeForce.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace OpenMM;
using namespace std;

void build_waterbox(System &system, int gridSize, double polarizability, vector<Vec3> & positions) {
    // Create a box of SWM4-NDP water molecules.  This involves constraints, virtual sites,
    // and Drude particles.
    const int numMolecules = gridSize*gridSize*gridSize;
    const double spacing = 0.8;
    const double boxSize = spacing*(gridSize+1);
    NonbondedForce* nonbonded = new NonbondedForce();
    DrudeForce* drude = new DrudeForce();
    CMMotionRemover* cmm = new CMMotionRemover(1);
    system.addForce(cmm);
    system.addForce(nonbonded);
    system.addForce(drude);
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(1.2);
    nonbonded->setSwitchingDistance(0.8);
    for (int i = 0; i < numMolecules; i++) {
        int startIndex = system.getNumParticles();
        system.addParticle(15.6); // O
        system.addParticle(0.4);  // D
        system.addParticle(1.0);  // H1
        system.addParticle(1.0);  // H2
        system.addParticle(0.0);  // M
        nonbonded->addParticle(1.71636, 0.318395, 0.21094*4.184);
        nonbonded->addParticle(-1.71636, 1, 0);
        nonbonded->addParticle(0.55733, 1, 0);
        nonbonded->addParticle(0.55733, 1, 0);
        nonbonded->addParticle(-1.11466, 1, 0);
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < j; k++)
                nonbonded->addException(startIndex+j, startIndex+k, 0, 1, 0);
        system.addConstraint(startIndex, startIndex+2, 0.09572);
        system.addConstraint(startIndex, startIndex+3, 0.09572);
        system.addConstraint(startIndex+2, startIndex+3, 0.15139);
        system.setVirtualSite(startIndex+4, new ThreeParticleAverageSite(startIndex, startIndex+2, startIndex+3, 0.786646558, 0.106676721, 0.106676721));
        drude->addParticle(startIndex+1, startIndex, -1, -1, -1, -1.71636, polarizability, 1, 1);
    }
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            for (int k = 0; k < gridSize; k++) {
                Vec3 pos(i*spacing, j*spacing, k*spacing);
                positions.push_back(pos);
                positions.push_back(pos);
                positions.push_back(pos+Vec3(0.09572, 0, 0));
                positions.push_back(pos+Vec3(-0.023999, 0.092663, 0));
                positions.push_back(pos);
            }
        }
    }
}

void testWaterBox() {
    // Create a box of SWM4-NDP water molecules.  This involves constraints, virtual sites,
    // and Drude particles.
    System system;
    const int gridSize = 3;
    vector<Vec3> positions;
    double polarizability = ONE_4PI_EPS0*1.71636*1.71636/(100000*4.184);

    build_waterbox(system, gridSize, polarizability, positions);

    const int numMolecules = gridSize*gridSize*gridSize;
    int numStandardDof = 3*3*numMolecules - system.getNumConstraints();
    int numDrudeDof = 3*numMolecules;
    int numDof = numStandardDof+numDrudeDof;
    const double temperature = 300.0;
    const double temperatureDrude = 10.0;

    // Simulate it and check the temperature.
    int chainLength = 4;
    int numMTS = 4;
    int numYS = 5;
    // N.B. These are higher frequencies than recommeded for production runs, but are used
    // here to achieve rapid equilibration to the target temperature, allowing a short run
    double frequency = 1000.0;
    double frequencyDrude = 1000.0;
    int randomSeed = 100;
    DrudeNoseHooverIntegrator integ(temperature, frequency,
                                    temperatureDrude, frequencyDrude, 0.0004,
                                    chainLength, numMTS, numYS);
    Context context(system, integ, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temperature, randomSeed);
    context.applyConstraints(1e-6);

    // Equilibrate
    integ.step(1500);

    double TOL = 1.5;
    try {
        if (platform.getPropertyValue(context, "Precision") != "double") {
            TOL = 2.0;
        }
    } catch(OpenMMException) {
        // The defaults above are for double precision, which is assumed in this case
    }

    // Compute the internal and center of mass temperatures.
    double totalKE = 0;
    const int numSteps = 500;
    double meanTemp = 0.0;
    double meanDrudeTemp = 0.0;
    double meanConserved = 0.0;
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        State state = context.getState(State::Energy);
        double KE = state.getKineticEnergy();
        double PE = state.getPotentialEnergy();
        double fullKE = integ.computeTotalKineticEnergy();
        double drudeKE = integ.computeDrudeKineticEnergy();
        double temp = KE/(0.5*numStandardDof*BOLTZ);
        double drudeTemp = drudeKE/(0.5*numDrudeDof*BOLTZ);
        meanTemp = (i*meanTemp + temp)/(i+1);
        meanDrudeTemp = (i*meanDrudeTemp + drudeTemp)/(i+1);
        double heatBathEnergy = integ.computeHeatBathEnergy();
        double conserved = PE + fullKE + heatBathEnergy;
        meanConserved = (i*meanConserved + conserved)/(i+1);
        totalKE += KE;
        ASSERT(fabs(meanConserved - conserved) < TOL);
    }
    totalKE /= numSteps;
    ASSERT_USUALLY_EQUAL_TOL(temperature, meanTemp,  0.03);
    ASSERT_USUALLY_EQUAL_TOL(temperatureDrude, meanDrudeTemp,  0.03);
}


double testWaterBoxWithHardWallConstraint(double hardWallConstraint){
    // Create a box of SWM4-NDP water molecules.  This involves constraints, virtual sites,
    // and Drude particles.
    System system;
    const int gridSize = 3;
    vector<Vec3> positions;

    double polarizability = 1e-2;
    build_waterbox(system, gridSize, polarizability, positions);

    const int numMolecules = gridSize*gridSize*gridSize;
    int numStandardDof = 3*3*numMolecules - system.getNumConstraints();
    int numDrudeDof = 3*numMolecules;
    int numDof = numStandardDof+numDrudeDof;
    const double temperature = 300.0;
    const double temperatureDrude = 10.0;

    // Simulate it and check the temperature.
    int chainLength = 4;
    int numMTS = 3;
    int numYS = 3;
    double frequency = 25.0;
    double frequencyDrude = 25.0;
    int randomSeed = 100;
    DrudeNoseHooverIntegrator integ(temperature, frequency,
                                    temperatureDrude, frequencyDrude, 0.0005,
                                    chainLength, numMTS, numYS);
    integ.setMaxDrudeDistance(hardWallConstraint);
    Context context(system, integ, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temperature, randomSeed);
    context.applyConstraints(1e-6);

    // Equilibrate.
    integ.step(50);

    // Compute the internal and center of mass temperatures.
    double totalKE = 0;
    const int numSteps = 500;
    double meanTemp = 0.0;
    double meanDrudeTemp = 0.0;
    double meanConserved = 0.0;
    double maxR = 0.0;
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        State state = context.getState(State::Energy | State::Positions);
        const auto & positions = state.getPositions();
        for(int mol = 0; mol < gridSize*gridSize*gridSize; ++mol) {
            auto dR = positions[5*mol+1] - positions[5*mol];
            maxR = std::max(maxR, std::sqrt(dR.dot(dR)));
        }
    }
    return maxR;
}

void testInitialTemperature() {
    // Check temperature initialization for a collection of randomly placed particles
    const int numRealParticles = 50000;
    const int numParticles = 2 * numRealParticles;
    const int nDoF = 3 * numRealParticles;
    const double targetTemperature = 300;
    const double drudeTemperature = 1;
    const double realMass = 10;
    const double drudeMass = 1;
    System system;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    std::vector<Vec3> positions(numParticles);
    DrudeForce* drude = new DrudeForce();

    for (int i = 0; i < numRealParticles; i++) {
        system.addParticle(realMass);
        system.addParticle(drudeMass);
        positions[2*i][0] = genrand_real2(sfmt);
        positions[2*i][1] = genrand_real2(sfmt);
        positions[2*i][2] = genrand_real2(sfmt);
        positions[2*i+1][0] = positions[2*i][0] + 0.01*genrand_real2(sfmt);
        positions[2*i+1][1] = positions[2*i][1] + 0.01*genrand_real2(sfmt);
        positions[2*i+1][2] = positions[2*i][2] + 0.01*genrand_real2(sfmt);
        drude->addParticle(2*i+1, 2*i, -1, -1, -1, -1.0, 0.001, 1, 1);
    }
    system.addForce(drude);
    CMMotionRemover* cmm = new CMMotionRemover(1);
    system.addForce(cmm);

    DrudeNoseHooverIntegrator integrator(targetTemperature, 25, drudeTemperature, 25, 0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(targetTemperature);
    auto velocities = context.getState(State::Velocities).getVelocities();
    double comKineticEnergy = 0;
    double relKineticEnergy = 0;
    for (int i = 0; i < numRealParticles; i++) {
        int m1 = realMass;
        int m2 = drudeMass;
        Vec3 v1 = velocities[2*i];
        Vec3 v2 = velocities[2*i + 1];
        double invMass = 1.0 / (m1 + m2);
        double redMass = m1 * m2 * invMass;
        double fracM1 = m1 * invMass;
        double fracM2 = m2 * invMass;
        Vec3 comVelocity = fracM1 * v1 + fracM2 * v2;
        Vec3 relVelocity = v2 - v1;

        comKineticEnergy += 0.5 * (m1 + m2) * comVelocity.dot(comVelocity);
        relKineticEnergy += 0.5 * redMass * relVelocity.dot(relVelocity);
    }
    double comTemperature = (2*comKineticEnergy / (nDoF*BOLTZ));
    double relTemperature = (2*relKineticEnergy / (nDoF*BOLTZ));
    ASSERT_USUALLY_EQUAL_TOL(targetTemperature, comTemperature, 0.01);
    ASSERT_USUALLY_EQUAL_TOL(drudeTemperature, relTemperature, 0.01);
}

void setupKernels(int argc, char* argv[]);
void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        setupKernels(argc, argv);

        testWaterBox();
        double maxDrudeDistance = 0.005;
        double observedDrudeDistance = testWaterBoxWithHardWallConstraint(0.0);
        ASSERT(observedDrudeDistance > maxDrudeDistance);
        observedDrudeDistance = testWaterBoxWithHardWallConstraint(maxDrudeDistance);
        // Remove later: just trying to find out why Jenkins is upset
        if(observedDrudeDistance >= maxDrudeDistance) printf("Max distance %16.10f\n", observedDrudeDistance);
        ASSERT(observedDrudeDistance < maxDrudeDistance);
        testInitialTemperature();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

