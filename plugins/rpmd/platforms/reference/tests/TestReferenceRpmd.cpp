/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2014 Stanford University and the Authors.      *
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
#include "openmm/CMMotionRemover.h"
#include "openmm/Context.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
#include "openmm/VirtualSite.h"
#include "SimTKOpenMMUtilities.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerRpmdReferenceKernelFactories();

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
    const double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        vector<State> state(numCopies);
        for (int j = 0; j < numCopies; j++)
            state[j] = integ.getState(j, State::Positions | State::Velocities, true);
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

Vec3 calcCM(const vector<Vec3>& values, System& system) {
    Vec3 cm;
    for (int j = 0; j < system.getNumParticles(); ++j) {
        cm[0] += values[j][0]*system.getParticleMass(j);
        cm[1] += values[j][1]*system.getParticleMass(j);
        cm[2] += values[j][2]*system.getParticleMass(j);
    }
    return cm;
}

void testCMMotionRemoval() {
    const int numParticles = 100;
    const int numCopies = 30;
    const double temperature = 300.0;
    const double mass = 1.0;
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(mass);
    system.addForce(new CMMotionRemover());
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
        Vec3 cmPos = calcCM(positions, system);
        for (int j = 0; j < numParticles; j++)
            positions[j] -= cmPos*(1/(mass*numParticles));
        integ.setPositions(i, positions);
    }
    
    // Make sure the CMMotionRemover is getting applied.
    
    for (int i = 0; i < 200; ++i) {
        integ.step(1);
        Vec3 pos;
        for (int j = 0; j < numCopies; j++) {
            State state = integ.getState(0, State::Positions | State::Velocities);
            pos += calcCM(state.getPositions(), system);
        }
        pos *= 1.0/numCopies;
        ASSERT_EQUAL_VEC(Vec3(0,0,0), pos, 0.5);
    }
}

void testVirtualSites() {
    const int gridSize = 3;
    const int numMolecules = gridSize*gridSize*gridSize;
    const int numParticles = numMolecules*3;
    const int numCopies = 10;
    const double spacing = 2.0;
    const double cutoff = 3.0;
    const double boxSize = spacing*(gridSize+1);
    const double temperature = 300.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    HarmonicBondForce* bonds = new HarmonicBondForce();
    system.addForce(bonds);
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    system.addForce(nonbonded);

    // Create a cloud of molecules.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numMolecules; i++) {
        system.addParticle(1.0);
        system.addParticle(1.0);
        system.addParticle(0.0);
        nonbonded->addParticle(-0.2, 0.2, 0.2);
        nonbonded->addParticle(0.1, 0.2, 0.2);
        nonbonded->addParticle(0.1, 0.2, 0.2);
        nonbonded->addException(3*i, 3*i+1, 0, 1, 0);
        nonbonded->addException(3*i, 3*i+2, 0, 1, 0);
        nonbonded->addException(3*i+1, 3*i+2, 0, 1, 0);
        bonds->addBond(3*i, 3*i+1, 1.0, 10000.0);
        system.setVirtualSite(3*i+2, new TwoParticleAverageSite(3*i, 3*i+1, 0.5, 0.5));
    }
    RPMDIntegrator integ(numCopies, temperature, 10.0, 0.001);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    for (int copy = 0; copy < numCopies; copy++) {
        for (int i = 0; i < gridSize; i++)
            for (int j = 0; j < gridSize; j++)
                for (int k = 0; k < gridSize; k++) {
                    Vec3 pos = Vec3(spacing*(i+0.02*genrand_real2(sfmt)), spacing*(j+0.02*genrand_real2(sfmt)), spacing*(k+0.02*genrand_real2(sfmt)));
                    int index = k+gridSize*(j+gridSize*i);
                    positions[3*index] = pos;
                    positions[3*index+1] = Vec3(pos[0]+1.0, pos[1], pos[2]);
                    positions[3*index+2] = Vec3();
                }
        integ.setPositions(copy, positions);
    }

    // Check the temperature and virtual site locations.
    
    const int numSteps = 1000;
    integ.step(1000);
    vector<double> ke(numCopies, 0.0);
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        vector<State> state(numCopies);
        for (int j = 0; j < numCopies; j++) {
            state[j] = integ.getState(j, State::Positions | State::Velocities | State::Forces, true);
            const vector<Vec3>& pos = state[j].getPositions();
            for (int k = 0; k < numMolecules; k++)
                ASSERT_EQUAL_VEC((pos[3*k]+pos[3*k+1])*0.5, pos[3*k+2], 1e-5);
        }
        for (int j = 0; j < numParticles; j++) {
            for (int k = 0; k < numCopies; k++) {
                Vec3 v = state[k].getVelocities()[j];
                ke[k] += 0.5*system.getParticleMass(j)*v.dot(v);
            }
        }
    }
    double meanKE = 0.0;
    for (int i = 0; i < numCopies; i++)
        meanKE += ke[i];
    meanKE /= numSteps*numCopies;
    double expectedKE = 0.5*numCopies*(2*numMolecules)*3*BOLTZ*temperature;
    ASSERT_USUALLY_EQUAL_TOL(expectedKE, meanKE, 1e-2);
}

void testContractions() {
    const int gridSize = 3;
    const int numMolecules = gridSize*gridSize*gridSize;
    const int numParticles = numMolecules*2;
    const int numCopies = 10;
    const double spacing = 2.0;
    const double cutoff = 3.0;
    const double boxSize = spacing*(gridSize+1);
    const double temperature = 300.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    HarmonicBondForce* bonds = new HarmonicBondForce();
    system.addForce(bonds);
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setForceGroup(1);
    nonbonded->setReciprocalSpaceForceGroup(2);
    system.addForce(nonbonded);

    // Create a cloud of molecules.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numMolecules; i++) {
        system.addParticle(1.0);
        system.addParticle(1.0);
        nonbonded->addParticle(-0.2, 0.2, 0.2);
        nonbonded->addParticle(0.2, 0.2, 0.2);
        nonbonded->addException(2*i, 2*i+1, 0, 1, 0);
        bonds->addBond(2*i, 2*i+1, 1.0, 10000.0);
    }
    map<int, int> contractions;
    contractions[1] = 3;
    contractions[2] = 1;
    RPMDIntegrator integ(numCopies, temperature, 50.0, 0.001, contractions);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    for (int copy = 0; copy < numCopies; copy++) {
        for (int i = 0; i < gridSize; i++)
            for (int j = 0; j < gridSize; j++)
                for (int k = 0; k < gridSize; k++) {
                    Vec3 pos = Vec3(spacing*(i+0.02*genrand_real2(sfmt)), spacing*(j+0.02*genrand_real2(sfmt)), spacing*(k+0.02*genrand_real2(sfmt)));
                    int index = k+gridSize*(j+gridSize*i);
                    positions[2*index] = pos;
                    positions[2*index+1] = Vec3(pos[0]+1.0, pos[1], pos[2]);
                }
        integ.setPositions(copy, positions);
    }

    // Check the temperature.
    
    const int numSteps = 1000;
    integ.step(1000);
    vector<double> ke(numCopies, 0.0);
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        vector<State> state(numCopies);
        for (int j = 0; j < numCopies; j++)
            state[j] = integ.getState(j, State::Velocities, true);
        for (int j = 0; j < numParticles; j++) {
            for (int k = 0; k < numCopies; k++) {
                Vec3 v = state[k].getVelocities()[j];
                ke[k] += 0.5*system.getParticleMass(j)*v.dot(v);
            }
        }
    }
    double meanKE = 0.0;
    for (int i = 0; i < numCopies; i++)
        meanKE += ke[i];
    meanKE /= numSteps*numCopies;
    double expectedKE = 0.5*numCopies*numParticles*3*BOLTZ*temperature;
    ASSERT_USUALLY_EQUAL_TOL(expectedKE, meanKE, 1e-2);
}

void testWithoutThermostat() {
    const int numParticles = 20;
    const int numCopies = 10;
    const double temperature = 300.0;
    const double mass = 2.0;
    
    // Create a chain of particles.
    
    System system;
    HarmonicBondForce* bonds = new HarmonicBondForce();
    system.addForce(bonds);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(mass);
        if (i > 0)
            bonds->addBond(i-1, i, 1.0, 1000.0);
    }
    RPMDIntegrator integ(numCopies, temperature, 1.0, 0.001);
    integ.setApplyThermostat(false);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<vector<Vec3> > positions(numCopies);
    for (int i = 0; i < numCopies; i++) {
        positions[i].resize(numParticles);
        for (int j = 0; j < numParticles; j++)
            positions[i][j] = Vec3(0.95*j, 0.01*genrand_real2(sfmt), 0.01*genrand_real2(sfmt));
        integ.setPositions(i, positions[i]);
    }
    
    // Simulate it and see if the energy remains constant.
    
    double initialEnergy;
    int numSteps = 100;
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        double energy = integ.getTotalEnergy();
        if (i == 0)
            initialEnergy = energy;
        else
            ASSERT_EQUAL_TOL(initialEnergy, energy, 1e-4);
    }
}

void testWithBarostat() {
    const int gridSize = 3;
    const int numMolecules = gridSize*gridSize*gridSize;
    const int numParticles = numMolecules*2;
    const int numCopies = 5;
    const double spacing = 2.0;
    const double cutoff = 3.0;
    const double boxSize = spacing*(gridSize+1);
    const double temperature = 300.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    HarmonicBondForce* bonds = new HarmonicBondForce();
    system.addForce(bonds);
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setForceGroup(1);
    nonbonded->setReciprocalSpaceForceGroup(2);
    system.addForce(nonbonded);
    system.addForce(new RPMDMonteCarloBarostat(0.5, 10));

    // Create a cloud of molecules.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numMolecules; i++) {
        system.addParticle(1.0);
        system.addParticle(1.0);
        nonbonded->addParticle(-0.2, 0.2, 0.2);
        nonbonded->addParticle(0.2, 0.2, 0.2);
        nonbonded->addException(2*i, 2*i+1, 0, 1, 0);
        bonds->addBond(2*i, 2*i+1, 1.0, 10000.0);
    }
    RPMDIntegrator integ(numCopies, temperature, 50.0, 0.001);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    for (int copy = 0; copy < numCopies; copy++) {
        for (int i = 0; i < gridSize; i++)
            for (int j = 0; j < gridSize; j++)
                for (int k = 0; k < gridSize; k++) {
                    Vec3 pos = Vec3(spacing*(i+0.02*genrand_real2(sfmt)), spacing*(j+0.02*genrand_real2(sfmt)), spacing*(k+0.02*genrand_real2(sfmt)));
                    int index = k+gridSize*(j+gridSize*i);
                    positions[2*index] = pos;
                    positions[2*index+1] = Vec3(pos[0]+1.0, pos[1], pos[2]);
                }
        integ.setPositions(copy, positions);
    }

    // Check the temperature.
    
    const int numSteps = 500;
    integ.step(100);
    vector<double> ke(numCopies, 0.0);
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        vector<State> state(numCopies);
        for (int j = 0; j < numCopies; j++)
            state[j] = integ.getState(j, State::Velocities, true);
        for (int j = 0; j < numParticles; j++) {
            for (int k = 0; k < numCopies; k++) {
                Vec3 v = state[k].getVelocities()[j];
                ke[k] += 0.5*system.getParticleMass(j)*v.dot(v);
            }
        }
    }
    double meanKE = 0.0;
    for (int i = 0; i < numCopies; i++)
        meanKE += ke[i];
    meanKE /= numSteps*numCopies;
    double expectedKE = 0.5*numCopies*numParticles*3*BOLTZ*temperature;
    ASSERT_USUALLY_EQUAL_TOL(expectedKE, meanKE, 1e-2);
}

int main() {
    try {
        registerRpmdReferenceKernelFactories();
        testFreeParticles();
        testCMMotionRemoval();
        testVirtualSites();
        testContractions();
        testWithoutThermostat();
        testWithBarostat();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
