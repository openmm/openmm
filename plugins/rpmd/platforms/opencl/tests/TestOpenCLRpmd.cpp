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
 * This tests the OpenCL implementation of RPMDIntegrator.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/Context.h"
#include "openmm/CustomNonbondedForce.h"
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

extern "C" OPENMM_EXPORT void registerRPMDOpenCLKernelFactories();

void testFreeParticles() {
    const int numParticles = 100;
    const int numCopies = 30;
    const double temperature = 300.0;
    const double mass = 1.0;
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(mass);
    RPMDIntegrator integ(numCopies, temperature, 10.0, 0.001);
    Platform& platform = Platform::getPlatformByName("OpenCL");
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

void testParaHydrogen() {
    const int numParticles = 32;
    const int numCopies = 12;
    const double temperature = 25.0;
    const double mass = 2.0;
    const double boxSize = 1.1896;
    const int numSteps = 1000;
    const int numBins = 200;
    const double reference[] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 4.932814042206152e-5, 1.244331241336431e-4, 4.052316284060125e-4,
        1.544810863683946e-3, 4.376197806690222e-3, 1.025847561714293e-2, 2.286702037465422e-2,
        4.371052180263602e-2, 7.518538770734748e-2, 0.122351534531647, 0.185758975626622,
        0.266399984652322, 0.363380262153250, 0.473696401293219, 0.595312098494172,
        0.726049519422861, 0.862264551954547, 0.991102029379444, 1.1147503922535,
        1.23587006992066, 1.33495411932817, 1.42208208736987, 1.49273884004107,
        1.54633319690403, 1.58714702233941, 1.60439217751355, 1.61804190608902,
        1.60680198476058, 1.58892222973695, 1.56387607986781, 1.52629494593350,
        1.48421439018970, 1.43656176771959, 1.38752775598872, 1.33310695719931,
        1.28363477223121, 1.23465642750248, 1.18874848666326, 1.14350496170519,
        1.10292486009936, 1.06107270157688, 1.02348927970441, 0.989729345271297,
        0.959273446941802, 0.932264875865758, 0.908818658748942, 0.890946420768315,
        0.869332737718165, 0.856401736350349, 0.842370069917020, 0.834386614237393,
        0.826268072171045, 0.821547250199453, 0.818786865315836, 0.819441757028076,
        0.819156933383128, 0.822275325148621, 0.828919078023881, 0.837233720599450,
        0.846961908186718, 0.855656955481099, 0.864520333201247, 0.876082425547566,
        0.886950044046000, 0.900275658318995
    };
    
    // Create a box of para-hydrogen.
    
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(mass);
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize,0,0), Vec3(0,boxSize,0), Vec3(0,0,boxSize));
    CustomNonbondedForce* nb = new CustomNonbondedForce("2625.49963*(exp(1.713-1.5671*p-0.00993*p*p)-(12.14/p^6+215.2/p^8-143.1/p^9+4813.9/p^10)*(step(rc-p)*exp(-(rc/p-1)^2)+1-step(rc-p))); p=r/0.05291772108; rc=8.32");
    nb->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    nb->setCutoffDistance(boxSize/2);
    vector<double> params;
    for (int i = 0; i < numParticles; i++)
        nb->addParticle(params);
    system.addForce(nb);
    RPMDIntegrator integ(numCopies, temperature, 10.0, 0.0005);
    Platform& platform = Platform::getPlatformByName("OpenCL");
    Context context(system, integ, platform);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++)
        positions[i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
    for (int i = 0; i < numCopies; i++)
        integ.setPositions(i, positions);
    integ.step(1000);
    
    // Simulate it.
    
    vector<int> counts(numBins, 0);
    const double invBoxSize = 1.0/boxSize;
    double meanKE = 0.0;
    const double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    for (int step = 0; step < numSteps; step++) {
        integ.step(20);
        vector<State> states(numCopies);
        for (int i = 0; i < numCopies; i++)
            states[i] = integ.getState(i, State::Positions | State::Forces);
        
        // Record the radial distribution function.
        
        const vector<Vec3>& pos = states[0].getPositions();
        for (int j = 0; j < numParticles; j++)
            for (int k = 0; k < j; k++) {
                Vec3 delta = pos[j]-pos[k];
                delta[0] -= floor(delta[0]*invBoxSize+0.5)*boxSize;
                delta[1] -= floor(delta[1]*invBoxSize+0.5)*boxSize;
                delta[2] -= floor(delta[2]*invBoxSize+0.5)*boxSize;
                double dist = sqrt(delta.dot(delta));
                int bin = (int) (numBins*(dist/boxSize));
                counts[bin]++;
            }
        
        // Calculate the quantum contribution to the kinetic energy.
        
        vector<Vec3> centroids(numParticles, Vec3());
        for (int i = 0; i < numCopies; i++) {
            const vector<Vec3>& pos = states[i].getPositions();
            for (int j = 0; j < numParticles; j++)
                centroids[j] += pos[j];
        }
        for (int j = 0; j < numParticles; j++)
            centroids[j] *= 1.0/numCopies;
        double ke = 0.0;
        for (int i = 0; i < numCopies; i++) {
            const vector<Vec3>& pos = states[i].getPositions();
            const vector<Vec3>& f = states[i].getForces();
            for (int j = 0; j < numParticles; j++) {
                Vec3 delta = centroids[j]-pos[j];
                ke += delta.dot(f[j]);
            }
        }
        meanKE += ke/(2*numCopies*numParticles);
    }
    
    // Check against expected values.
    
    double scale = (boxSize*boxSize*boxSize)/(numSteps*0.5*numParticles*numParticles);
    for (int i = 0; i < numBins/2; i++) {
        double r1 = i*boxSize/numBins;
        double r2 = (i+1)*boxSize/numBins;
        double volume = (4.0/3.0)*M_PI*(r2*r2*r2-r1*r1*r1);
        ASSERT_USUALLY_EQUAL_TOL(reference[i], scale*counts[i]/volume, 0.1);
    }
    meanKE /= numSteps*BOLTZ;
    ASSERT_USUALLY_EQUAL_TOL(60.0, 1.5*temperature+meanKE, 0.02);
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
    Platform& platform = Platform::getPlatformByName("OpenCL");
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
    Platform& platform = Platform::getPlatformByName("OpenCL");
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
    Platform& platform = Platform::getPlatformByName("OpenCL");
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
    Platform& platform = Platform::getPlatformByName("OpenCL");
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
    Platform& platform = Platform::getPlatformByName("OpenCL");
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

int main(int argc, char* argv[]) {
    try {
        registerRPMDOpenCLKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("OpenCL").setPropertyDefaultValue("Precision", string(argv[1]));
        testFreeParticles();
        testParaHydrogen();
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
