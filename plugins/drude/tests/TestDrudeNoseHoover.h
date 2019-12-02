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
#include <unistd.h>

using namespace OpenMM;
using namespace std;
extern "C" OPENMM_EXPORT void registerDrudeReferenceKernelFactories();
extern "C" OPENMM_EXPORT void registerKernelFactories();
Platform& initializePlatform(int argc, char* argv[]);

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

void testWaterBox(Platform& platform) {
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
    double frequency = 800.0;
    double frequencyDrude = 2000.0;
    int randomSeed = 100;
    DrudeNoseHooverIntegrator integ(temperature, frequency, 
                                    temperatureDrude, frequencyDrude, 0.0005, 
                                    chainLength, numMTS, numYS);
    Context context(system, integ, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temperature, randomSeed);
    std::vector<Vec3> velocities = context.getState(State::Velocities).getVelocities();
    for (int i = 0; i < numMolecules; i++){
        Vec3 noize;
        for (int j = 0; j < 3; j++){
            noize[j] = float(((i+18311)*(j+18253) * 313419097822414) % 18313) / float(18313);
            noize[j] *= sqrt(3 * BOLTZ * temperatureDrude / 0.4);
        }
        velocities[5*i+1] = velocities[5*i] + noize;
    }
    context.setVelocities(velocities);
    context.applyConstraints(1e-6);

    // Equilibrate.
    integ.step(500);

    // Compute the internal and center of mass temperatures.

    double totalKE = 0;
    const int numSteps = 4000;
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
#define DEBUG 0
#if DEBUG
        if(i%10 == 0)
        std::cout << std::setw(6) << i
                  << std::setprecision(8) << std::setw(16) << KE
                  << std::setprecision(8) << std::setw(16) << drudeKE
                  << std::setprecision(8) << std::setw(16) << meanTemp
                  << std::setprecision(8) << std::setw(16) << meanDrudeTemp
                  << std::setprecision(8) << std::setw(16) << heatBathEnergy
                  << std::setprecision(8) << std::setw(16) << fullKE
                  << std::setprecision(8) << std::setw(16) << conserved
                  << std::endl;
#endif
        totalKE += KE;
        ASSERT(fabs(meanConserved - conserved) < 0.6);
    }
    totalKE /= numSteps;
    ASSERT_USUALLY_EQUAL_TOL(temperature, meanTemp,  0.004);
    ASSERT_USUALLY_EQUAL_TOL(temperatureDrude, meanDrudeTemp,  0.004);
}


double testWaterBoxWithHardWallConstraint(Platform& platform, double hardWallConstraint){
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
    std::vector<Vec3> velocities = context.getState(State::Velocities).getVelocities();
    for (int i = 0; i < numMolecules; i++){
        Vec3 noize;
        for (int j = 0; j < 3; j++){
            noize[j] = float(((i+18311)*(j+18253) * 313419097822414) % 18313) / float(18313);
            noize[j] *= sqrt(3 * BOLTZ * temperatureDrude / 0.4);
        }
        velocities[5*i+1] = velocities[5*i] + noize;
    }
    context.setVelocities(velocities);
    context.applyConstraints(1e-6);
    // Equilibrate.
    
    integ.step(10);
    // Compute the internal and center of mass temperatures.
    
    double totalKE = 0;
    const int numSteps = 10;
    double meanTemp = 0.0;
    double meanDrudeTemp = 0.0;
    double meanConserved = 0.0;
    double maxR = 0.0;
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        State state = context.getState(State::Energy | State::Positions);
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
        const auto & positions = state.getPositions();
        for(int mol = 0; mol < gridSize*gridSize*gridSize; ++mol) {
            auto dR = positions[5*mol+1] - positions[5*mol];
            maxR = std::max(maxR, std::sqrt(dR.dot(dR)));
        }
#if DEBUG
        if(i%1 == 0)
        std::cout << std::setw(6) << i
                  << std::setprecision(8) << std::setw(16) << KE
                  << std::setprecision(8) << std::setw(16) << drudeKE
                  << std::setprecision(8) << std::setw(16) << meanTemp
                  << std::setprecision(8) << std::setw(16) << meanDrudeTemp
                  << std::setprecision(8) << std::setw(16) << heatBathEnergy
                  << std::setprecision(8) << std::setw(16) << fullKE
                  << std::setprecision(8) << std::setw(16) << conserved
                  << std::setprecision(8) << std::setw(16) << maxR
                  << std::endl;
#endif
        totalKE += KE;
    }
    totalKE /= numSteps;
    return maxR;
}

double testPositionsAfterShortRun(Platform &platform) {   
    // Make sure the different platforms are consistent after a short run

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
    Context context(system, integ, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temperature, randomSeed);
    std::vector<Vec3> velocities = context.getState(State::Velocities).getVelocities();
    for (int i = 0; i < numMolecules; i++){
        Vec3 noize;
        for (int j = 0; j < 3; j++){
            noize[j] = float(((i+18311)*(j+18253) * 313419097822414) % 18313) / float(18313);
            noize[j] *= sqrt(3 * BOLTZ * temperatureDrude / 0.4);
        }
        velocities[5*i+1] = velocities[5*i] + noize;
    }
    context.setVelocities(velocities);
    context.applyConstraints(1e-6);
    // Equilibrate.
    integ.step(10);
    const auto &state = context.getState(State::Velocities | State::Positions);
    const auto &v = state.getVelocities();
    const auto &p = state.getPositions();

    const double TOL = 5e-4;
    ASSERT_EQUAL_VEC(Vec3(      0.35675068,       0.65083786,       0.27473552), v[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(     0.099853537,        0.2840745,      0.069048963), v[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.24183664,       0.12975303,       -1.1057667), v[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.99646242,       0.73498627,        1.3332999), v[3], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.22331845,     -0.071958999,       0.34057062), v[5], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.52415795,      -0.39856875,       0.34114842), v[6], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.80329166,       0.91595337,        -2.784552), v[7], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -1.3914175,      -0.45328548,       0.18652167), v[8], TOL);
    ASSERT_EQUAL_VEC(Vec3(       0.4929438,      0.043788117,      0.096765579), v[10], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.20558536,      -0.20498271,       0.38598319), v[11], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.42455333,      -0.83355484,       0.83189942), v[12], TOL);
    ASSERT_EQUAL_VEC(Vec3(       1.2999188,       0.19737524,       0.67573237), v[13], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.27898695,       0.33332618,       0.09044357), v[15], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.38330415,      -0.41558649,      -0.08833055), v[16], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.29949837,      0.019242996,       -0.4435692), v[17], TOL);
    ASSERT_EQUAL_VEC(Vec3(     0.064409407,        0.3854447,       0.84016629), v[18], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -1.0318902,      -0.53333331,       0.43815812), v[20], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -1.1861828,       -1.3165392,       0.43687045), v[21], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -1.1182471,       -1.7061111,       0.95895761), v[22], TOL);
    ASSERT_EQUAL_VEC(Vec3(     0.090421672,      -0.32791357,       0.94539742), v[23], TOL);
    ASSERT_EQUAL_VEC(Vec3(   -0.0047757094,       0.65850882,      -0.12948891), v[25], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.21030035,       -0.1260296,     0.0062573642), v[26], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.54136664,       -2.2377762,        1.1346483), v[27], TOL);
    ASSERT_EQUAL_VEC(Vec3(       2.8681838,       0.94407246,       0.16030045), v[28], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.28275107,      -0.26207416,      -0.21379774), v[30], TOL);
    ASSERT_EQUAL_VEC(Vec3(       0.1932711,      -0.50089485,      -0.34898139), v[31], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.17557494,      -0.41978113,        1.2024479), v[32], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.43345702,      -0.22569473,      -0.37391908), v[33], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.31114587,      -0.20169499,      -0.24499386), v[35], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.19272819,      -0.37087802,      -0.22718454), v[36], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.16162891,       0.86008276,        1.0320346), v[37], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.93698344,      -0.77471996,          1.44348), v[38], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.34928463,     -0.041834318,      -0.14345014), v[40], TOL);
    ASSERT_EQUAL_VEC(Vec3(       0.2098324,      -0.25200593,     -0.097569375), v[41], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.28597108,       -1.0748928,      -0.53100758), v[42], TOL);
    ASSERT_EQUAL_VEC(Vec3(       1.2903163,      0.035927006,       -1.6031631), v[43], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.16759809,     -0.013783732,       -0.4579853), v[45], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.46506268,      -0.25388673,      -0.59233902), v[46], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.11648621,     -0.050055908,       0.52736643), v[47], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.10540497,     -0.066070446,       0.35401521), v[48], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.57336726,       -0.5475818,      0.019023218), v[50], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -1.1881536,      -0.75891282,     -0.053482097), v[51], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.60733307,       -1.1807411,       0.53263676), v[52], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.15417722,      -0.67460933,         -2.29721), v[53], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -0.4956203,      0.048324709,     -0.094323454), v[55], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -1.1523683,       -0.1901027,       0.14323752), v[56], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.58710395,       0.35404589,       -1.3777924), v[57], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.81763406,     -0.041441416,      -0.12522791), v[58], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -0.2572238,      -0.14896639,       -0.8763284), v[60], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.93036275,      -0.97084853,      -0.98645054), v[61], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.30721487,        0.5514243,      -0.27501022), v[62], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -1.1273817,       -1.2297732,        2.8713348), v[63], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.31531418,      0.021452941,       0.42715776), v[65], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -1.0044077,      -0.75542034,       0.33564906), v[66], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.33976679,       0.11366743,        1.1047837), v[67], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.44083608,     -0.021987305,       0.85722954), v[68], TOL);
    ASSERT_EQUAL_VEC(Vec3(    -0.031147356,       0.33417734,      0.091803041), v[70], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.77436714,      -0.42999344,       0.19097312), v[71], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -0.1035208,      -0.82460701,       0.30823973), v[72], TOL);
    ASSERT_EQUAL_VEC(Vec3(       1.1116936,        0.5557477,     -0.083490269), v[73], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.58286536,        0.4940384,       0.15641722), v[75], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.11641305,        0.1695193,       0.12513766), v[76], TOL);
    ASSERT_EQUAL_VEC(Vec3(       0.4535942,        0.2204254,       -1.3901558), v[77], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.97038426,       0.54050455,        1.0798394), v[78], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.45515399,       0.26922391,      -0.81838449), v[80], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.25238199,      0.033498606,      -0.83556535), v[81], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.51678441,      -0.13017143,        3.3787147), v[82], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.45904589,       0.26453112,      -0.51102248), v[83], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.14722412,       0.07568654,      -0.22156668), v[85], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.81605979,      -0.19865477,        -0.266775), v[86], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.23665969,       0.54751246,        1.0088569), v[87], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.37752407,      -0.18085199,        -2.089343), v[88], TOL);
    ASSERT_EQUAL_VEC(Vec3(        0.380595,     -0.086932235,       0.67363332), v[90], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.14066565,      -0.28964581,       0.60712732), v[91], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.10235887,     -0.015072241,        2.9681353), v[92], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.59316071,      -0.14691037,      -0.77151166), v[93], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.28464737,       0.35107925,      -0.64974579), v[95], TOL);
    ASSERT_EQUAL_VEC(Vec3(        0.129441,       0.16055839,      -0.73824348), v[96], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.23882253,        0.4929275,        0.2715657), v[97], TOL);
    ASSERT_EQUAL_VEC(Vec3(     0.035100963,       0.22233184,       0.40791401), v[98], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.40030227,       0.50888684,      -0.26677383), v[100], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.55297771,       0.25903267,      -0.16886074), v[101], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.61181856,      -0.93279371,        1.1527765), v[102], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.90962545,       0.73132217,       0.35544777), v[103], TOL);
    ASSERT_EQUAL_VEC(Vec3(     0.038858756,       0.50588223,    -0.0079598324), v[105], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.24049741,       -0.2213647,      0.048064937), v[106], TOL);
    ASSERT_EQUAL_VEC(Vec3(    -0.091382572,       0.19851193,        1.5364968), v[107], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.10782317,       0.44417826,        1.1980488), v[108], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.14934041,      -0.23923058,       0.17286555), v[110], TOL);
    ASSERT_EQUAL_VEC(Vec3(    -0.066611113,      -0.99264022,      0.079926972), v[111], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.10303023,      -0.72206024,       0.97196124), v[112], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.72708003,       -0.2131092,       -1.2298523), v[113], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.10761756,      -0.34030673,      -0.53235022), v[115], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.34035525,       -1.1046927,       -0.5478589), v[116], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.26262173,      -0.40764237,        1.1786035), v[117], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.15472318,       -0.3614884,      -0.13415629), v[118], TOL);
    ASSERT_EQUAL_VEC(Vec3(     -0.85446696,       0.34659134,       0.83936102), v[120], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -1.2544238,      -0.12368226,       0.92987729), v[121], TOL);
    ASSERT_EQUAL_VEC(Vec3(      -1.2424201,        -1.949643,        2.2131077), v[122], TOL);
    ASSERT_EQUAL_VEC(Vec3(       1.7943001,       0.18094597,       -2.1686384), v[123], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.40146263,      -0.14102695,       -0.4071045), v[125], TOL);
    ASSERT_EQUAL_VEC(Vec3(     0.048158634,      -0.56057396,      -0.39355628), v[126], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.31440051,       -1.4227727,      -0.53995029), v[127], TOL);
    ASSERT_EQUAL_VEC(Vec3(       1.6724331,      -0.00827732,         1.004507), v[128], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.33843467,       -1.0216178,       0.37209769), v[130], TOL);
    ASSERT_EQUAL_VEC(Vec3(    -0.018097152,       -1.4187292,       0.24976993), v[131], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.33741779,       -1.0772242,       0.24605084), v[132], TOL);
    ASSERT_EQUAL_VEC(Vec3(      0.39823074,       -1.0175166,       0.83231671), v[133], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0018222952,   0.0031807229,    0.001323597), p[0], TOL);
    ASSERT_EQUAL_VEC(Vec3( -0.0012379191,   0.0044980753,    0.002644199), p[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(   0.097253311,  0.00048281192,  -0.0056017853), p[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(  -0.019228065,    0.096406261,   0.0066339568), p[3], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0097569797,    0.012837913,   0.0011513117), p[4], TOL);
    ASSERT_EQUAL_VEC(Vec3( -0.0012305501, -0.00044218817,     0.80172705), p[5], TOL);
    ASSERT_EQUAL_VEC(Vec3( -0.0016287879,  0.00016207453,      0.8012577), p[6], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.09302479,   0.0048340898,      0.7859027), p[7], TOL);
    ASSERT_EQUAL_VEC(Vec3(  -0.030836062,    0.090580612,     0.80090122), p[8], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0056660816,   0.0098306817,     0.79995087), p[9], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0024241203,  0.00025219236,      1.6005443), p[10], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0022792801,  -0.0029879077,      1.5990044), p[11], TOL);
    ASSERT_EQUAL_VEC(Vec3(   0.097973697,  -0.0040832734,      1.6042592), p[12], TOL);
    ASSERT_EQUAL_VEC(Vec3(  -0.017484864,    0.093833162,      1.6034688), p[13], TOL);
    ASSERT_EQUAL_VEC(Vec3(   0.010493211,   0.0097726101,      1.6012525), p[14], TOL);
    ASSERT_EQUAL_VEC(Vec3( -0.0013502532,     0.80165229,  0.00043054169), p[15], TOL);
    ASSERT_EQUAL_VEC(Vec3( -0.0039218748,       0.798544,  0.00069798185), p[16], TOL);
    ASSERT_EQUAL_VEC(Vec3(   0.094317285,     0.79996678,  -0.0022526919), p[17], TOL);
    ASSERT_EQUAL_VEC(Vec3(  -0.023618485,     0.89467115,    0.004163616), p[18], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0064797441,     0.81139543,  0.00054253525), p[19], TOL);
    ASSERT_EQUAL_VEC(Vec3(  -0.005133847,     0.79726965,     0.80223214), p[20], TOL);
    ASSERT_EQUAL_VEC(Vec3( -0.0075640917,     0.79573852,      0.8005436), p[21], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.09037011,     0.79136871,     0.80478002), p[22], TOL);
    ASSERT_EQUAL_VEC(Vec3(  -0.023450282,     0.89118779,     0.80472417), p[23], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0031002647,     0.80665903,     0.80276978), p[24], TOL);
    ASSERT_EQUAL_VEC(Vec3(-0.00012486134,     0.80321813,      1.5993745), p[25], TOL);
    ASSERT_EQUAL_VEC(Vec3(  4.813237e-05,     0.80026365,       1.598746), p[26], TOL);
    ASSERT_EQUAL_VEC(Vec3(   0.094244383,     0.78852317,      1.6057638), p[27], TOL);
    ASSERT_EQUAL_VEC(Vec3( -0.0097386928,     0.89844269,      1.6008495), p[28], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0089165681,     0.81180876,      1.6002135), p[29], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0014410392,       1.598711,  -0.0011258767), p[30], TOL);
    ASSERT_EQUAL_VEC(Vec3(-0.00098756021,      1.5967208,  0.00063645275), p[31], TOL);
    ASSERT_EQUAL_VEC(Vec3(   0.096891446,       1.597886,   0.0060055707), p[32], TOL);
    ASSERT_EQUAL_VEC(Vec3(  -0.021766997,       1.691572,    -0.00187122), p[33], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0091476185,      1.6085291, -0.00044462805), p[34], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0015647814,      1.5989168,      0.7987792), p[35], TOL);
    ASSERT_EQUAL_VEC(Vec3(-0.00070865882,      1.5988095,     0.79842537), p[36], TOL);
    ASSERT_EQUAL_VEC(Vec3(   0.096898851,      1.6046072,     0.80520987), p[37], TOL);
    ASSERT_EQUAL_VEC(Vec3(  -0.028503489,      1.6893951,     0.80725858), p[38], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0085271228,      1.6091758,     0.80036976), p[39], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0017073728,      1.5997787,      1.5993442), p[40], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0021442366,      1.5980895,      1.5972258), p[41], TOL);
    ASSERT_EQUAL_VEC(Vec3(   0.097270115,      1.5946811,      1.5973205), p[42], TOL);
    ASSERT_EQUAL_VEC(Vec3(  -0.017501011,      1.6932589,       1.591941), p[43], TOL);
    ASSERT_EQUAL_VEC(Vec3(  0.0098526054,       1.609207,      1.5983386), p[44], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80085021, -5.9180924e-05,  -0.0023143796), p[45], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79678957,   -0.002615688,  -0.0018961084), p[46], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.89644148,  3.0747659e-05,   0.0026473648), p[47], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.77651987,    0.092426548,   0.0017774772), p[48], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80845209,   0.0098164867,  -0.0013485711), p[49], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79716733,  -0.0028446492,     0.80005672), p[50], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79332157,  -0.0010797458,     0.80142981), p[51], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.89280431,  -0.0058722341,     0.80264894), p[52], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.77639763,    0.089870925,     0.78844809), p[53], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80515392,   0.0067229714,     0.79909488), p[54], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79746637,   0.0001884067,      1.5995574), p[55], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79631426,   0.0010217036,      1.5992362), p[56], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.89295697,   0.0018736447,      1.5931524), p[57], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.77177432,    0.092395895,      1.5994309), p[58], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80491225,    0.010204575,      1.5988606), p[59], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79866334,     0.79912882,   -0.004466586), p[60], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79764067,     0.79406125,  -0.0025853199), p[61], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.89424305,     0.80319359,   -0.001255052), p[62], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.77018574,     0.88848395,    0.014690841), p[63], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80582157,     0.80909455,  -0.0020803386), p[64], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79843536,     0.80003067,     0.80212268), p[65], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79445277,     0.79927926,     0.80224222), p[66], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.89409415,     0.80061518,     0.80549497), p[67], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.77378082,     0.89249621,      0.8042665), p[68], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80600986,     0.80995695,     0.80271112), p[69], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79978452,     0.80160908,      1.6005167), p[70], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79847858,     0.80007324,      1.5986022), p[71], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.89532462,     0.79583362,      1.6015424), p[72], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.78142937,     0.89554812,      1.5995876), p[73], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80801836,     0.81101408,       1.600527), p[74], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80285082,      1.6024665,  0.00075658535), p[75], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80189366,      1.6010496,   0.0016871161), p[76], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.89824818,      1.6009935,   -0.006956627), p[77], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.78058391,      1.6954436,   0.0054224398), p[78], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.81065214,      1.6122279,  0.00043150321), p[79], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80212713,      1.6012727,     0.79583605), p[80], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79816623,      1.6017889,     0.79738371), p[81], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.89539911,      1.5997901,      0.8172943), p[82], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.77857869,       1.694035,      0.7975505), p[83], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80956501,      1.6110101,     0.79830804), p[84], TOL);
    ASSERT_EQUAL_VEC(Vec3(     0.7992246,      1.6003456,      1.5989158), p[85], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.79816797,      1.5998807,      1.5979427), p[86], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.89472482,      1.6025392,      1.6050162), p[87], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.77366676,      1.6921049,       1.589464), p[88], TOL);
    ASSERT_EQUAL_VEC(Vec3(    0.80668582,      1.6103682,      1.5985583), p[89], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6018461, -0.00041508242,   0.0033721932), p[90], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6024467,  -0.0020852138,   0.0028925909), p[91], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6968702, -0.00036335958,    0.014893386), p[92], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5785034,    0.092129595,  -0.0039021765), p[93], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6094929,    0.009462798,   0.0038252304), p[94], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6014409,   0.0017635208,     0.79673787), p[95], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5993737, -0.00068558795,     0.79694993), p[96], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6970452,   0.0026916419,      0.8013502), p[97], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5762634,    0.093961686,     0.80202057), p[98], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6089538,    0.011697928,     0.79779344), p[99], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5979928,    0.002452089,      1.5987178), p[100], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5960886,   0.0040539637,      1.5969176), p[101], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6931939,  -0.0045515828,      1.6057901), p[102], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5805546,    0.096519704,      1.6018017), p[103], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6062883,    0.011739785,      1.5998012), p[104], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6001142,     0.80250915, -7.7946734e-06), p[105], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6006025,     0.79846627,  -0.0013688965), p[106], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6955081,     0.80123118,   0.0077822584), p[107], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5767785,     0.89514216,   0.0060660484), p[108], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6078011,     0.81225461,   0.0014711603), p[109], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6007701,     0.79882407,     0.80081999), p[110], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5985657,     0.79386367,     0.80234737), p[111], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6963725,     0.79626947,     0.80481689), p[112], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5795222,     0.89189174,     0.79380174), p[113], TOL);
    ASSERT_EQUAL_VEC(Vec3(      1.608702,     0.80847971,     0.80049768), p[114], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5994516,      0.7982202,      1.5972756), p[115], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5973651,      0.7975011,      1.5995656), p[116], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6947817,     0.79802711,      1.6059044), p[117], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5753545,     0.89083413,      1.5993514), p[118], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6070505,     0.80807935,      1.5984175), p[119], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5956633,      1.6016008,   0.0042085051), p[120], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5958491,      1.6011046,   0.0043274523), p[121], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6903928,      1.5897705,    0.011186024), p[122], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5842502,       1.695412,   -0.011007606), p[123], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6045512,      1.6103463,   0.0033296391), p[124], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6019429,      1.5992615,     0.79795125), p[125], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6024851,      1.5975188,     0.79844451), p[126], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6974452,      1.5928404,      0.7973163), p[127], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5842436,      1.6930615,     0.80506611), p[128], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6102427,      1.6085828,      0.7986425), p[129], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6017066,      1.5948476,      1.6018311), p[130], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5993111,      1.5947857,       1.602598), p[131], TOL);
    ASSERT_EQUAL_VEC(Vec3(      1.697424,       1.594611,      1.6011635), p[132], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.5779522,      1.6875457,      1.6040914), p[133], TOL);
    ASSERT_EQUAL_VEC(Vec3(     1.6093834,      1.6047111,       1.602001), p[134], TOL);
}

int main(int argc, char* argv[]) {
    try {
        registerKernelFactories();

        Platform& platform = initializePlatform(argc, argv);
        testWaterBox(platform);
        double maxDrudeDistance = 0.005;
        double observedDrudeDistance = testWaterBoxWithHardWallConstraint(platform, 0.0);
        ASSERT(observedDrudeDistance > maxDrudeDistance);
        observedDrudeDistance = testWaterBoxWithHardWallConstraint(platform, maxDrudeDistance);
        ASSERT(observedDrudeDistance < maxDrudeDistance);
        testPositionsAfterShortRun(platform);
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

