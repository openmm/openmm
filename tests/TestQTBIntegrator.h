/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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
#include "openmm/CustomNonbondedForce.h"
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

    int numFrequencies = 10;
    int numReplicas = 10;
    int numParticles = numFrequencies*numReplicas;
    double temperature = 300.0;
    double mass = 1.0;
    System system;
    vector<Vec3> positions(numParticles);
    vector<double> k;
    CustomExternalForce* force = new CustomExternalForce("0.5*k*(x*x+y*y+z*z)");
    system.addForce(force);
    force->addPerParticleParameter("k");
    QTBIntegrator integrator(temperature, 10.0, 0.001);
    integrator.setCutoffFrequency(500.0);
    for (int j = 0; j < numReplicas; j++) {
        int base = system.getNumParticles();
        for (int i = 0; i < numFrequencies; i++) {
            system.addParticle(mass);
            k.push_back(1000.0*(i+1)*(i+1));
            force->addParticle(i+base, {k[i]});
            integrator.setParticleType(i+base, i);
        }
    }
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temperature);

    // Compute the average energy of each particle over a simulation.

    integrator.step(10000);
    vector<double> energy(numParticles, 0.0);
    int numSteps = 10000;
    for (int i = 0; i < numSteps; i++) {
        integrator.step(10);
        State state = context.getState(State::Positions);
        for (int j = 0; j < numParticles; j++) {
            Vec3 p = state.getPositions()[j];
            energy[j%numFrequencies] += 0.5*k[j]*p.dot(p);
        }
    }
    for (int i = 0; i < numParticles; i++)
        energy[i] /= numSteps*numReplicas;

    // Compare to the expected distribution.

    for (int i = 0; i < numFrequencies; i++) {
        double w = sqrt(k[i]/mass);
        double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
        double kT = BOLTZ*temperature;
        double expected = 1.5*hbar*w*(0.5+1/(exp(hbar*w/kT)-1));
        ASSERT_USUALLY_EQUAL_TOL(expected, energy[i], 0.07);
    }
}

void testCoupledHarmonic() {
    // Create a collection of weakly coupled harmonic oscillators.

    int numFrequencies = 8;
    int numReplicas = 10;
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
    integrator.setDefaultAdaptationRate(1.0);
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

    // Equilibrate with a high adaptation rate to let the spectrum converge.

    for (int i = 0; i < 15; i++) {
        integrator.step(10000);
        // OpenCL on Mac hangs if we go too long without downloading anything.  I don't know why.
        context.getState(State::Positions);
    }
    integrator.setDefaultAdaptationRate(0.05);
    context.reinitialize(true);

    // Compute the average energy of each particle over a simulation.

    vector<double> energy(numFrequencies, 0.0);
    int numSteps = 20000;
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
        ASSERT_USUALLY_EQUAL_TOL(expected, energy[i], 0.15);
    }
}

void testParaHydrogen() {
    const int numParticles = 32;
    const double temperature = 25.0;
    const double mass = 2.0;
    const double boxSize = 1.1896;
    const int numSteps = 2000;
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
    QTBIntegrator integ(temperature, 40.0, 0.001);
    for (int i = 0; i < numParticles; i++)
        integ.setParticleType(i, 0);
    integ.setDefaultAdaptationRate(0.5);
    integ.setSegmentLength(0.5);
    Context context(system, integ, platform);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++)
        positions[i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
    context.setPositions(positions);
    integ.step(50000);

    // Simulate it.

    vector<int> counts(numBins, 0);
    const double invBoxSize = 1.0/boxSize;
    for (int step = 0; step < numSteps; step++) {
        integ.step(20);
        State state = context.getState(State::Positions);

        // Record the radial distribution function.

        const vector<Vec3>& pos = state.getPositions();
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
    }

    // Check against expected values.

    double scale = (boxSize*boxSize*boxSize)/(numSteps*0.5*numParticles*numParticles);
    for (int i = 0; i < numBins/2; i++) {
        double r1 = i*boxSize/numBins;
        double r2 = (i+1)*boxSize/numBins;
        double volume = (4.0/3.0)*M_PI*(r2*r2*r2-r1*r1*r1);
        ASSERT_USUALLY_EQUAL_TOL(reference[i], scale*counts[i]/volume, 0.2);
    }
}

void testConstraints() {
    const int numParticles = 8;
    const int numConstraints = 5;
    const double temp = 100.0;
    System system;
    QTBIntegrator integrator(temp, 2.0, 0.01);
    integrator.setConstraintTolerance(1e-5);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(10.0);
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

    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Positions);
        for (int j = 0; j < numConstraints; ++j) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(j, particle1, particle2, distance);
            Vec3 p1 = state.getPositions()[particle1];
            Vec3 p2 = state.getPositions()[particle2];
            double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
            ASSERT_EQUAL_TOL(distance, dist, 1e-4);
        }
        integrator.step(1);
    }
}

void testConstrainedMasslessParticles() {
    System system;
    system.addParticle(0.0);
    system.addParticle(1.0);
    system.addConstraint(0, 1, 1.5);
    vector<Vec3> positions(2);
    positions[0] = Vec3(-1, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    QTBIntegrator integrator(300.0, 2.0, 0.01);
    bool failed = false;
    try {
        // This should throw an exception.
        
        Context context(system, integrator, platform);
    }
    catch (exception& ex) {
        failed = true;
    }
    ASSERT(failed);
    
    // Now make both particles massless, which should work.
    
    system.setParticleMass(1, 0.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(300.0);
    integrator.step(1);
    State state = context.getState(State::Velocities);
    ASSERT_EQUAL(0.0, state.getVelocities()[0][0]);
}

void testRandomSeed() {
    const int numParticles = 8;
    const double temp = 100.0;
    System system;
    QTBIntegrator integrator(temp, 2.0, 0.01);
    NonbondedForce* forceField = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(2.0);
        forceField->addParticle((i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(forceField);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
        velocities[i] = Vec3(0, 0, 0);
    }

    // Try twice with the same random seed.

    integrator.setRandomNumberSeed(5);
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

    integrator.setRandomNumberSeed(10);
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
            ASSERT_EQUAL_TOL(state1.getPositions()[i][j], state2.getPositions()[i][j], 1e-6);
            ASSERT_EQUAL_TOL(state3.getPositions()[i][j], state4.getPositions()[i][j], 1e-6);
            ASSERT(state1.getPositions()[i][j] != state3.getPositions()[i][j]);
        }
    }
}

void testInitialTemperature() {
    // Check temperature initialization for a collection of randomly placed particles
    const int numParticles = 50000;
    const int nDoF = 3 * numParticles;
    const double targetTemperature = 300;
    System system;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    std::vector<Vec3> positions(numParticles);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions[i][0] = genrand_real2(sfmt);
        positions[i][1] = genrand_real2(sfmt);
        positions[i][2] = genrand_real2(sfmt);
    }

    QTBIntegrator integrator(300, 25, 0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(targetTemperature);
    auto velocities = context.getState(State::Velocities).getVelocities();
    double kineticEnergy = 0;
    for(const auto &v : velocities) kineticEnergy += 0.5 * v.dot(v);
    double temperature = (2*kineticEnergy / (nDoF*BOLTZ));
    ASSERT_USUALLY_EQUAL_TOL(targetTemperature, temperature, 0.01);
}

void testSerializeParameters() {
    // Test serializing the integrator's parameters.

    int numParticles = 10;
    double temperature = 300.0;
    double mass = 2.0;
    System system;
    vector<Vec3> positions(numParticles);
    CustomExternalForce* force = new CustomExternalForce("0.5*k*(x*x+y*y+z*z)");
    system.addForce(force);
    force->addPerParticleParameter("k");
    QTBIntegrator integrator(temperature, 10.0, 0.001);
    integrator.setCutoffFrequency(500.0);
    integrator.setDefaultAdaptationRate(0.1);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(mass);
        force->addParticle(i, {1000.0*(i+1)*(i+1)});
        if (i < 5)
            integrator.setParticleType(i, i%3);
    }
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(temperature);

    // Run for a little while, then record a State.

    integrator.step(5000);
    State state = context.getState(State::IntegratorParameters);

    // Create a new Integrator and Context, set the State, and see if the adapted
    // friction coefficients were set correctly.

    QTBIntegrator integrator2(temperature, 10.0, 0.001);
    for (auto type : integrator.getParticleTypes())
        integrator2.setParticleType(type.first, type.second);
    Context context2(system, integrator2, platform);
    context2.setPositions(positions);
    context2.setState(state);
    for (int i = 0; i < numParticles; i++) {
        vector<double> f1, f2;
        integrator.getAdaptedFriction(i, f1);
        integrator2.getAdaptedFriction(i, f2);
        ASSERT_EQUAL_CONTAINERS(f1, f2);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testHarmonic();
        testCoupledHarmonic();
        testParaHydrogen();
        testConstraints();
        testConstrainedMasslessParticles();
        testRandomSeed();
        testInitialTemperature();
        testSerializeParameters();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
