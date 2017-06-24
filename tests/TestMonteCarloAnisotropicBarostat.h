/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
#include "openmm/CustomExternalForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/MonteCarloAnisotropicBarostat.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
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
    const double pressure = 1.5;
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
    MonteCarloAnisotropicBarostat* barostat = new MonteCarloAnisotropicBarostat(Vec3(pressure, pressure, pressure), temp[0], true, true, true, frequency);
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
        
        integrator.step(10000);
        
        // Now run it for a while and see if the volume is correct.
        
        double volume = 0.0;
        for (int j = 0; j < steps; ++j) {
            Vec3 box[3];
            context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
            volume += box[0][0]*box[1][1]*box[2][2];
            integrator.step(frequency);
        }
        volume /= steps;
        double expected = (numParticles+1)*BOLTZ*temp[i]/pressureInMD;
        ASSERT_USUALLY_EQUAL_TOL(expected, volume, 3/std::sqrt((double) steps));
    }
}

void testIdealGasAxis(int axis) {
    // Test scaling just one axis.
    const int numParticles = 64;
    const int frequency = 10;
    const int steps = 1000;
    const double pressure = 1.5;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double temp[] = {300.0, 600.0, 1000.0};
    const double initialVolume = numParticles*BOLTZ*temp[1]/pressureInMD;
    const double initialLength = std::pow(initialVolume, 1.0/3.0);
    const bool scaleX = (axis == 0);
    const bool scaleY = (axis == 1);
    const bool scaleZ = (axis == 2);
    double boxX;
    double boxY;
    double boxZ;
    
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
    MonteCarloAnisotropicBarostat* barostat = new MonteCarloAnisotropicBarostat(Vec3(pressure, pressure, pressure), temp[0], scaleX, scaleY, scaleZ, frequency);
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
        
        integrator.step(10000);
        
        // Now run it for a while and see if the volume is correct.
        
        double volume = 0.0;
        for (int j = 0; j < steps; ++j) {
            Vec3 box[3];
            context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
            boxX = box[0][0];
            boxY = box[1][1];
            boxZ = box[2][2];
            volume += box[0][0]*box[1][1]*box[2][2];
            integrator.step(frequency);
        }
        volume /= steps;
        double expected = (numParticles+1)*BOLTZ*temp[i]/pressureInMD;
        ASSERT_USUALLY_EQUAL_TOL(expected, volume, 3/std::sqrt((double) steps));
        if (!scaleX) {
            ASSERT(boxX == initialLength);
        }
        if (!scaleY) {
            ASSERT(boxY == 0.5*initialLength);
        }
        if (!scaleZ) {
            ASSERT(boxZ == 2*initialLength);
        }
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
    MonteCarloAnisotropicBarostat* barostat = new MonteCarloAnisotropicBarostat(Vec3(pressure, pressure, pressure), temp, true, true, true, 1);
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

void testTriclinic() {
    const int numParticles = 64;
    const int frequency = 10;
    const int steps = 1000;
    const double pressure = 1.5;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double temperature = 300.0;
    const double initialVolume = numParticles*BOLTZ*temperature/pressureInMD;
    const double initialLength = std::pow(initialVolume, 1.0/3.0);

    // Create a gas of noninteracting particles.

    System system;
    Vec3 initialBox[3];
    initialBox[0] = Vec3(initialLength, 0, 0);
    initialBox[1] = Vec3(0.2*initialLength, initialLength, 0);
    initialBox[2] = Vec3(0.1*initialLength, 0.3*initialLength, initialLength);
    system.setDefaultPeriodicBoxVectors(initialBox[0], initialBox[1], initialBox[2]);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(1.0);
        positions[i] = Vec3(initialLength*genrand_real2(sfmt), initialLength*genrand_real2(sfmt), initialLength*genrand_real2(sfmt));
    }
    MonteCarloAnisotropicBarostat* barostat = new MonteCarloAnisotropicBarostat(Vec3(pressure, pressure, pressure), temperature, true, true, true, frequency);
    system.addForce(barostat);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->setUsesPeriodicBoundaryConditions(true);
    system.addForce(bonds); // So it won't complain the system is non-periodic.

    // Run a simulation

    LangevinIntegrator integrator(temperature, 0.1, 0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    // Let it equilibrate.

    integrator.step(10000);

    // Now run it for a while and see if the volume is correct.

    double volume = 0.0;
    for (int j = 0; j < steps; ++j) {
        Vec3 box[3];
        context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
        volume += box[0][0]*box[1][1]*box[2][2];
        integrator.step(frequency);
    }
    volume /= steps;
    double expected = (numParticles+1)*BOLTZ*temperature/pressureInMD;
    ASSERT_USUALLY_EQUAL_TOL(expected, volume, 3/std::sqrt((double) steps));

    // Make sure the box vectors have been scaled consistently.

    State state = context.getState(State::Positions);
    Vec3 box[3];
    state.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double xscale = box[2][0]/(0.1*initialLength);
    double yscale = box[2][1]/(0.3*initialLength);
    double zscale = box[2][2]/(1.0*initialLength);
    for (int i = 0; i < 3; i++) {
        ASSERT_EQUAL_VEC(Vec3(xscale*initialBox[i][0], yscale*initialBox[i][1], zscale*initialBox[i][2]), box[i], 1e-5);
    }

    // The barostat should have put all particles inside the first periodic box.  One integration step
    // has happened since then, so they may have moved slightly outside it.

    for (int i = 0; i < numParticles; i++) {
        Vec3 pos = state.getPositions()[i];
        ASSERT(pos[2]/box[2][2] > -1 && pos[2]/box[2][2] < 2);
        pos -= box[2]*floor(pos[2]/box[2][2]);
        ASSERT(pos[1]/box[1][1] > -1 && pos[1]/box[1][1] < 2);
        pos -= box[1]*floor(pos[1]/box[1][1]);
        ASSERT(pos[0]/box[0][0] > -1 && pos[0]/box[0][0] < 2);
    }
}

/**
 * Run a constant pressure simulation on an anisotropic Einstein crystal
 * using isotropic and anisotropic barostats.  There are a total of 15 simulations:
 *
 * 1) 3 pressures: 9.0, 10.0, 11.0 bar, for each of the following groups:
 * 2) 3 groups of simulations that scale just one axis: x, y, z
 * 3) 1 group of simulations that scales all three axes in the anisotropic barostat
 * 4) 1 group of simulations that scales all three axes in the isotropic barostat
 *
 * Results that we will check:
 *
 * a) In each group of simulations, the volume should decrease with increasing pressure
 * b) In the three simulation groups that scale just one axis, the compressibility (i.e. incremental volume change
 * with increasing pressure) should go like kx > ky > kz (because the spring constant is largest in the z-direction)
 * c) The anisotropic barostat should produce the same result as the isotropic barostat when all three axes are scaled
 */
void testEinsteinCrystal() {
    const int numParticles = 64;
    const int frequency = 2;
    const int equil = 10000;
    const int steps = 5000;
    const double pressure = 10.0;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double temp = 300.0; // Only test one temperature since we're looking at three pressures.
    const double pres3[] = {2.0, 8.0, 15.0};
    const double initialVolume = numParticles*BOLTZ*temp/pressureInMD;
    const double initialLength = std::pow(initialVolume, 1.0/3.0);
    vector<double> initialPositions(3);
    vector<double> results;
    // Run four groups of anisotropic simulations; scaling just x, y, z, then all three.
    for (int a = 0; a < 4; a++) {
        // Test barostat for three different pressures.
        for (int p = 0; p < 3; p++) {
            // Create a system of noninteracting particles attached by harmonic springs to their initial positions.
            System system;
            system.setDefaultPeriodicBoxVectors(Vec3(initialLength, 0, 0), Vec3(0, initialLength, 0), Vec3(0, 0, initialLength));
            vector<Vec3> positions(numParticles);
            OpenMM_SFMT::SFMT sfmt;
            init_gen_rand(0, sfmt);
            // Anisotropic force constants.
            CustomExternalForce* force = new CustomExternalForce("0.005*(x-x0)^2 + 0.01*(y-y0)^2 + 0.02*(z-z0)^2");
            force->addPerParticleParameter("x0");
            force->addPerParticleParameter("y0");
            force->addPerParticleParameter("z0");
            NonbondedForce* nb = new NonbondedForce();
            nb->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
            for (int i = 0; i < numParticles; ++i) {
                system.addParticle(1.0);
                positions[i] = Vec3(((i/16)%4+0.5)*initialLength/4, ((i/4)%4+0.5)*initialLength/4, (i%4+0.5)*initialLength/4);
                initialPositions[0] = positions[i][0];
                initialPositions[1] = positions[i][1];
                initialPositions[2] = positions[i][2];
                force->addParticle(i, initialPositions);
                nb->addParticle(0, initialLength/6, 0.1);
            }
            system.addForce(force);
            system.addForce(nb);
            // Create the barostat.
            MonteCarloAnisotropicBarostat* barostat = new MonteCarloAnisotropicBarostat(Vec3(pres3[p], pres3[p], pres3[p]), temp, (a==0||a==3), (a==1||a==3), (a==2||a==3), frequency);
            system.addForce(barostat);
            barostat->setDefaultTemperature(temp);
            LangevinIntegrator integrator(temp, 0.1, 0.01);
            Context context(system, integrator, platform);
            context.setPositions(positions);
            // Let it equilibrate.
            integrator.step(equil);
            // Now run it for a while and see if the volume is correct.
            double volume = 0.0;
            for (int j = 0; j < steps; ++j) {
                Vec3 box[3];
                context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
                volume += box[0][0]*box[1][1]*box[2][2];
                integrator.step(frequency);
            }
            volume /= steps;
            results.push_back(volume);
        }
    }
    for (int p = 0; p < 3; p++) {
        // Create a system of noninteracting particles attached by harmonic springs to their initial positions.
        System system;
        system.setDefaultPeriodicBoxVectors(Vec3(initialLength, 0, 0), Vec3(0, initialLength, 0), Vec3(0, 0, initialLength));
        vector<Vec3> positions(numParticles);
        OpenMM_SFMT::SFMT sfmt;
        init_gen_rand(0, sfmt);
        // Anisotropic force constants.
        CustomExternalForce* force = new CustomExternalForce("0.005*(x-x0)^2 + 0.01*(y-y0)^2 + 0.02*(z-z0)^2");
        force->addPerParticleParameter("x0");
        force->addPerParticleParameter("y0");
        force->addPerParticleParameter("z0");
        NonbondedForce* nb = new NonbondedForce();
        nb->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
        for (int i = 0; i < numParticles; ++i) {
            system.addParticle(1.0);
            positions[i] = Vec3(((i/16)%4+0.5)*initialLength/4, ((i/4)%4+0.5)*initialLength/4, (i%4+0.5)*initialLength/4);
            initialPositions[0] = positions[i][0];
            initialPositions[1] = positions[i][1];
            initialPositions[2] = positions[i][2];
            force->addParticle(i, initialPositions);
            nb->addParticle(0, initialLength/6, 0.1);
        }
        system.addForce(force);
        system.addForce(nb);
        // Create the barostat.
        MonteCarloBarostat* barostat = new MonteCarloBarostat(pres3[p], temp, frequency);
        system.addForce(barostat);
        barostat->setDefaultTemperature(temp);
        LangevinIntegrator integrator(temp, 0.1, 0.001);
        Context context(system, integrator, platform);
        context.setPositions(positions);
        // Let it equilibrate.
        integrator.step(equil);
        // Now run it for a while and see if the volume is correct.
        double volume = 0.0;
        for (int j = 0; j < steps; ++j) {
            Vec3 box[3];
            context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
            volume += box[0][0]*box[1][1]*box[2][2];
            integrator.step(frequency);
        }
        volume /= steps;
        results.push_back(volume);
    }
    
    // Check to see if volumes decrease with increasing pressure.
    ASSERT_USUALLY_TRUE(results[0] > results[1]);
    ASSERT_USUALLY_TRUE(results[1] > results[2]);
    ASSERT_USUALLY_TRUE(results[3] > results[4]);
    ASSERT_USUALLY_TRUE(results[4] > results[5]);
    ASSERT_USUALLY_TRUE(results[6] > results[7]);
    ASSERT_USUALLY_TRUE(results[7] > results[8]);

    // Check to see if incremental volume changes with increasing pressure go like kx > ky > kz.
    ASSERT_USUALLY_TRUE((results[0] - results[1]) > (results[3] - results[4]));
    ASSERT_USUALLY_TRUE((results[1] - results[2]) > (results[4] - results[5]));
    ASSERT_USUALLY_TRUE((results[3] - results[4]) > (results[6] - results[7]));
    ASSERT_USUALLY_TRUE((results[4] - results[5]) > (results[7] - results[8]));
    
    // Check to see if the volumes are equal for isotropic and anisotropic (all axis).
    ASSERT_USUALLY_EQUAL_TOL(results[9], results[12], 3/std::sqrt((double) steps));
    ASSERT_USUALLY_EQUAL_TOL(results[10], results[13], 3/std::sqrt((double) steps));
    ASSERT_USUALLY_EQUAL_TOL(results[11], results[14], 3/std::sqrt((double) steps));
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testIdealGas();
        testIdealGasAxis(0);
        testIdealGasAxis(1);
        testIdealGasAxis(2);
        testRandomSeed();
        testTriclinic();
        //testEinsteinCrystal();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

