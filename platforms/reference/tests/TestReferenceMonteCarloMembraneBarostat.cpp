/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
 * This tests the reference implementation of MonteCarloMembraneBarostat.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/MonteCarloMembraneBarostat.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/HarmonicBondForce.h"
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

ReferencePlatform platform;

void testIdealGas(MonteCarloMembraneBarostat::XYMode xymode, MonteCarloMembraneBarostat::ZMode zmode) {
    const int numParticles = 64;
    const int frequency = 1;
    const int steps = 5000;
    const double pressure = 1.5;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double tension = (zmode == MonteCarloMembraneBarostat::ZFixed ? 0.2 : 0.0);
    const double tensionInMD = tension*(AVOGADRO*1e-25); // surface tension in kJ/mol/nm^2
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
    MonteCarloMembraneBarostat* barostat = new MonteCarloMembraneBarostat(pressure, tension, temp[0], xymode, zmode, frequency);
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

        integrator.step(1000);

        // Now run it for a while and see if the volume is correct.

        double volume = 0.0, zsize = 0.0;
        for (int j = 0; j < steps; ++j) {
            Vec3 box[3];
            context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
            volume += box[0][0]*box[1][1]*box[2][2];
            zsize += box[2][2];
            if (xymode == MonteCarloMembraneBarostat::XYIsotropic)
                ASSERT_EQUAL_TOL(0.5*box[0][0], box[1][1], 1e-5);
            if (zmode == MonteCarloMembraneBarostat::ZFixed)
                ASSERT_EQUAL_TOL(2*initialLength, box[2][2], 1e-5);
            if (zmode == MonteCarloMembraneBarostat::ConstantVolume)
                ASSERT_EQUAL_TOL(initialVolume, box[0][0]*box[1][1]*box[2][2], 1e-5);
            integrator.step(frequency);
        }
        volume /= steps;
        zsize /= steps;
        if (zmode != MonteCarloMembraneBarostat::ConstantVolume) {
            double effectivePressure = pressureInMD-tensionInMD/zsize;
            double expected = (numParticles+1)*BOLTZ*temp[i]/effectivePressure;
            ASSERT_USUALLY_EQUAL_TOL(expected, volume, 0.05);
        }
    }
}

void testRandomSeed() {
    const int numParticles = 8;
    const double temp = 100.0;
    const double pressure = 1.5;
    const double tension = 0.3;
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
    MonteCarloMembraneBarostat* barostat = new MonteCarloMembraneBarostat(pressure, tension, temp, MonteCarloMembraneBarostat::XYAnisotropic, MonteCarloMembraneBarostat::ZFree, 1);
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

int main() {
    try {
        testIdealGas(MonteCarloMembraneBarostat::XYIsotropic, MonteCarloMembraneBarostat::ZFree);
        testIdealGas(MonteCarloMembraneBarostat::XYIsotropic, MonteCarloMembraneBarostat::ZFixed);
        testIdealGas(MonteCarloMembraneBarostat::XYIsotropic, MonteCarloMembraneBarostat::ConstantVolume);
        testIdealGas(MonteCarloMembraneBarostat::XYAnisotropic, MonteCarloMembraneBarostat::ZFree);
        testIdealGas(MonteCarloMembraneBarostat::XYAnisotropic, MonteCarloMembraneBarostat::ZFixed);
        testIdealGas(MonteCarloMembraneBarostat::XYAnisotropic, MonteCarloMembraneBarostat::ConstantVolume);
        testRandomSeed();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
