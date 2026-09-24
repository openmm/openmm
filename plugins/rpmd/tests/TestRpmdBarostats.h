/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2026 Stanford University and the Authors.           *
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
 * This tests barostats for RPMD.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloAnisotropicBarostat.h"
#include "openmm/RPMDMonteCarloBarostat.h"
#include "openmm/RPMDMonteCarloFlexibleBarostat.h"
#include "openmm/RPMDMonteCarloMembraneBarostat.h"
#include "SimTKOpenMMUtilities.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;


void testIdealGas() {
    const int numCopies = 3;
    const int numParticles = 64;
    const int frequency = 1;
    const int steps = 1000;
    const double pressure = 1.5;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double temp[] = {300.0, 600.0};
    const double initialVolume = numParticles*BOLTZ*0.5*(temp[0]+temp[1])/pressureInMD;
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
    RPMDMonteCarloBarostat* barostat = new RPMDMonteCarloBarostat(pressure, frequency);
    system.addForce(barostat);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->setUsesPeriodicBoundaryConditions(true);
    system.addForce(bonds); // So it won't complain the system is non-periodic.

    // Test it for two different temperatures.

    for (int i = 0; i < 2; i++) {
        RPMDIntegrator integrator(numCopies, temp[i], 0.1, 0.01);
        Context context(system, integrator, platform);
        for (int copy = 0; copy < numCopies; copy++)
            integrator.setPositions(copy, positions);

        // Let it equilibrate.

        integrator.step(1000);

        // Now run it for a while and see if the volume is correct.

        double volume = 0.0;
        for (int j = 0; j < steps; ++j) {
            Vec3 box[3];
            integrator.getState(0, 0).getPeriodicBoxVectors(box[0], box[1], box[2]);
            volume += box[0][0]*box[1][1]*box[2][2];
            ASSERT_EQUAL_TOL(0.5*box[0][0], box[1][1], 1e-5);
            ASSERT_EQUAL_TOL(2*box[0][0], box[2][2], 1e-5);
            integrator.step(frequency);
        }
        volume /= steps;
        double expected = (numParticles+1)*BOLTZ*temp[i]/pressureInMD;
        ASSERT_USUALLY_EQUAL_TOL(expected, volume, 0.05);
    }
}

void testAnisotropicIdealGas() {
    const int numCopies = 3;
    const int numParticles = 64;
    const int frequency = 1;
    const int steps = 1000;
    const double pressure = 3.0;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double temp = 300.0;
    const double initialVolume = numParticles*BOLTZ*temp/pressureInMD;
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
    RPMDMonteCarloAnisotropicBarostat* barostat = new RPMDMonteCarloAnisotropicBarostat(Vec3(pressure, pressure, pressure), true, true, true, frequency);
    system.addForce(barostat);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->setUsesPeriodicBoundaryConditions(true);
    system.addForce(bonds); // So it won't complain the system is non-periodic.

    // Try simulating it.

    RPMDIntegrator integrator(numCopies, temp, 0.1, 0.01);
    Context context(system, integrator, platform);
    for (int copy = 0; copy < numCopies; copy++)
        integrator.setPositions(copy, positions);

    // Let it equilibrate.

    integrator.step(1000);

    // Now run it for a while and see if the volume is correct.

    double volume = 0.0;
    for (int j = 0; j < steps; ++j) {
        Vec3 box[3];
        context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
        volume += box[0][0]*box[1][1]*box[2][2];
        integrator.step(frequency);
    }
    volume /= steps;
    double expected = (numParticles+1)*BOLTZ*temp/pressureInMD;
    ASSERT_USUALLY_EQUAL_TOL(expected, volume, 0.05);
}

void testIdealGasAxis(int axis) {
    const int numCopies = 3;
    const int numParticles = 64;
    const int frequency = 1;
    const int steps = 1000;
    const double pressure = 3.0;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double temp = 300.0;
    const double initialVolume = numParticles*BOLTZ*temp/pressureInMD;
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
    RPMDMonteCarloAnisotropicBarostat* barostat = new RPMDMonteCarloAnisotropicBarostat(Vec3(pressure, pressure, pressure), scaleX, scaleY, scaleZ, frequency);
    system.addForce(barostat);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->setUsesPeriodicBoundaryConditions(true);
    system.addForce(bonds); // So it won't complain the system is non-periodic.

    // Try simulating it.

    RPMDIntegrator integrator(numCopies, temp, 0.1, 0.01);
    Context context(system, integrator, platform);
    for (int copy = 0; copy < numCopies; copy++)
        integrator.setPositions(copy, positions);

    // Let it equilibrate.

    integrator.step(1000);

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
    double expected = (numParticles+1)*BOLTZ*temp/pressureInMD;
    ASSERT_USUALLY_EQUAL_TOL(expected, volume, 0.05);
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

void testFlexibleIdealGas() {
    const int numCopies = 3;
    const int numParticles = 64;
    const int frequency = 1;
    const int steps = 1000;
    const double pressure = 3.0;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double temp = 300.0;
    const double initialVolume = numParticles*BOLTZ*temp/pressureInMD;
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
    RPMDMonteCarloFlexibleBarostat* barostat = new RPMDMonteCarloFlexibleBarostat(pressure, frequency);
    system.addForce(barostat);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->setUsesPeriodicBoundaryConditions(true);
    system.addForce(bonds); // So it won't complain the system is non-periodic.

    // Try simulating it.

    RPMDIntegrator integrator(numCopies, temp, 0.1, 0.01);
    Context context(system, integrator, platform);
    for (int copy = 0; copy < numCopies; copy++)
        integrator.setPositions(copy, positions);

    // Let it equilibrate.

    integrator.step(1000);

    // Now run it for a while and see if the volume is correct.

    double volume = 0.0;
    for (int j = 0; j < steps; ++j) {
        Vec3 box[3];
        context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
        volume += box[0][0]*box[1][1]*box[2][2];
        integrator.step(frequency);
    }
    volume /= steps;
    double expected = (numParticles+1)*BOLTZ*temp/pressureInMD;
    ASSERT_USUALLY_EQUAL_TOL(expected, volume, 0.05);
}

void testMembraneIdealGas(RPMDMonteCarloMembraneBarostat::XYMode xymode, RPMDMonteCarloMembraneBarostat::ZMode zmode) {
    const int numCopies = 3;
    const int numParticles = 64;
    const int frequency = 1;
    const int steps = 1000;
    const double pressure = 1.5;
    const double pressureInMD = pressure*(AVOGADRO*1e-25); // pressure in kJ/mol/nm^3
    const double tension = (zmode == RPMDMonteCarloMembraneBarostat::ZFixed ? 0.2 : 0.0);
    const double tensionInMD = tension*(AVOGADRO*1e-25); // surface tension in kJ/mol/nm^2
    const double temp = 300.0;
    const double initialVolume = numParticles*BOLTZ*temp/pressureInMD;
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
    RPMDMonteCarloMembraneBarostat* barostat = new RPMDMonteCarloMembraneBarostat(pressure, tension, xymode, zmode, frequency);
    system.addForce(barostat);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->setUsesPeriodicBoundaryConditions(true);
    system.addForce(bonds); // So it won't complain the system is non-periodic.

    // Try simulating it.

    RPMDIntegrator integrator(numCopies, temp, 0.1, 0.01);
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
        if (xymode == RPMDMonteCarloMembraneBarostat::XYIsotropic)
            ASSERT_EQUAL_TOL(0.5*box[0][0], box[1][1], 1e-5);
        if (zmode == RPMDMonteCarloMembraneBarostat::ZFixed)
            ASSERT_EQUAL_TOL(2*initialLength, box[2][2], 1e-5);
        if (zmode == RPMDMonteCarloMembraneBarostat::ConstantVolume)
            ASSERT_EQUAL_TOL(initialVolume, box[0][0]*box[1][1]*box[2][2], 1e-5);
        integrator.step(frequency);
    }
    volume /= steps;
    zsize /= steps;
    if (zmode != RPMDMonteCarloMembraneBarostat::ConstantVolume) {
        double effectivePressure = pressureInMD-tensionInMD/zsize;
        double expected = (numParticles+1)*BOLTZ*temp/effectivePressure;
        ASSERT_USUALLY_EQUAL_TOL(expected, volume, 0.05);
    }
}

void testWater() {
    const int numCopies = 8;
    const int gridSize = 8;
    const int numMolecules = gridSize*gridSize*gridSize;
    const int frequency = 10;
    const int steps = 400;
    const double temp = 273.15;
    const double pressure = 3;
    const double spacing = 0.31;
    const double angle = 112*M_PI/180;
    const double dOH = 0.1;

    // Create a box of q-SPC/Fw water molecules.

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(gridSize*spacing, 0, 0), Vec3(0, gridSize*spacing, 0), Vec3(0, 0, gridSize*spacing));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setUseDispersionCorrection(true);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    HarmonicAngleForce* angles = new HarmonicAngleForce();
    vector<Vec3> positions;
    Vec3 offset1(dOH, 0, 0);
    Vec3 offset2(dOH*std::cos(angle), dOH*std::sin(angle), 0);
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            for (int k = 0; k < gridSize; ++k) {
                int firstParticle = system.getNumParticles();
                system.addParticle(16.0);
                system.addParticle(1.0);
                system.addParticle(1.0);
                nonbonded->addParticle(-0.84, 0.3165492, 0.650143);
                nonbonded->addParticle(0.42, 1, 0);
                nonbonded->addParticle(0.42, 1, 0);
                Vec3 pos = Vec3(spacing*i, spacing*j, spacing*k);
                positions.push_back(pos);
                positions.push_back(pos+offset1);
                positions.push_back(pos+offset2);
                bonds->addBond(firstParticle, firstParticle+1, dOH, 443153.38);
                bonds->addBond(firstParticle, firstParticle+2, dOH, 443153.38);
                angles->addAngle(firstParticle+1, firstParticle, firstParticle+2, angle, 317.5656);
                nonbonded->addException(firstParticle, firstParticle+1, 0, 1, 0);
                nonbonded->addException(firstParticle, firstParticle+2, 0, 1, 0);
                nonbonded->addException(firstParticle+1, firstParticle+2, 0, 1, 0);
            }
        }
    }
    system.addForce(nonbonded);
    system.addForce(bonds);
    system.addForce(angles);
    RPMDMonteCarloBarostat* barostat = new RPMDMonteCarloBarostat(pressure, frequency);
    system.addForce(barostat);

    // Simulate it and see if the density matches the expected value (1 g/mL).

    RPMDIntegrator integrator(numCopies, temp, 1.0, 0.001);
    Context context(system, integrator, platform);
    for (int copy = 0; copy < numCopies; copy++)
        integrator.setPositions(copy, positions);
    integrator.step(3000);
    double volume = 0.0;
    for (int j = 0; j < steps; ++j) {
        Vec3 box[3];
        context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
        volume += box[0][0]*box[1][1]*box[2][2];
        integrator.step(frequency);
    }
    volume /= steps;
    double density = numMolecules*18/(AVOGADRO*volume*1e-21);
    ASSERT_USUALLY_EQUAL_TOL(1.0, density, 0.04);
}

void testAnisotropicWater() {
    const int numCopies = 8;
    const int gridSize = 8;
    const int numMolecules = gridSize*gridSize*gridSize;
    const int frequency = 10;
    const int steps = 400;
    const double temp = 273.15;
    const double pressure = 3;
    const double spacing = 0.31;
    const double angle = 112*M_PI/180;
    const double dOH = 0.1;

    // Create a box of q-SPC/Fw water molecules.

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(gridSize*spacing, 0, 0), Vec3(0, gridSize*spacing, 0), Vec3(0, 0, gridSize*spacing));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setUseDispersionCorrection(true);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    HarmonicAngleForce* angles = new HarmonicAngleForce();
    vector<Vec3> positions;
    Vec3 offset1(dOH, 0, 0);
    Vec3 offset2(dOH*std::cos(angle), dOH*std::sin(angle), 0);
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            for (int k = 0; k < gridSize; ++k) {
                int firstParticle = system.getNumParticles();
                system.addParticle(16.0);
                system.addParticle(1.0);
                system.addParticle(1.0);
                nonbonded->addParticle(-0.84, 0.3165492, 0.650143);
                nonbonded->addParticle(0.42, 1, 0);
                nonbonded->addParticle(0.42, 1, 0);
                Vec3 pos = Vec3(spacing*i, spacing*j, spacing*k);
                positions.push_back(pos);
                positions.push_back(pos+offset1);
                positions.push_back(pos+offset2);
                bonds->addBond(firstParticle, firstParticle+1, dOH, 443153.38);
                bonds->addBond(firstParticle, firstParticle+2, dOH, 443153.38);
                angles->addAngle(firstParticle+1, firstParticle, firstParticle+2, angle, 317.5656);
                nonbonded->addException(firstParticle, firstParticle+1, 0, 1, 0);
                nonbonded->addException(firstParticle, firstParticle+2, 0, 1, 0);
                nonbonded->addException(firstParticle+1, firstParticle+2, 0, 1, 0);
            }
        }
    }
    system.addForce(nonbonded);
    system.addForce(bonds);
    system.addForce(angles);
    RPMDMonteCarloAnisotropicBarostat* barostat = new RPMDMonteCarloAnisotropicBarostat(Vec3(pressure, pressure, pressure), false, true, false, frequency);
    system.addForce(barostat);

    // Simulate it and see if the density matches the expected value (1 g/mL).

    RPMDIntegrator integrator(numCopies, temp, 1.0, 0.001);
    Context context(system, integrator, platform);
    for (int copy = 0; copy < numCopies; copy++)
        integrator.setPositions(copy, positions);
    integrator.step(3000);
    double volume = 0.0;
    for (int j = 0; j < steps; ++j) {
        Vec3 box[3];
        context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
        volume += box[0][0]*box[1][1]*box[2][2];
        integrator.step(frequency);
        ASSERT_EQUAL(gridSize*spacing, box[0][0]);
        ASSERT_EQUAL(gridSize*spacing, box[2][2]);
    }
    volume /= steps;
    double density = numMolecules*18/(AVOGADRO*volume*1e-21);
    ASSERT_USUALLY_EQUAL_TOL(1.0, density, 0.04);
}

void testFlexibleWater() {
    const int numCopies = 8;
    const int gridSize = 8;
    const int numMolecules = gridSize*gridSize*gridSize;
    const int frequency = 10;
    const int steps = 400;
    const double temp = 273.15;
    const double pressure = 3;
    const double spacing = 0.31;
    const double angle = 112*M_PI/180;
    const double dOH = 0.1;

    // Create a box of q-SPC/Fw water molecules.

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(gridSize*spacing, 0, 0), Vec3(0, gridSize*spacing, 0), Vec3(0, 0, gridSize*spacing));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setUseDispersionCorrection(true);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    HarmonicAngleForce* angles = new HarmonicAngleForce();
    vector<Vec3> positions;
    Vec3 offset1(dOH, 0, 0);
    Vec3 offset2(dOH*std::cos(angle), dOH*std::sin(angle), 0);
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            for (int k = 0; k < gridSize; ++k) {
                int firstParticle = system.getNumParticles();
                system.addParticle(16.0);
                system.addParticle(1.0);
                system.addParticle(1.0);
                nonbonded->addParticle(-0.84, 0.3165492, 0.650143);
                nonbonded->addParticle(0.42, 1, 0);
                nonbonded->addParticle(0.42, 1, 0);
                Vec3 pos = Vec3(spacing*i, spacing*j, spacing*k);
                positions.push_back(pos);
                positions.push_back(pos+offset1);
                positions.push_back(pos+offset2);
                bonds->addBond(firstParticle, firstParticle+1, dOH, 443153.38);
                bonds->addBond(firstParticle, firstParticle+2, dOH, 443153.38);
                angles->addAngle(firstParticle+1, firstParticle, firstParticle+2, angle, 317.5656);
                nonbonded->addException(firstParticle, firstParticle+1, 0, 1, 0);
                nonbonded->addException(firstParticle, firstParticle+2, 0, 1, 0);
                nonbonded->addException(firstParticle+1, firstParticle+2, 0, 1, 0);
            }
        }
    }
    system.addForce(nonbonded);
    system.addForce(bonds);
    system.addForce(angles);
    RPMDMonteCarloFlexibleBarostat* barostat = new RPMDMonteCarloFlexibleBarostat(pressure, frequency);
    system.addForce(barostat);

    // Simulate it and see if the density matches the expected value (1 g/mL).

    RPMDIntegrator integrator(numCopies, temp, 1.0, 0.001);
    Context context(system, integrator, platform);
    for (int copy = 0; copy < numCopies; copy++)
        integrator.setPositions(copy, positions);
    integrator.step(4000);
    double volume = 0.0;
    for (int j = 0; j < steps; ++j) {
        Vec3 box[3];
        context.getState(0).getPeriodicBoxVectors(box[0], box[1], box[2]);
        volume += box[0][0]*box[1][1]*box[2][2];
        integrator.step(frequency);
    }
    volume /= steps;
    double density = numMolecules*18/(AVOGADRO*volume*1e-21);
    ASSERT_USUALLY_EQUAL_TOL(1.0, density, 0.07);
}

void setupKernels(int argc, char* argv[]);
void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        setupKernels(argc, argv);
        testIdealGas();
        testAnisotropicIdealGas();
        testIdealGasAxis(0);
        testIdealGasAxis(1);
        testIdealGasAxis(2);
        testFlexibleIdealGas();
        testMembraneIdealGas(RPMDMonteCarloMembraneBarostat::XYIsotropic, RPMDMonteCarloMembraneBarostat::ZFree);
        testMembraneIdealGas(RPMDMonteCarloMembraneBarostat::XYIsotropic, RPMDMonteCarloMembraneBarostat::ZFixed);
        testMembraneIdealGas(RPMDMonteCarloMembraneBarostat::XYIsotropic, RPMDMonteCarloMembraneBarostat::ConstantVolume);
        testMembraneIdealGas(RPMDMonteCarloMembraneBarostat::XYAnisotropic, RPMDMonteCarloMembraneBarostat::ZFree);
        testMembraneIdealGas(RPMDMonteCarloMembraneBarostat::XYAnisotropic, RPMDMonteCarloMembraneBarostat::ZFixed);
        testMembraneIdealGas(RPMDMonteCarloMembraneBarostat::XYAnisotropic, RPMDMonteCarloMembraneBarostat::ConstantVolume);
        runPlatformTests();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
