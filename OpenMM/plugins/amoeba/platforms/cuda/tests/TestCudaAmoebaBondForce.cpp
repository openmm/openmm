/* -------------------------------------------------------------------------- *
 *                                   OpenMMAmoeba                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs                                                   *
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
 * This tests the Cuda implementation of AmoebaBondForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "CudaPlatform.h"
#include "openmm/Context.h"
#include "openmm/CustomBondForce.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;

extern "C" void registerAmoebaCudaKernelFactories();

const double TOL = 1e-5;
static void computeAmoebaBondForce(int bondIndex,  std::vector<Vec3>& positions, AmoebaBondForce& amoebaBondForce,
                                           std::vector<Vec3>& forces, double* energy) {

    int particle1, particle2;
    double bondLength;
    double quadraticK;
    double cubicK    = amoebaBondForce.getAmoebaGlobalBondCubic();
    double quarticK  = amoebaBondForce.getAmoebaGlobalBondQuartic();
    amoebaBondForce.getBondParameters(bondIndex, particle1, particle2,  bondLength,  quadraticK);

    double deltaR[3];
    double r2 = 0.0;
    for (int ii = 0; ii < 3; ii++) {
           deltaR[ii]    = positions[particle2][ii] - positions[particle1][ii];
           r2           += deltaR[ii]*deltaR[ii];
    }
    double r                   = sqrt(r2);

    double bondDelta           = (r - bondLength);
    double bondDelta2          = bondDelta*bondDelta;
    double dEdR                = 1.0 + 1.5*cubicK*bondDelta + 2.0*quarticK*bondDelta2;

           dEdR               *= (r > 0.0) ? (2.0*quadraticK*bondDelta)/r : 0.0;

   forces[particle1][0]       += dEdR*deltaR[0];
   forces[particle1][1]       += dEdR*deltaR[1];
   forces[particle1][2]       += dEdR*deltaR[2];

   forces[particle2][0]       -= dEdR*deltaR[0];
   forces[particle2][1]       -= dEdR*deltaR[1];
   forces[particle2][2]       -= dEdR*deltaR[2];

   *energy                    += (1.0f + cubicK*bondDelta + quarticK*bondDelta2)*quadraticK*bondDelta2;

}

static void computeAmoebaBondForces(Context& context, AmoebaBondForce& amoebaBondForce,
                                             std::vector<Vec3>& expectedForces, double* expectedEnergy) {

    // get positions and zero forces

    State state = context.getState(State::Positions);
    std::vector<Vec3> positions = state.getPositions();
    expectedForces.resize(positions.size());
    
    for (unsigned int ii = 0; ii < expectedForces.size(); ii++) {
        expectedForces[ii][0] = expectedForces[ii][1] = expectedForces[ii][2] = 0.0;
    }

    // calculates forces/energy

    *expectedEnergy = 0.0;
    for (int ii = 0; ii < amoebaBondForce.getNumBonds(); ii++) {
        computeAmoebaBondForce(ii, positions, amoebaBondForce, expectedForces, expectedEnergy);
    }
}

void compareWithExpectedForceAndEnergy(Context& context, AmoebaBondForce& amoebaBondForce, double tolerance, const std::string& idString) {

    std::vector<Vec3> expectedForces;
    double expectedEnergy;
    computeAmoebaBondForces(context, amoebaBondForce, expectedForces, &expectedEnergy);
   
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();

    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        ASSERT_EQUAL_VEC(expectedForces[ii], forces[ii], tolerance);
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), tolerance);
}

void testOneBond() {

    System system;

    system.addParticle(1.0);
    system.addParticle(1.0);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaBondForce* amoebaBondForce = new AmoebaBondForce();

    double bondLength = 1.5;
    double quadraticK = 1.0;
    double cubicK     = 2.0;
    double quarticicK = 3.0;
    amoebaBondForce->setAmoebaGlobalBondCubic(cubicK);
    amoebaBondForce->setAmoebaGlobalBondQuartic(quarticicK);
    amoebaBondForce->addBond(0, 1, bondLength, quadraticK);

    system.addForce(amoebaBondForce);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    std::vector<Vec3> positions(2);

    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy(context, *amoebaBondForce, TOL, "testOneBond");
}

void testTwoBond() {

    System system;

    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaBondForce* amoebaBondForce = new AmoebaBondForce();

    double bondLength = 1.5;
    double quadraticK = 1.0;
    double cubicK     = 2.0;
    double quarticicK = 3.0;
    amoebaBondForce->setAmoebaGlobalBondCubic(cubicK);
    amoebaBondForce->setAmoebaGlobalBondQuartic(quarticicK);
    amoebaBondForce->addBond(0, 1, bondLength, quadraticK);
    amoebaBondForce->addBond(1, 2, bondLength, quadraticK);

    system.addForce(amoebaBondForce);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    //Context context(system, integrator, platform);
    std::vector<Vec3> positions(3);

    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 1);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy(context, *amoebaBondForce, TOL, "testTwoBond");
    
    // Try changing the bond parameters and make sure it's still correct.
    
    amoebaBondForce->setBondParameters(0, 0, 1, 1.1*bondLength, 1.4*quadraticK);
    amoebaBondForce->setBondParameters(1, 1, 2, 1.2*bondLength, 0.9*quadraticK);
    bool exceptionThrown = false;
    try {
        // This should throw an exception.
        compareWithExpectedForceAndEnergy(context, *amoebaBondForce, TOL, "testTwoBond");
    }
    catch (std::exception ex) {
        exceptionThrown = true;
    }
    ASSERT(exceptionThrown);
    amoebaBondForce->updateParametersInContext(context);
    compareWithExpectedForceAndEnergy(context, *amoebaBondForce, TOL, "testTwoBond");
}

void testPeriodic() {
    // Create a force that uses periodic boundary conditions, then compare to an identical custom force.
    
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    int numParticles = 2;
    for (int ii = 0; ii < numParticles; ii++)
        system.addParticle(1.0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    AmoebaBondForce* amoebaBondForce = new AmoebaBondForce();
    double bondLength = 1.5;
    double quadraticK = 1.0;
    double cubicK     = 2.0;
    double quarticK   = 3.0;
    amoebaBondForce->setAmoebaGlobalBondCubic(cubicK);
    amoebaBondForce->setAmoebaGlobalBondQuartic(quarticK);
    amoebaBondForce->addBond(0, 1, bondLength, quadraticK);
    amoebaBondForce->setUsesPeriodicBoundaryConditions(true);
    system.addForce(amoebaBondForce);
    CustomBondForce* customForce = new CustomBondForce("k2*delta^2 + k3*delta^3 + k4*delta^4; delta=r-r0");
    customForce->addGlobalParameter("r0", bondLength);
    customForce->addGlobalParameter("k2", quadraticK);
    customForce->addGlobalParameter("k3", cubicK);
    customForce->addGlobalParameter("k4", quarticK);
    customForce->addBond(0, 1);
    customForce->setUsesPeriodicBoundaryConditions(true);
    customForce->setForceGroup(1);
    system.addForce(customForce);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));

    std::vector<Vec3> positions(numParticles);

    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 0);

    context.setPositions(positions);
    State s1 = context.getState(State::Forces | State::Energy, true, 1);
    State s2 = context.getState(State::Forces | State::Energy, true, 2);
    ASSERT_EQUAL_TOL(s2.getPotentialEnergy(), s1.getPotentialEnergy(), 1e-5);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(s2.getForces()[i], s1.getForces()[i], 1e-5);
}

int main(int argc, char* argv[]) {
    try {
        std::cout << "TestCudaAmoebaBondForce running test..." << std::endl;
        registerAmoebaCudaKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("Precision", std::string(argv[1]));
        testTwoBond();
        testPeriodic();
    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
