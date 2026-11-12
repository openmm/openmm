/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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

#include "CudaTests.h"
#include "TestNonbondedForce.h"
#include <cuda.h>
#include <string>

void testDeterministicForces() {
    // Check that the CudaDeterministicForces property works correctly.
    
    const int numParticles = 1000;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(6, 0, 0), Vec3(2.1, 6, 0), Vec3(-1.5, -0.5, 6));
    NonbondedForce *nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    system.addForce(nonbonded);
    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        nonbonded->addParticle(i%2 == 0 ? 1 : -1, 1, 0);
        positions.push_back(Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5)*6);
    }
    VerletIntegrator integrator(0.001);
    map<string, string> properties;
    properties[CudaPlatform::CudaDeterministicForces()] = "true";
    Context context(system, integrator, platform, properties);
    context.setPositions(positions);
    State state1 = context.getState(State::Forces);
    State state2 = context.getState(State::Forces);
    
    // All forces should be *exactly* equal.
    
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL(state1.getForces()[i][0], state2.getForces()[i][0]);
        ASSERT_EQUAL(state1.getForces()[i][1], state2.getForces()[i][1]);
        ASSERT_EQUAL(state1.getForces()[i][2], state2.getForces()[i][2]);
    }
}

void testLargeBlockNeighborList() {
    // Create a >90,000-atom triclinic system with swapped endpoints so both compact superblocks and dispersed blocks are exercised.

    const int gridX = 46, gridY = 45, gridZ = 44, numParticles = gridX*gridY*gridZ;
    Vec3 boxX(18.4, 0, 0), boxY(2.0, 18.0, 0), boxZ(-1.5, 1.5, 17.6);
    System system;
    system.setDefaultPeriodicBoxVectors(boxX, boxY, boxZ);
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(0.9);
    nonbonded->setReactionFieldDielectric(1e6);
    nonbonded->setUseSwitchingFunction(true);
    nonbonded->setSwitchingDistance(0.75);
    system.addForce(nonbonded);
    vector<Vec3> positions(numParticles), velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(10.0);
        nonbonded->addParticle(i%2 == 0 ? 0.2 : -0.2, 0.15, 0.2);
        positions[i] = boxX*((i/(gridY*gridZ)+0.1*genrand_real2(sfmt))/gridX) + boxY*(((i/gridZ)%gridY+0.1*genrand_real2(sfmt))/gridY) + boxZ*((i%gridZ+0.1*genrand_real2(sfmt))/gridZ);
        velocities[i] = Vec3(2.0*(genrand_real2(sfmt)-0.5), 2.0*(genrand_real2(sfmt)-0.5), 2.0*(genrand_real2(sfmt)-0.5));
    }
    for (int i = 0; i < 16; i++)
        swap(positions[i], positions[numParticles-1-i]);
    for (int i = 0; i+1 < numParticles; i += 3)
        nonbonded->addException(i, i+1, 0.0, 0.15, 0.0);

    // Compare forces and energies against the Reference platform across reuse steps, box changes, and parameter updates.

    ReferencePlatform reference;
    VerletIntegrator integrator(0.002), referenceIntegrator(0.002);
    Context context(system, integrator, platform), referenceContext(system, referenceIntegrator, reference);
    context.setPositions(positions);
    context.setVelocities(velocities);
    auto verifyAgainstReference = [&]() {
        State cudaState = context.getState(State::Positions | State::Forces | State::Energy);
        referenceContext.setPositions(cudaState.getPositions());
        State referenceState = referenceContext.getState(State::Forces | State::Energy);
        ASSERT_EQUAL_TOL(referenceState.getPotentialEnergy(), cudaState.getPotentialEnergy(), 1e-3);
        for (int i = 0; i < numParticles; i++)
            ASSERT_EQUAL_VEC(referenceState.getForces()[i], cudaState.getForces()[i], 1e-3);
    };
    verifyAgainstReference();
    integrator.step(5);
    verifyAgainstReference();
    context.setPeriodicBoxVectors(boxX*1.01, boxY*1.01, boxZ*1.01);
    referenceContext.setPeriodicBoxVectors(boxX*1.01, boxY*1.01, boxZ*1.01);
    verifyAgainstReference();
    integrator.step(25);
    verifyAgainstReference();
    for (int i = 0; i < numParticles; i += 7)
        nonbonded->setParticleParameters(i, i%2 == 0 ? 0.15 : -0.15, 0.15, 0.25);
    nonbonded->updateParametersInContext(context);
    nonbonded->updateParametersInContext(referenceContext);
    verifyAgainstReference();
}

bool canRunHugeTest() {
    // Create a minimal context just to see which device is being used.

    System system;
    system.addParticle(1.0);
    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    int deviceIndex = stoi(platform.getPropertyValue(context, CudaPlatform::CudaDeviceIndex()));

    // Find out how much memory the device has.

    CUdevice device;
    cuDeviceGet(&device, deviceIndex);
    size_t memory;
    cuDeviceTotalMem(&memory, device);

    // Only run the huge test if the device has at least 4 GB of memory.

    return (memory >= 4L*(1<<30));
}

void runPlatformTests() {
    testParallelComputation(NonbondedForce::NoCutoff);
    testParallelComputation(NonbondedForce::Ewald);
    testParallelComputation(NonbondedForce::PME);
    testParallelComputation(NonbondedForce::LJPME);
    testReordering();
    testDeterministicForces();
    if (canRunHugeTest())
        testHugeSystem();
    testLargeBlockNeighborList();
}
