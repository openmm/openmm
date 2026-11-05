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
#include "openmm/CustomNonbondedForce.h"
#include <cuda.h>
#include <string>

void testSubtileWarpEvaluation() {
    // Exercise >3,000-atom neighbor-list evaluation, 4/8/16/32 sub-tile compaction, dynamic
    // work-stealing across multiple MD steps, half-diagonal exclusion tiles, single pairs,
    // energy-only kernels, and high-parameter CustomNonbondedForce launch bounds.
    const int numParticles = 3456;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(12.0, 0, 0), Vec3(0, 12.0, 0), Vec3(0, 0, 12.0));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(1.0);
    CustomNonbondedForce* custom = new CustomNonbondedForce("(a1+a2)*(b1+b2)*(c1+c2)*(d1+d2)*(e1+e2)/(r^2+0.2)");
    custom->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    custom->setCutoffDistance(1.0);
    for (const char* name : {"a", "b", "c", "d", "e"})
        custom->addPerParticleParameter(name);
    system.addForce(nonbonded);
    system.addForce(custom);
    vector<Vec3> positions(numParticles), velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(42, sfmt);
    vector<pair<int, int> > bonds;
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(12.0);
        nonbonded->addParticle((i%2 == 0 ? 0.25 : -0.25), 0.20, 0.12);
        custom->addParticle({0.1 + 0.01*(i%5), 0.2, 0.15, 0.05*(1 + i%3), 0.1});
        int ix = i / 288, iy = (i / 24) % 12, iz = i % 24;
        positions[i] = Vec3(0.95*ix + 0.2*genrand_real2(sfmt), 0.95*iy + 0.2*genrand_real2(sfmt), 0.48*iz + 0.2*genrand_real2(sfmt));
        velocities[i] = Vec3(0.3*(genrand_real2(sfmt)-0.5), 0.3*(genrand_real2(sfmt)-0.5), 0.3*(genrand_real2(sfmt)-0.5));
        if (i%4 != 0) {
            bonds.push_back({i-1, i});
            for (int d = 1; d <= i%4; d++)
                custom->addExclusion(i-d, i);
        }
    }
    nonbonded->createExceptionsFromBonds(bonds, 0.5, 0.5);
    VerletIntegrator intCuda(0.002), intRef(0.002);
    ReferencePlatform refPlatform;
    Context ctxCuda(system, intCuda, platform), ctxRef(system, intRef, refPlatform);
    ctxCuda.setPositions(positions);
    ctxCuda.setVelocities(velocities);
    ctxRef.setPositions(positions);
    State fullRef = ctxRef.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(fullRef.getPotentialEnergy(), ctxCuda.getState(State::Energy).getPotentialEnergy(), 1e-3);
    State fullCuda = ctxCuda.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(fullRef.getPotentialEnergy(), fullCuda.getPotentialEnergy(), 1e-3);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(fullRef.getForces()[i], fullCuda.getForces()[i], 1e-3);

    intCuda.step(5);
    State stepCuda = ctxCuda.getState(State::Positions | State::Forces | State::Energy);
    ctxRef.setPositions(stepCuda.getPositions());
    State stepRef = ctxRef.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(stepRef.getPotentialEnergy(), stepCuda.getPotentialEnergy(), 1e-3);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(stepRef.getForces()[i], stepCuda.getForces()[i], 1e-3);
}

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
    testSubtileWarpEvaluation();
    if (canRunHugeTest())
        testHugeSystem();
}
