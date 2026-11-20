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

void testExclusionTiles() {
    // Verify half-diagonal exclusion tile evaluation, partial atom blocks,
    // cross-block exclusions, and energy parameter derivatives.

    const int numParticles = 100;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(5, 0, 0), Vec3(0, 5, 0), Vec3(0, 0, 5));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(1.5);
    CustomNonbondedForce* custom = new CustomNonbondedForce("scale*(a1+a2)/(r^2+0.2)");
    custom->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    custom->setCutoffDistance(1.5);
    custom->addPerParticleParameter("a");
    custom->addGlobalParameter("scale", 1.5);
    custom->addEnergyParameterDerivative("scale");
    system.addForce(nonbonded);
    system.addForce(custom);
    vector<Vec3> positions(numParticles);
    vector<pair<int, int> > bonds;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(12.0);
        nonbonded->addParticle(i%2 == 0 ? 0.3 : -0.3, 0.2, 0.15);
        custom->addParticle({0.1 + 0.02*(i%5)});
        positions[i] = Vec3(1.0*(i/20) + 0.1*genrand_real2(sfmt), 1.0*((i/4)%5) + 0.1*genrand_real2(sfmt), 1.2*(i%4) + 0.1*genrand_real2(sfmt));
        if (i > 0 && i%3 != 0)
            bonds.push_back({i-1, i});
    }
    nonbonded->createExceptionsFromBonds(bonds, 0.5, 0.5);
    custom->createExclusionsFromBonds(bonds, 3);
    VerletIntegrator cudaIntegrator(0.001), referenceIntegrator(0.001);
    ReferencePlatform referencePlatform;
    Context cudaContext(system, cudaIntegrator, platform);
    Context referenceContext(system, referenceIntegrator, referencePlatform);
    cudaContext.setPositions(positions);
    referenceContext.setPositions(positions);
    int infoMask = State::Forces | State::Energy | State::ParameterDerivatives;
    State cudaState = cudaContext.getState(infoMask);
    State referenceState = referenceContext.getState(infoMask);
    ASSERT_EQUAL_TOL(referenceState.getPotentialEnergy(), cudaState.getPotentialEnergy(), 1e-4);
    ASSERT_EQUAL_TOL(referenceState.getEnergyParameterDerivatives().at("scale"), cudaState.getEnergyParameterDerivatives().at("scale"), 1e-4);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(referenceState.getForces()[i], cudaState.getForces()[i], 1e-4);
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
    testExclusionTiles();
    testParallelComputation(NonbondedForce::NoCutoff);
    testParallelComputation(NonbondedForce::Ewald);
    testParallelComputation(NonbondedForce::PME);
    testParallelComputation(NonbondedForce::LJPME);
    testReordering();
    testDeterministicForces();
    if (canRunHugeTest())
        testHugeSystem();
}
