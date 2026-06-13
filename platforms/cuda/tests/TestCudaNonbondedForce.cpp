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

void getPmeGridForKernel(NonbondedForce::ReciprocalSpaceKernelType kernel, int& nx, int& ny, int& nz) {
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(4, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 4));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    nonbonded->setReciprocalSpaceKernelType(kernel);
    nonbonded->setCutoffDistance(1.0);
    nonbonded->setEwaldErrorTolerance(1e-4);
    system.addForce(nonbonded);
    vector<Vec3> positions;
    for (int i = 0; i < 4; i++) {
        system.addParticle(1.0);
        nonbonded->addParticle(i%2 == 0 ? 1.0 : -1.0, 1.0, 0.0);
        positions.push_back(Vec3(0.3+0.4*i, 0.2+0.3*i, 0.1+0.2*i));
    }
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    double alpha;
    nonbonded->getPMEParametersInContext(context, alpha, nx, ny, nz);
}

void testEspSelectsDistinctPmeGrid() {
    int pmeNx, pmeNy, pmeNz, espNx, espNy, espNz;
    getPmeGridForKernel(NonbondedForce::PMEKernel, pmeNx, pmeNy, pmeNz);
    getPmeGridForKernel(NonbondedForce::ESPKernel, espNx, espNy, espNz);
    ASSERT(pmeNx != espNx || pmeNy != espNy || pmeNz != espNz);
}

void testEspRequiresPme() {
    System system;
    system.addParticle(1.0);
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setReciprocalSpaceKernelType(NonbondedForce::ESPKernel);
    nonbonded->addParticle(1.0, 1.0, 0.0);
    system.addForce(nonbonded);
    VerletIntegrator integrator(0.001);
    bool threw = false;
    try {
        Context context(system, integrator, platform);
    }
    catch (const OpenMMException& e) {
        ASSERT(string(e.what()).find("only supported for PME") != string::npos);
        threw = true;
    }
    ASSERT(threw);
}

void testEspDirectSpaceUsesPswfSplit() {
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(2, 0, 0), Vec3(0, 2, 0), Vec3(0, 0, 2));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    nonbonded->setReciprocalSpaceKernelType(NonbondedForce::ESPKernel);
    nonbonded->setReciprocalSpaceForceGroup(1);
    nonbonded->setCutoffDistance(0.5);
    nonbonded->setEwaldErrorTolerance(1e-4);
    system.addForce(nonbonded);
    system.addParticle(1.0);
    system.addParticle(1.0);
    nonbonded->addParticle(1.0, 1.0, 0.0);
    nonbonded->addParticle(1.0, 1.0, 0.0);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions({Vec3(0, 0, 0), Vec3(0.4, 0, 0)});

    const double directEnergy = context.getState(State::Energy, true, 1<<0).getPotentialEnergy();
    const double espShortRangeScale = 1.847619872734e-3;
    ASSERT_EQUAL_TOL(ONE_4PI_EPS0*espShortRangeScale/0.4, directEnergy, 1e-5);
}

void computePmeKernelForces(NonbondedForce::ReciprocalSpaceKernelType kernel, bool useChargeOffsets,
        double& energy, vector<Vec3>& forces) {
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    nonbonded->setReciprocalSpaceKernelType(kernel);
    nonbonded->setCutoffDistance(0.9);
    nonbonded->setEwaldErrorTolerance(1e-4);
    system.addForce(nonbonded);
    vector<double> charges = {1.0, -1.0, 0.7, -0.7, 0.4, -0.4, 0.2, -0.2};
    for (int i = 0; i < charges.size(); i++) {
        system.addParticle(1.0);
        nonbonded->addParticle(charges[i], 1.0, 0.0);
    }
    nonbonded->addException(0, 1, 0.0, 1.0, 0.0);
    if (useChargeOffsets) {
        nonbonded->addGlobalParameter("scale", 0.3);
        nonbonded->addParticleParameterOffset("scale", 0, 0.5, 0.0, 0.0);
        nonbonded->addParticleParameterOffset("scale", 1, -0.4, 0.0, 0.0);
        nonbonded->addParticleParameterOffset("scale", 2, 0.2, 0.0, 0.0);
    }
    vector<Vec3> positions = {
        Vec3(0.1, 0.2, 0.3),
        Vec3(0.8, 0.4, 0.2),
        Vec3(1.4, 0.9, 0.7),
        Vec3(2.1, 0.6, 1.2),
        Vec3(0.5, 1.8, 1.6),
        Vec3(1.7, 2.2, 0.4),
        Vec3(2.5, 1.5, 2.0),
        Vec3(1.2, 2.7, 2.4)
    };
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    if (useChargeOffsets)
        context.setParameter("scale", 0.7);
    State state = context.getState(State::Forces | State::Energy);
    energy = state.getPotentialEnergy();
    forces = state.getForces();
}

void testEspMatchesPmeCoulomb() {
    double pmeEnergy, espEnergy;
    vector<Vec3> pmeForces, espForces;
    computePmeKernelForces(NonbondedForce::PMEKernel, false, pmeEnergy, pmeForces);
    computePmeKernelForces(NonbondedForce::ESPKernel, false, espEnergy, espForces);
    ASSERT_EQUAL_TOL(pmeEnergy, espEnergy, 5e-3);
    for (int i = 0; i < pmeForces.size(); i++)
        ASSERT_EQUAL_VEC(pmeForces[i], espForces[i], 5e-3);
}

void testEspMatchesPmeCoulombWithChargeOffsets() {
    double pmeEnergy, espEnergy;
    vector<Vec3> pmeForces, espForces;
    computePmeKernelForces(NonbondedForce::PMEKernel, true, pmeEnergy, pmeForces);
    computePmeKernelForces(NonbondedForce::ESPKernel, true, espEnergy, espForces);
    ASSERT_EQUAL_TOL(pmeEnergy, espEnergy, 5e-3);
    for (int i = 0; i < pmeForces.size(); i++)
        ASSERT_EQUAL_VEC(pmeForces[i], espForces[i], 5e-3);
}

void runPlatformTests() {
    testParallelComputation(NonbondedForce::NoCutoff);
    testParallelComputation(NonbondedForce::Ewald);
    testParallelComputation(NonbondedForce::PME);
    testParallelComputation(NonbondedForce::LJPME);
    testEspSelectsDistinctPmeGrid();
    testEspRequiresPme();
    testEspDirectSpaceUsesPswfSplit();
    testEspMatchesPmeCoulomb();
    testEspMatchesPmeCoulombWithChargeOffsets();
    testReordering();
    testDeterministicForces();
    if (canRunHugeTest())
        testHugeSystem();
}
