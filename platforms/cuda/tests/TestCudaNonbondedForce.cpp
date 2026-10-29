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

void testNeighborListRebuildAndOverflow() {
    // Verify neighbor list rebuilds (>3,000 atoms so useNeighborList is active), checkpoint
    // restores, periodic box shrinkage, triclinic culling (including skewed #4343 nearest-image
    // edge cases), neighbor-list capacity overflow growth, and >90,000-atom useLargeBlocks.
    const int numParticles = 3456;
    ReferencePlatform refPlatform;

    for (int boxMode = 0; boxMode < 3; boxMode++) {
        System system;
        if (boxMode == 0)
            system.setDefaultPeriodicBoxVectors(Vec3(12.0, 0, 0), Vec3(0, 12.0, 0), Vec3(0, 0, 12.0));
        else if (boxMode == 1)
            system.setDefaultPeriodicBoxVectors(Vec3(12.0, 0, 0), Vec3(3.5, 11.5, 0), Vec3(-2.5, 3.0, 11.0));
        else
            system.setDefaultPeriodicBoxVectors(Vec3(12.0, 0, 0), Vec3(5.0, 5.2, 0), Vec3(0, 0, 12.0));
        NonbondedForce* nonbonded = new NonbondedForce();
        nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
        nonbonded->setCutoffDistance(1.0);
        system.addForce(nonbonded);
        vector<Vec3> positions(numParticles), velocities(numParticles);
        OpenMM_SFMT::SFMT sfmt;
        init_gen_rand(42 + boxMode, sfmt);
        for (int i = 0; i < numParticles; i++) {
            system.addParticle(10.0);
            nonbonded->addParticle((i%2 == 0 ? 0.2 : -0.2), 0.18, (boxMode == 2 ? 0.01 : 0.1));
            if (i < 32) {
                positions[i] = Vec3(0.35 * (i % 8) * 3.5, (boxMode == 2 ? 0.15 : 0.35) * ((i / 8) % 4) * 7.0, 0.35 * i);
            }
            else if (boxMode == 2) {
                int b = i / 32, lane = i % 32;
                int bx = b % 6, by = (b / 6) % 3, bz = b / 18;
                int lx = lane % 4, ly = (lane / 4) % 4, lz = lane / 16;
                positions[i] = Vec3(2.0 * bx + 0.23 * by + 0.35 * lx, 1.73 * by + 0.56 * ly, 2.0 * bz + 0.23 * by + 0.45 * lz);
            }
            else {
                int ix = i / 288, iy = (i / 24) % 12, iz = i % 24;
                positions[i] = Vec3(0.85 * ix + 0.05 * genrand_real2(sfmt), 0.85 * iy + 0.05 * genrand_real2(sfmt), 0.44 * iz + 0.05 * genrand_real2(sfmt));
            }
            velocities[i] = Vec3(2.5 * (genrand_real2(sfmt) - 0.5), 2.5 * (genrand_real2(sfmt) - 0.5), 2.5 * (genrand_real2(sfmt) - 0.5));
        }
        for (int i = 0; i + 1 < numParticles; i += 3)
            nonbonded->addException(i, i + 1, 0.0, 1.0, 0.0);

        VerletIntegrator integrator(0.002), refIntegrator(0.002);
        Context context(system, integrator, platform), refContext(system, refIntegrator, refPlatform);
        context.setPositions(positions);
        context.setVelocities(velocities);
        refContext.setPositions(positions);
        State initialRef = refContext.getState(State::Forces | State::Energy);
        State initialCuda = context.getState(State::Forces | State::Energy);
        ASSERT_EQUAL_TOL(initialRef.getPotentialEnergy(), initialCuda.getPotentialEnergy(), 1e-3);
        for (int i = 0; i < numParticles; i++)
            ASSERT_EQUAL_VEC(initialRef.getForces()[i], initialCuda.getForces()[i], 1e-3);

        stringstream checkpoint;
        context.createCheckpoint(checkpoint);

        integrator.step(35);
        State cudaState = context.getState(State::Positions | State::Forces | State::Energy);
        refContext.setPositions(cudaState.getPositions());
        State refState = refContext.getState(State::Forces | State::Energy);
        ASSERT_EQUAL_TOL(refState.getPotentialEnergy(), cudaState.getPotentialEnergy(), 1e-3);
        for (int i = 0; i < numParticles; i++)
            ASSERT_EQUAL_VEC(refState.getForces()[i], cudaState.getForces()[i], 1e-3);

        if (boxMode == 0) {
            context.loadCheckpoint(checkpoint);
            State restoredState = context.getState(State::Forces | State::Energy);
            ASSERT_EQUAL_TOL(initialRef.getPotentialEnergy(), restoredState.getPotentialEnergy(), 1e-3);
            for (int i = 0; i < numParticles; i++)
                ASSERT_EQUAL_VEC(initialRef.getForces()[i], restoredState.getForces()[i], 1e-3);

            context.setPeriodicBoxVectors(Vec3(11.4, 0, 0), Vec3(0, 11.4, 0), Vec3(0, 0, 11.4));
            refContext.setPositions(positions);
            refContext.setPeriodicBoxVectors(Vec3(11.4, 0, 0), Vec3(0, 11.4, 0), Vec3(0, 0, 11.4));
            State shrunkCuda = context.getState(State::Forces | State::Energy);
            State shrunkRef = refContext.getState(State::Forces | State::Energy);
            ASSERT_EQUAL_TOL(shrunkRef.getPotentialEnergy(), shrunkCuda.getPotentialEnergy(), 1e-3);
            for (int i = 0; i < numParticles; i++)
                ASSERT_EQUAL_VEC(shrunkRef.getForces()[i], shrunkCuda.getForces()[i], 1e-3);
        }
    }

    // Exercise neighbor-list overflow and resize (~4,950 interacting tiles > initial maxTiles = 2,000).
    const int denseParticles = 3200;
    System denseSystem;
    denseSystem.setDefaultPeriodicBoxVectors(Vec3(6.0, 0, 0), Vec3(0, 6.0, 0), Vec3(0, 0, 6.0));
    NonbondedForce* denseNb = new NonbondedForce();
    denseNb->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    denseNb->setCutoffDistance(1.5);
    denseSystem.addForce(denseNb);
    vector<Vec3> densePos(denseParticles);
    for (int i = 0; i < denseParticles; i++) {
        denseSystem.addParticle(10.0);
        denseNb->addParticle((i%2 == 0 ? 0.05 : -0.05), 0.05, 1e-4);
        densePos[i] = Vec3(0.08 * (i % 16), 0.08 * ((i / 16) % 16), 0.08 * (i / 256));
    }
    VerletIntegrator denseInt(0.001), denseRefInt(0.001);
    Context denseCtx(denseSystem, denseInt, platform), denseRefCtx(denseSystem, denseRefInt, refPlatform);
    denseCtx.setPositions(densePos);
    denseRefCtx.setPositions(densePos);
    State denseCuda = denseCtx.getState(State::Forces | State::Energy);
    State denseRef = denseRefCtx.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(denseRef.getPotentialEnergy(), denseCuda.getPotentialEnergy(), 1e-3);
    for (int i = 0; i < denseParticles; i++)
        ASSERT_EQUAL_VEC(denseRef.getForces()[i], denseCuda.getForces()[i], 1e-3);

    // Exercise >90,000-atom useLargeBlocks 32-block superblocks, dispersed largeBlocksGlobal
    // blocks, non-rebuild skips, periodic box shrink, and #4343 triclinic superblock overrides
    // against ReferencePlatform via exact 36-fold periodic supercell replication (110,592 atoms).
    const int baseParticles = 3072;
    for (int triclinicMode = 0; triclinicMode < 2; triclinicMode++) {
        const int rx = (triclinicMode == 0 ? 3 : 6);
        const int ry = (triclinicMode == 0 ? 3 : 6);
        const int rz = (triclinicMode == 0 ? 4 : 1);
        const int numReplicas = rx * ry * rz;
        const int superParticles = numReplicas * baseParticles;
        Vec3 baseA = (triclinicMode == 0 ? Vec3(8.0, 0, 0) : Vec3(6.0, 0, 0));
        Vec3 baseB = (triclinicMode == 0 ? Vec3(2.0, 7.8, 0) : Vec3(1.5, 5.8, 0));
        Vec3 baseC = (triclinicMode == 0 ? Vec3(-1.5, 1.8, 7.6) : Vec3(0, 2.2, 5.2));

        System baseSystem, superSystem;
        baseSystem.setDefaultPeriodicBoxVectors(baseA, baseB, baseC);
        superSystem.setDefaultPeriodicBoxVectors(baseA * rx, baseB * ry, baseC * rz);
        NonbondedForce* baseNb = new NonbondedForce();
        NonbondedForce* superNb = new NonbondedForce();
        baseNb->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
        superNb->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
        baseNb->setCutoffDistance(1.0);
        superNb->setCutoffDistance(1.0);
        baseNb->setReactionFieldDielectric(1e6);
        superNb->setReactionFieldDielectric(1e6);
        baseNb->setUseSwitchingFunction(true);
        superNb->setUseSwitchingFunction(true);
        baseNb->setSwitchingDistance(0.8);
        superNb->setSwitchingDistance(0.8);
        baseNb->setUseDispersionCorrection(false);
        superNb->setUseDispersionCorrection(false);
        baseSystem.addForce(baseNb);
        superSystem.addForce(superNb);

        vector<Vec3> basePos(baseParticles), baseVel(baseParticles);
        OpenMM_SFMT::SFMT sfmt;
        init_gen_rand(100 + triclinicMode, sfmt);
        for (int i = 0; i < baseParticles; i++) {
            baseSystem.addParticle(10.0);
            baseNb->addParticle((i % 3 == 0 ? 0.08 : (i % 3 == 1 ? -0.05 : -0.03)), 0.10, 0.02);
            int b = i / 32, lane = i % 32;
            int bx = b % 6, by = (b / 6) % 4, bz = b / 24;
            int lx = lane % 4, ly = (lane / 4) % 4, lz = lane / 16;
            double zSpan = (triclinicMode == 0 ? 0.45 * lz : 1.68 * lz);
            double u0 = genrand_real2(sfmt), u1 = genrand_real2(sfmt), u2 = genrand_real2(sfmt);
            double v0 = genrand_real2(sfmt), v1 = genrand_real2(sfmt), v2 = genrand_real2(sfmt);
            Vec3 frac((bx + 0.24 * lx + 0.03 * u0) / 6.0,
                      (by + 0.24 * ly + 0.03 * u1) / 4.0,
                      (bz + 0.03 * u2) / 4.0);
            basePos[i] = baseA * frac[0] + baseB * frac[1] + baseC * frac[2] + Vec3(0, 0, zSpan);
            baseVel[i] = Vec3(1.8 * (v0 - 0.5), 1.8 * (v1 - 0.5), 1.8 * (v2 - 0.5));
        }
        for (int i = 0; i + 1 < baseParticles; i += 4)
            baseNb->addException(i, i + 1, 0.0, 1.0, 0.0);

        vector<int> superToOrig(superParticles), origToSuper(superParticles);
        for (int s = 0; s < superParticles; s++)
            superToOrig[s] = s;
        for (int r = 0; r < 32; r++)
            swap(superToOrig[r], superToOrig[r * baseParticles]);
        for (int s = 0; s < superParticles; s++)
            origToSuper[superToOrig[s]] = s;

        vector<Vec3> superPos(superParticles), superVel(superParticles);
        for (int s = 0; s < superParticles; s++) {
            int orig = superToOrig[s];
            int rep = orig / baseParticles, baseIdx = orig % baseParticles;
            int ix = rep % rx, iy = (rep / rx) % ry, iz = rep / (rx * ry);
            superSystem.addParticle(10.0);
            double q, sig, eps;
            baseNb->getParticleParameters(baseIdx, q, sig, eps);
            superNb->addParticle(q, sig, eps);
            superPos[s] = basePos[baseIdx] + baseA * ix + baseB * iy + baseC * iz;
            superVel[s] = baseVel[baseIdx];
        }
        for (int r = 0; r < numReplicas; r++)
            for (int i = 0; i + 1 < baseParticles; i += 4)
                superNb->addException(origToSuper[r * baseParticles + i], origToSuper[r * baseParticles + i + 1], 0.0, 1.0, 0.0);

        VerletIntegrator baseRefInt(0.002), baseCudaInt(0.002), superInt(0.002);
        Context baseRefCtx(baseSystem, baseRefInt, refPlatform);
        Context baseCudaCtx(baseSystem, baseCudaInt, platform);
        Context superCtx(superSystem, superInt, platform);
        baseRefCtx.setPositions(basePos);
        baseCudaCtx.setPositions(basePos);
        baseCudaCtx.setVelocities(baseVel);
        superCtx.setPositions(superPos);
        superCtx.setVelocities(superVel);

        State initRef = baseRefCtx.getState(State::Forces | State::Energy);
        State initSuper = superCtx.getState(State::Forces | State::Energy);
        ASSERT_EQUAL_TOL(numReplicas * initRef.getPotentialEnergy(), initSuper.getPotentialEnergy(), 1e-2);
        for (int s = 0; s < superParticles; s++)
            ASSERT_EQUAL_VEC(initRef.getForces()[superToOrig[s] % baseParticles], initSuper.getForces()[s], 1e-3);

        baseCudaInt.step(35);
        superInt.step(35);
        State steppedBase = baseCudaCtx.getState(State::Positions);
        baseRefCtx.setPositions(steppedBase.getPositions());
        State steppedRef = baseRefCtx.getState(State::Forces | State::Energy);
        State steppedSuper = superCtx.getState(State::Positions | State::Forces | State::Energy);
        ASSERT_EQUAL_TOL(numReplicas * steppedRef.getPotentialEnergy(), steppedSuper.getPotentialEnergy(), 2e-2);
        VerletIntegrator verifyStepInt(0.002);
        Context verifyStepCtx(superSystem, verifyStepInt, platform);
        verifyStepCtx.setPositions(steppedSuper.getPositions());
        State freshStepSuper = verifyStepCtx.getState(State::Forces | State::Energy);
        ASSERT_EQUAL_TOL(freshStepSuper.getPotentialEnergy(), steppedSuper.getPotentialEnergy(), 1e-2);
        for (int s = 0; s < superParticles; s++)
            ASSERT_EQUAL_VEC(freshStepSuper.getForces()[s], steppedSuper.getForces()[s], 1e-3);

        if (triclinicMode == 0) {
            Vec3 shrunkA = baseA * (0.992 * rx), shrunkB = baseB * (0.992 * ry), shrunkC = baseC * (0.992 * rz);
            superCtx.setPeriodicBoxVectors(shrunkA, shrunkB, shrunkC);
            State shrunkSuper = superCtx.getState(State::Forces | State::Energy);
            VerletIntegrator freshInt(0.002);
            Context freshSuperCtx(superSystem, freshInt, platform);
            freshSuperCtx.setPeriodicBoxVectors(shrunkA, shrunkB, shrunkC);
            freshSuperCtx.setPositions(steppedSuper.getPositions());
            State freshSuper = freshSuperCtx.getState(State::Forces | State::Energy);
            ASSERT_EQUAL_TOL(freshSuper.getPotentialEnergy(), shrunkSuper.getPotentialEnergy(), 1e-2);
            for (int s = 0; s < superParticles; s++)
                ASSERT_EQUAL_VEC(freshSuper.getForces()[s], shrunkSuper.getForces()[s], 1e-3);
        }
    }

    // Exercise >90,000-atom CutoffNonPeriodic when a dispersed block's half-width exceeds FP16 max (65,504 nm -> +Inf).
    const int npParticles = 90050;
    System npSystem, npRefSystem;
    NonbondedForce* npNb = new NonbondedForce();
    NonbondedForce* npRefNb = new NonbondedForce();
    npNb->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
    npRefNb->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
    npNb->setCutoffDistance(1.0);
    npRefNb->setCutoffDistance(1.0);
    npSystem.addForce(npNb);
    npRefSystem.addForce(npRefNb);
    vector<Vec3> npPos(npParticles);
    for (int i = 0; i < npParticles; i++) {
        npSystem.addParticle(10.0);
        if (i == 0) {
            npNb->addParticle(0.51, 0.2, 0.0);
            npPos[i] = Vec3(-70000.0, 0, 0);
        }
        else if (i == 1) {
            npNb->addParticle(0.49, 0.2, 0.1);
            npPos[i] = Vec3(70000.0, 0, 0);
        }
        else if (i == 32) {
            npNb->addParticle(-0.50, 0.2, 0.1);
            npPos[i] = Vec3(70000.375, 0, 0);
        }
        else {
            npNb->addParticle(0.0, 0.1, 0.0);
            npPos[i] = Vec3(3.0 * (i % 64), 3.0 * ((i / 64) % 64), 3.0 * (i / 4096) + 10.0);
        }
    }
    npNb->addException(0, 1, 0.0, 1.0, 0.0);
    const int refIndices[3] = {0, 1, 32};
    vector<Vec3> npRefPos(3);
    for (int k = 0; k < 3; k++) {
        npRefSystem.addParticle(10.0);
        double q, sig, eps;
        npNb->getParticleParameters(refIndices[k], q, sig, eps);
        npRefNb->addParticle(q, sig, eps);
        npRefPos[k] = npPos[refIndices[k]];
    }
    npRefNb->addException(0, 1, 0.0, 1.0, 0.0);
    VerletIntegrator npInt(0.001), npRefInt(0.001);
    Context npCtx(npSystem, npInt, platform), npRefCtx(npRefSystem, npRefInt, refPlatform);
    npCtx.setPositions(npPos);
    npRefCtx.setPositions(npRefPos);
    State npCuda = npCtx.getState(State::Forces | State::Energy);
    State npRef = npRefCtx.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(npRef.getPotentialEnergy(), npCuda.getPotentialEnergy(), 1e-4);
    for (int k = 0; k < 3; k++)
        ASSERT_EQUAL_VEC(npRef.getForces()[k], npCuda.getForces()[refIndices[k]], 1e-4);
}

void runPlatformTests() {
    testNeighborListRebuildAndOverflow();
    testParallelComputation(NonbondedForce::NoCutoff);
    testParallelComputation(NonbondedForce::Ewald);
    testParallelComputation(NonbondedForce::PME);
    testParallelComputation(NonbondedForce::LJPME);
    testReordering();
    testDeterministicForces();
    if (canRunHugeTest())
        testHugeSystem();
}
