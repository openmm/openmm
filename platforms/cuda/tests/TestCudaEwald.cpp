/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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
#include "TestEwald.h"
#include <cuda.h>

void testLargeSystemFusedPmeExclusions() {
    // Compare large single-context PME (fusePmeExclusions == true when paddedNumAtoms*40 > L2CacheSize)
    // against an identical active cluster below the L2 threshold and ReferencePlatform, covering water
    // triangles, pair exclusions (including coincident r == 0), 1-4 exceptions, and self-exclusion.
    System probe;
    probe.addParticle(1.0);
    VerletIntegrator probeInt(1.0);
    Context probeCtx(probe, probeInt, platform);
    CUdevice device;
    cuDeviceGet(&device, stoi(platform.getPropertyValue(probeCtx, CudaPlatform::CudaDeviceIndex())));
    size_t memory = 0;
    int l2Bytes = 0;
    cuDeviceTotalMem(&memory, device);
    cuDeviceGetAttribute(&l2Bytes, CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE, device);
    if (memory < 6ULL*(1<<30) || l2Bytes <= 0 || l2Bytes > 64*(1<<20))
        return;

    const int numDummy = ((l2Bytes/40 + 64 + 31)/32)*32;
    const int numDummyBlocks = numDummy/32;
    const int gridSide = (int) ceil(cbrt((double) numDummyBlocks)) + 1;
    const double spacing = 1.5;
    const double boxSize = 10.0 + gridSide*spacing;

    auto buildCluster = [&](System& sys, NonbondedForce*& nb, vector<Vec3>& pos) {
        sys.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
        nb = new NonbondedForce();
        nb->setNonbondedMethod(NonbondedForce::PME);
        nb->setCutoffDistance(0.8);
        nb->setPMEParameters(2.0, 32, 32, 32);
        nb->setUseDispersionCorrection(false);
        sys.addForce(nb);
        for (int w = 0; w < 8; w++) {
            int base = 3*w;
            Vec3 o(2.0 + 0.35*(w%3), 2.0 + 0.35*(w/3), 2.0);
            for (int k = 0; k < 3; k++) {
                sys.addParticle(k == 0 ? 16.0 : 1.0);
                nb->addParticle(k == 0 ? -0.834 : 0.417, k == 0 ? 0.315 : 0.1, 0.0);
            }
            pos.insert(pos.end(), {o, o + Vec3(0.09572, 0.0, 0.0), o + Vec3(-0.024, 0.0927, 0.0)});
            nb->addException(base, base+1, 0.0, 0.1, 0.0);
            nb->addException(base, base+2, 0.0, 0.1, 0.0);
            nb->addException(base+1, base+2, 0.0, 0.1, 0.0);
        }
        int c0 = sys.getNumParticles();
        const double charges[8] = {0.25, -0.25, 0.15, -0.25, 0.10, 0.0, 0.0, 0.0};
        for (int i = 0; i < 8; i++) {
            sys.addParticle(12.0);
            nb->addParticle(charges[i], 0.3, 0.2);
            pos.push_back(Vec3(2.0 + 0.14*(i == 4 ? 3 : i), 3.2, 2.0));
        }
        for (int i = 0; i < 7; i++)
            nb->addException(c0+i, c0+i+1, 0.0, 0.3, 0.0);
        for (int i = 0; i < 6; i++)
            nb->addException(c0+i, c0+i+2, 0.0, 0.3, 0.0);
        nb->addException(c0, c0+3, -0.03125, 0.3, 0.1);
        nb->addException(c0+1, c0+4, -0.0125, 0.3, 0.1);
        nb->addException(c0+7, c0+7, 0.0, 0.3, 0.0);
    };

    System smallSys, largeSys;
    NonbondedForce* smallNb = nullptr;
    NonbondedForce* largeNb = nullptr;
    vector<Vec3> smallPos, largePos;
    buildCluster(smallSys, smallNb, smallPos);
    largePos.reserve(32 + numDummy);
    buildCluster(largeSys, largeNb, largePos);
    for (int b = 0; b < numDummyBlocks; b++) {
        Vec3 p(8.0 + (b%gridSide)*spacing, 8.0 + ((b/gridSide)%gridSide)*spacing, 8.0 + (b/(gridSide*gridSide))*spacing);
        for (int k = 0; k < 32; k++) {
            largeSys.addParticle(1.0);
            largeNb->addParticle(0.0, 0.1, 0.0);
            largePos.push_back(p + Vec3(0.02*(k%4), 0.02*((k/4)%4), 0.02*(k/16)));
        }
    }

    ReferencePlatform refPlatform;
    VerletIntegrator intRef(0.001);
    VerletIntegrator intSmall(0.001);
    VerletIntegrator intLarge(0.001);
    Context ctxRef(smallSys, intRef, refPlatform);
    Context ctxSmall(smallSys, intSmall, platform);
    Context ctxLarge(largeSys, intLarge, platform);
    ctxRef.setPositions(smallPos);
    ctxSmall.setPositions(smallPos);
    ctxLarge.setPositions(largePos);
    State sRef = ctxRef.getState(State::Forces | State::Energy);
    State sSmall = ctxSmall.getState(State::Forces | State::Energy);
    State sLarge = ctxLarge.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(sRef.getPotentialEnergy(), sLarge.getPotentialEnergy(), 1e-3*fabs(sRef.getPotentialEnergy()) + 1e-2);
    ASSERT_EQUAL_TOL(sSmall.getPotentialEnergy(), sLarge.getPotentialEnergy(), 1e-4*fabs(sSmall.getPotentialEnergy()) + 1e-3);
    for (int i = 0; i < 32; i++) {
        ASSERT_EQUAL_VEC(sSmall.getForces()[i], sLarge.getForces()[i], 1e-3);
        ASSERT_EQUAL_VEC(sRef.getForces()[i], sLarge.getForces()[i], 2e-2);
    }
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), sLarge.getForces()[32], 1e-6);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), sLarge.getForces().back(), 1e-6);

    // Verify a coincident 1-4 exception (r == 0) still yields non-finite energy like ReferencePlatform.
    smallPos[24] = smallPos[27];
    largePos[24] = largePos[27];
    ctxRef.setPositions(smallPos);
    ctxLarge.setPositions(largePos);
    ASSERT(!isfinite(ctxRef.getState(State::Energy).getPotentialEnergy()));
    ASSERT(!isfinite(ctxLarge.getState(State::Energy).getPotentialEnergy()));
}

void runPlatformTests() {
    testLargeSystemFusedPmeExclusions();
}
