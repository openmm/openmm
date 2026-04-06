/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- */

// This file requires AVX2. It is only compiled on x86 with AVX2 flags.
#if defined(__AVX2__) || (defined(_MSC_VER) && defined(__AVX2__))

#include "CpuNonbondedForceCluster.h"
#include "SimTKOpenMMRealType.h"
#include "openmm/internal/hardware.h"
#include "openmm/internal/vectorizeAvx.h"
#include <algorithm>
#include <cmath>
#include <immintrin.h>

using namespace std;
using namespace OpenMM;

bool CpuNonbondedForceCluster::isSupported() {
    // Runtime AVX2 check via CPUID, cached after first call.
    static const bool supported = isAvx2Supported();
    return supported;
}

// Fast 8-wide horizontal sum: shuffle+add (6-cycle throughput vs 15 for dp_ps).
static inline float fastReduceAdd(fvec8 v) {
    __m128 lo = _mm256_castps256_ps128(v);
    __m128 hi = _mm256_extractf128_ps(v, 1);
    __m128 sum128 = _mm_add_ps(lo, hi);
    __m128 shuf = _mm_movehdup_ps(sum128);
    __m128 sum64 = _mm_add_ps(sum128, shuf);
    __m128 hi64 = _mm_movehl_ps(sum64, sum64);
    return _mm_cvtss_f32(_mm_add_ss(sum64, hi64));
}

/**
 * Flush accumulated i-cluster forces to the per-thread force buffer.
 */
static inline void flushIForces(const AtomCluster& ci, fvec8& fix, fvec8& fiy, fvec8& fiz, float* forces) {
    float f8x[8], f8y[8], f8z[8];
    fix.store(f8x); fiy.store(f8y); fiz.store(f8z);
    for (int ii = 0; ii < ci.size; ii++) {
        int iAtom = ci.atomIndex[ii];
        forces[4*iAtom]   += f8x[ii];
        forces[4*iAtom+1] += f8y[ii];
        forces[4*iAtom+2] += f8z[ii];
    }
}

/**
 * Process one j-cluster against pre-loaded i-cluster data.
 * This is the inner kernel that processes the j-atom loop only.
 * i-data (ix,iy,iz,iq,isigma,ieps) and i-force accumulators (fix,fiy,fiz) are
 * managed by the caller to enable grouped processing across consecutive pairs.
 */
template<bool USE_PERIODIC, bool USE_CUTOFF, bool USE_EWALD, bool USE_SWITCH, bool COMPUTE_ENERGY>
static inline void processJCluster(
    const fvec8& ix, const fvec8& iy, const fvec8& iz, const fvec8& iq,
    const fvec8& isigma, const fvec8& ieps,
    fvec8& fix, fvec8& fiy, fvec8& fiz,
    double& energy,
    const AtomCluster& __restrict cj,
    uint64_t exclMask, float* __restrict forces,
    float cutoff2, float bsx, float bsy, float bsz, float ibsx, float ibsy, float ibsz,
    float krf, float crf, float twoKrf,
    float alphaEwald,
    const float* erfcTable, const float* ewaldScaleTable,
    float ewaldDXInv, float erfcDXInv, int numTablePoints,
    uint8_t cp_jAtomMask,
    float switchingDistance, float invSwitchingInterval)
{
    double localEnergy = 0;

    // Pre-extract exclusion columns. Skip when exclMask==0 (98%+ of pairs).
    uint8_t exclColumns[CLUSTER_SIZE] = {};
    if (exclMask != 0) {
        const uint64_t colMask = 0x0101010101010101ULL;
        for (int jj = 0; jj < CLUSTER_SIZE; jj++)
            exclColumns[jj] = (uint8_t)_pext_u64(exclMask, colMask << jj);
    }

    fvec8 cutoff2Vec;
    if constexpr (USE_CUTOFF || USE_EWALD) {
        cutoff2Vec = fvec8(cutoff2);
    }

    // Process each active j-atom.
    unsigned int jMask = cp_jAtomMask;
    while (jMask) {
        int jj = _tzcnt_u32(jMask);
        jMask &= jMask - 1;

        fvec8 jx(cj.x[jj]), jy(cj.y[jj]), jz(cj.z[jj]);
        float jq_scalar = cj.q[jj];
        float jeps_scalar = cj.epsilon[jj];

        fvec8 dx = ix - jx;
        fvec8 dy = iy - jy;
        fvec8 dz = iz - jz;
        if constexpr (USE_PERIODIC) {
            dx -= round(dx * ibsx) * bsx;
            dy -= round(dy * ibsy) * bsy;
            dz -= round(dz * ibsz) * bsz;
        }
        fvec8 r2 = dx*dx + dy*dy + dz*dz;

        int exclBits = exclColumns[jj];

        // r2pos (r2>0) check only needed for exclusion path (self-pairs can have r2=0).
        // The 98% fast path (exclBits==0) never has r2=0 — saves 2 ops per j-atom.
        fvec8 include;
        if constexpr (USE_CUTOFF || USE_EWALD) {
            if (exclBits == 0) {
                include = fvec8(r2 < cutoff2Vec);
            } else {
                fvec8 r2pos(_mm256_cmp_ps(r2, _mm256_setzero_ps(), _CMP_GT_OQ));
                include = fvec8(_mm256_and_ps(blendZero(r2 < cutoff2Vec,
                    fvec8::expandBitsToMask(~exclBits)), r2pos));
            }
        } else {
            if (exclBits == 0) {
                include = fvec8(_mm256_set1_ps(-1.0f)); // all-true mask
            } else {
                fvec8 r2pos(_mm256_cmp_ps(r2, _mm256_setzero_ps(), _CMP_GT_OQ));
                include = fvec8(_mm256_and_ps(fvec8::expandBitsToMask(~exclBits), r2pos));
            }
        }
        // Note: early exit removed. blendZero(dEdR*inverseR2, include) handles
        // correctness for excluded/out-of-range pairs. Removing the branch
        // eliminates a 14% misprediction penalty (~2 cycles/j-atom average).
        // The 14% of j-atoms that would have been skipped now compute uselessly
        // but cost less than the branch misprediction overhead for all 100%.

        fvec8 inverseR = rsqrt(r2);
        fvec8 inverseR2 = inverseR * inverseR;

        fvec8 dEdR(0.0f), sig6(0.0f), epsSig6(0.0f);
        [[maybe_unused]] fvec8 pairEnergy(0.0f);
        if (jeps_scalar != 0.0f) {
            fvec8 sig = isigma + cj.sigma[jj];
            fvec8 eps = ieps * jeps_scalar;
            fvec8 sig2 = sig * sig * inverseR2;
            sig6 = sig2 * sig2 * sig2;
            sig6 = blend(fvec8(0.0f), sig6,
                         fvec8(_mm256_cmp_ps(eps, _mm256_setzero_ps(), _CMP_NEQ_OQ)));
            epsSig6 = eps * sig6;
            dEdR = epsSig6 * (12.0f * sig6 - 6.0f);
            if constexpr (COMPUTE_ENERGY)
                pairEnergy = epsSig6 * (sig6 - 1.0f);
        }

        fvec8 chargeProd = iq * jq_scalar;
        fvec8 r;
        if constexpr (USE_EWALD || USE_SWITCH) {
            r = r2 * inverseR;
        }
        if constexpr (USE_EWALD) {
            fvec8 ewaldScaleIdx = r * ewaldDXInv;
            ivec8 ewaldIntIdx = min(floor(ewaldScaleIdx), fvec8((float)numTablePoints));
            fvec8 ewaldFrac = ewaldScaleIdx - fvec8(_mm256_cvtepi32_ps(ewaldIntIdx));
            __m256i gIdx = ewaldIntIdx;
            __m256i gIdx1 = _mm256_add_epi32(gIdx, _mm256_set1_epi32(1));
            fvec8 es0(_mm256_i32gather_ps(ewaldScaleTable, gIdx, 4));
            fvec8 es1(_mm256_i32gather_ps(ewaldScaleTable, gIdx1, 4));
            fvec8 ewaldScale = es0 + ewaldFrac * (es1 - es0);
            dEdR = dEdR + chargeProd * inverseR * ewaldScale;
            if constexpr (COMPUTE_ENERGY) {
                // Only look up erfc table when energy is actually needed.
                fvec8 ef0(_mm256_i32gather_ps(erfcTable, gIdx, 4));
                fvec8 ef1(_mm256_i32gather_ps(erfcTable, gIdx1, 4));
                fvec8 erfc = ef0 + ewaldFrac * (ef1 - ef0);
                pairEnergy = pairEnergy + chargeProd * inverseR * erfc;
            }
        } else if constexpr (USE_CUTOFF) {
            dEdR = dEdR + chargeProd * (inverseR - twoKrf*r2);
            if constexpr (COMPUTE_ENERGY)
                pairEnergy = pairEnergy + chargeProd * (inverseR + krf*r2 - crf);
        } else {
            dEdR = dEdR + chargeProd * inverseR;
            if constexpr (COMPUTE_ENERGY)
                pairEnergy = pairEnergy + chargeProd * inverseR;
        }

        if constexpr (USE_SWITCH) {
            fvec8 sw_t = blendZero((r - switchingDistance) * invSwitchingInterval, r > fvec8(switchingDistance));
            fvec8 sw_t2 = sw_t * sw_t;
            fvec8 sw_t3 = sw_t2 * sw_t;
            fvec8 sw_S = 1.0f + sw_t3 * (-10.0f + sw_t * (15.0f - 6.0f * sw_t));
            fvec8 sw_dS = sw_t2 * (-30.0f + sw_t * (60.0f - 30.0f * sw_t)) * invSwitchingInterval;
            fvec8 sw_ljE = epsSig6 * (sig6 - 1.0f);
            fvec8 sw_ljF = epsSig6 * (12.0f * sig6 - 6.0f);
            dEdR = dEdR - sw_ljF + sw_ljF * sw_S - sw_ljE * sw_dS * r;
            if constexpr (COMPUTE_ENERGY)
                pairEnergy = pairEnergy - sw_ljE + sw_ljE * sw_S;
        }

        dEdR = blendZero(dEdR * inverseR2, include);

        fvec8 fx = dx * dEdR;
        fvec8 fy = dy * dEdR;
        fvec8 fz = dz * dEdR;

        fix = fix + fx;
        fiy = fiy + fy;
        fiz = fiz + fz;

        // Newton's third law on j-atom.
        // Note: jj < cj.size is guaranteed by jAtomMask (padding atoms have bit=0).
        {
            int jAtom = cj.atomIndex[jj];
            fvec4 jForce = fvec4(forces+4*jAtom) - reduceToVec3(fx, fy, fz);
            jForce.store(forces+4*jAtom);
        }

        if constexpr (COMPUTE_ENERGY) {
            pairEnergy = blendZero(pairEnergy, include);
            localEnergy += fastReduceAdd(pairEnergy);
        }
    }

    if constexpr (COMPUTE_ENERGY)
        energy += localEnergy;
}

/**
 * Grouped pair processing loop. Pairs are sorted by clusterI, so consecutive
 * pairs share the same i-cluster. This function:
 * - Loads i-data once per clusterI group (not per pair)
 * - Accumulates i-forces in registers across the group
 * - Stores i-forces once per group
 * This eliminates ~99% of i-data loads and i-force stores (163 pairs/group avg).
 */
template<bool P, bool C, bool E, bool S, bool EN>
static void runGroupedLoop(
    const ClusterPair* __restrict pairs, const AtomCluster* __restrict clusters,
    int numPairs, float* __restrict forces, double& energy,
    float cutoff2, float bsx, float bsy, float bsz, float ibsx, float ibsy, float ibsz,
    float krf, float crf,
    float alphaEwald,
    const float* erfcTable, const float* ewaldScaleTable,
    float ewaldDXInv, float erfcDXInv, int numTablePoints,
    float switchingDistance, float invSwitchInterval,
    atomic<int>& atomicCounter)
{
    const int CHUNK_SIZE = 128;
    const float twoKrf = 2.0f * krf;

    int prevI = -1;
    fvec8 ix, iy, iz, iq, isigma, ieps;
    fvec8 fix(0.0f), fiy(0.0f), fiz(0.0f);
    const AtomCluster* ciPtr = nullptr;

    while (true) {
        int chunkStart = atomicCounter.fetch_add(CHUNK_SIZE, std::memory_order_relaxed);
        if (chunkStart >= numPairs) break;
        int chunkEnd = min(chunkStart + CHUNK_SIZE, numPairs);

        for (int p = chunkStart; p < chunkEnd; p++) {
            const ClusterPair& cp = pairs[p];
            int curI = cp.clusterI;

            if (curI != prevI) {
                // Flush previous group's i-forces.
                if (ciPtr != nullptr)
                    flushIForces(*ciPtr, fix, fiy, fiz, forces);

                // Load new i-cluster data into registers.
                ciPtr = &clusters[curI];
                ix = fvec8(ciPtr->x); iy = fvec8(ciPtr->y); iz = fvec8(ciPtr->z);
                iq = fvec8(ciPtr->q) * (float)ONE_4PI_EPS0;
                isigma = fvec8(ciPtr->sigma); ieps = fvec8(ciPtr->epsilon);
                fix = fvec8(0.0f); fiy = fvec8(0.0f); fiz = fvec8(0.0f);
                prevI = curI;
            }

            processJCluster<P,C,E,S,EN>(
                ix, iy, iz, iq, isigma, ieps, fix, fiy, fiz, energy,
                clusters[cp.clusterJ], cp.exclusionMask, forces,
                cutoff2, bsx, bsy, bsz, ibsx, ibsy, ibsz,
                krf, crf, twoKrf, alphaEwald, erfcTable, ewaldScaleTable,
                ewaldDXInv, erfcDXInv, numTablePoints,
                cp.jAtomMask, switchingDistance, invSwitchInterval);
        }
    }

    // Flush final group.
    if (ciPtr != nullptr)
        flushIForces(*ciPtr, fix, fiy, fiz, forces);
}

/**
 * Legacy per-pair kernel (used by computeClusterPairIxn for tests).
 * Loads i-data, processes j-cluster, stores i-forces - all in one call.
 */
template<bool USE_PERIODIC, bool USE_CUTOFF, bool USE_EWALD, bool USE_SWITCH, bool COMPUTE_ENERGY>
static void computeClusterPairIxnImpl(
    const AtomCluster& __restrict ci, const AtomCluster& __restrict cj,
    uint64_t exclMask, float* __restrict forces, double& energy,
    float cutoff2, const fvec4& boxSize, const fvec4& invBoxSize,
    float krf, float crf,
    float alphaEwald,
    const float* erfcTable, const float* ewaldScaleTable,
    float ewaldDXInv, float erfcDXInv, int numTablePoints,
    uint8_t cp_jAtomMask,
    float switchingDistance, float invSwitchingInterval)
{
    fvec8 ix(ci.x), iy(ci.y), iz(ci.z), iq(ci.q);
    fvec8 isigma(ci.sigma), ieps(ci.epsilon);
    iq = iq * (float)ONE_4PI_EPS0;
    fvec8 fix(0.0f), fiy(0.0f), fiz(0.0f);

    float bsx = 0, bsy = 0, bsz = 0, ibsx = 0, ibsy = 0, ibsz = 0;
    if constexpr (USE_PERIODIC) {
        bsx = boxSize[0]; bsy = boxSize[1]; bsz = boxSize[2];
        ibsx = invBoxSize[0]; ibsy = invBoxSize[1]; ibsz = invBoxSize[2];
    }

    processJCluster<USE_PERIODIC, USE_CUTOFF, USE_EWALD, USE_SWITCH, COMPUTE_ENERGY>(
        ix, iy, iz, iq, isigma, ieps, fix, fiy, fiz, energy,
        cj, exclMask, forces, cutoff2, bsx, bsy, bsz, ibsx, ibsy, ibsz,
        krf, crf, 2.0f*krf, alphaEwald, erfcTable, ewaldScaleTable,
        ewaldDXInv, erfcDXInv, numTablePoints,
        cp_jAtomMask, switchingDistance, invSwitchingInterval);

    flushIForces(ci, fix, fiy, fiz, forces);
}

// Dispatch: select the correct template instantiation at runtime.
void CpuNonbondedForceCluster::computeClusterPairIxn(
    const AtomCluster& ci, const AtomCluster& cj,
    uint64_t exclMask, float* forces, double& energy,
    float cutoff2, const fvec4& boxSize, const fvec4& invBoxSize,
    bool usePeriodic, bool useCutoff, float krf, float crf,
    bool useEwald, float alphaEwald,
    const float* erfcTable, const float* ewaldScaleTable,
    float ewaldDXInv, float erfcDXInv, int numTablePoints,
    uint8_t cp_jAtomMask)
{
    if (useEwald) {
        computeClusterPairIxnImpl<true, true, true, false, true>(
            ci, cj, exclMask, forces, energy, cutoff2, boxSize, invBoxSize,
            krf, crf, alphaEwald, erfcTable, ewaldScaleTable, ewaldDXInv, erfcDXInv, numTablePoints, cp_jAtomMask, 0, 0);
    } else if (useCutoff && usePeriodic) {
        computeClusterPairIxnImpl<true, true, false, false, true>(
            ci, cj, exclMask, forces, energy, cutoff2, boxSize, invBoxSize,
            krf, crf, 0, nullptr, nullptr, 0, 0, 0, cp_jAtomMask, 0, 0);
    } else if (useCutoff) {
        computeClusterPairIxnImpl<false, true, false, false, true>(
            ci, cj, exclMask, forces, energy, cutoff2, boxSize, invBoxSize,
            krf, crf, 0, nullptr, nullptr, 0, 0, 0, cp_jAtomMask, 0, 0);
    } else {
        computeClusterPairIxnImpl<false, false, false, false, true>(
            ci, cj, exclMask, forces, energy, cutoff2, boxSize, invBoxSize,
            0, 0, 0, nullptr, nullptr, 0, 0, 0, cp_jAtomMask, 0, 0);
    }
}

void CpuNonbondedForceCluster::calculateDirectIxn(
    const CpuClusterPairList& clusterPairList,
    const float* posq,
    vector<AlignedArray<float>>& threadForce,
    double* totalEnergy,
    ThreadPool& threads,
    float cutoffDistance,
    bool usePeriodic,
    const Vec3* periodicBoxVectors,
    bool useCutoff,
    float krf,
    float crf,
    bool useEwald,
    float alphaEwald,
    const float* erfcTable,
    const float* ewaldScaleTable,
    float ewaldDXInv,
    float erfcDXInv,
    int numTablePoints,
    bool useSwitch,
    float switchingDistance)
{
    float cutoff2 = cutoffDistance * cutoffDistance;
    float bsx = 0, bsy = 0, bsz = 0, ibsx = 0, ibsy = 0, ibsz = 0;
    if (usePeriodic) {
        bsx = (float)periodicBoxVectors[0][0]; bsy = (float)periodicBoxVectors[1][1];
        bsz = (float)periodicBoxVectors[2][2];
        ibsx = 1.0f/bsx; ibsy = 1.0f/bsy; ibsz = 1.0f/bsz;
    }

    const vector<ClusterPair>& pairs = clusterPairList.getClusterPairs();
    const vector<AtomCluster>& clusters = clusterPairList.getClusters();
    int numPairs = pairs.size();
    int numThreads = threads.getNumThreads();

    vector<double> threadEnergy(numThreads, 0.0);
    atomic<int> atomicCounter(0);

    bool needEnergy = (totalEnergy != NULL);
    float invSwitchInterval = (cutoffDistance > switchingDistance && switchingDistance > 0) ?
                               1.0f / (cutoffDistance - switchingDistance) : 0;

    // Select the correct grouped loop template. The loop processes consecutive pairs
    // sharing the same clusterI as a group, keeping i-data in registers and accumulating
    // i-forces without writing to memory until the group changes.
    using LoopFn = void(*)(const ClusterPair*, const AtomCluster*, int, float*, double&,
                            float, float, float, float, float, float, float,
                            float, float, float,
                            const float*, const float*, float, float, int,
                            float, float, atomic<int>&);

    LoopFn loopFn;
    if (needEnergy) {
        if (useEwald && useSwitch)
            loopFn = &runGroupedLoop<true, true, true, true, true>;
        else if (useEwald)
            loopFn = &runGroupedLoop<true, true, true, false, true>;
        else if (useCutoff && usePeriodic && useSwitch)
            loopFn = &runGroupedLoop<true, true, false, true, true>;
        else if (useCutoff && usePeriodic)
            loopFn = &runGroupedLoop<true, true, false, false, true>;
        else if (useCutoff)
            loopFn = &runGroupedLoop<false, true, false, false, true>;
        else
            loopFn = &runGroupedLoop<false, false, false, false, true>;
    } else {
        if (useEwald && useSwitch)
            loopFn = &runGroupedLoop<true, true, true, true, false>;
        else if (useEwald)
            loopFn = &runGroupedLoop<true, true, true, false, false>;
        else if (useCutoff && usePeriodic && useSwitch)
            loopFn = &runGroupedLoop<true, true, false, true, false>;
        else if (useCutoff && usePeriodic)
            loopFn = &runGroupedLoop<true, true, false, false, false>;
        else if (useCutoff)
            loopFn = &runGroupedLoop<false, true, false, false, false>;
        else
            loopFn = &runGroupedLoop<false, false, false, false, false>;
    }

    threads.execute([&](ThreadPool& threads, int threadIndex) {
        float* forces = &threadForce[threadIndex][0];
        double myEnergy = 0.0;
        loopFn(pairs.data(), clusters.data(), numPairs, forces, myEnergy,
               cutoff2, bsx, bsy, bsz, ibsx, ibsy, ibsz,
               krf, crf, alphaEwald, erfcTable, ewaldScaleTable,
               ewaldDXInv, erfcDXInv, numTablePoints,
               switchingDistance, invSwitchInterval, atomicCounter);
        threadEnergy[threadIndex] = myEnergy;
    });
    threads.waitForThreads();

    if (totalEnergy != NULL)
        for (int t = 0; t < numThreads; t++)
            *totalEnergy += threadEnergy[t];
}

#endif // __AVX2__
