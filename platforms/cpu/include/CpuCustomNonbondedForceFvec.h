
/* Portions copyright (c) 2009-2022 Stanford University and Simbios.
 * Contributors: Peter Eastman
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef OPENMM_CPU_CUSTOM_NONBONDED_FORCE_FVEC_H__
#define OPENMM_CPU_CUSTOM_NONBONDED_FORCE_FVEC_H__

#include "CpuCustomNonbondedForce.h"

namespace OpenMM {

enum PeriodicType {NoCutoff, NoPeriodic, PeriodicPerAtom, PeriodicPerInteraction, PeriodicTriclinic};

template <typename FVEC, int BLOCK_SIZE>
class CpuCustomNonbondedForceFvec : public CpuCustomNonbondedForce {
public:
    CpuCustomNonbondedForceFvec(ThreadPool& threads, const CpuNeighborList& neighbors) : CpuCustomNonbondedForce(threads, neighbors) {
    }

protected:
    /**
     * Calculate all the interactions for one block of atoms.
     * 
     * @param data            workspace for the current thread
     * @param blockIndex      the index of the atom block
     * @param forces          force array (forces added)
     * @param totalEnergy     total energy
     * @param boxSize         the size of the periodic box
     * @param invBoxSize       the inverse size of the periodic box
     */
    void calculateBlockIxn(ThreadData& data, int blockIndex, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize);

    template <int PERIODIC_TYPE>
    void calculateBlockIxnImpl(ThreadData& data, int blockIndex, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize, const fvec4& blockCenter);

    /**
     * Compute the displacement and squared distance between a collection of points, optionally using
     * periodic boundary conditions.
     */
    template <int PERIODIC_TYPE>
    void getDeltaR(const fvec4& posI, const FVEC& x, const FVEC& y, const FVEC& z, FVEC& dx, FVEC& dy, FVEC& dz, FVEC& r2, const fvec4& boxSize, const fvec4& invBoxSize) const;
};

template<typename FVEC, int BLOCK_SIZE>
void CpuCustomNonbondedForceFvec<FVEC, BLOCK_SIZE>::calculateBlockIxn(ThreadData& data, int blockIndex, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Determine whether we need to apply periodic boundary conditions.

    PeriodicType periodicType;
    fvec4 blockCenter;
    if (!periodic) {
        periodicType = NoPeriodic;
        blockCenter = 0.0f;
    }
    else {
        const int32_t* blockAtom = &neighborList->getSortedAtoms()[BLOCK_SIZE*blockIndex];
        float minx, maxx, miny, maxy, minz, maxz;
        minx = maxx = posq[4*blockAtom[0]];
        miny = maxy = posq[4*blockAtom[0]+1];
        minz = maxz = posq[4*blockAtom[0]+2];
        for (int i = 1; i < BLOCK_SIZE; i++) {
            minx = std::min(minx, posq[4*blockAtom[i]]);
            maxx = std::max(maxx, posq[4*blockAtom[i]]);
            miny = std::min(miny, posq[4*blockAtom[i]+1]);
            maxy = std::max(maxy, posq[4*blockAtom[i]+1]);
            minz = std::min(minz, posq[4*blockAtom[i]+2]);
            maxz = std::max(maxz, posq[4*blockAtom[i]+2]);
        }
        blockCenter = fvec4(0.5f*(minx+maxx), 0.5f*(miny+maxy), 0.5f*(minz+maxz), 0.0f);
        if (!(minx < cutoffDistance || miny < cutoffDistance || minz < cutoffDistance ||
                maxx > boxSize[0]-cutoffDistance || maxy > boxSize[1]-cutoffDistance || maxz > boxSize[2]-cutoffDistance))
            periodicType = NoPeriodic;
        else if (triclinic)
            periodicType = PeriodicTriclinic;
        else if (0.5f*(boxSize[0]-(maxx-minx)) >= cutoffDistance &&
                 0.5f*(boxSize[1]-(maxy-miny)) >= cutoffDistance &&
                 0.5f*(boxSize[2]-(maxz-minz)) >= cutoffDistance)
            periodicType = PeriodicPerAtom;
        else
            periodicType = PeriodicPerInteraction;
    }

    // Call the appropriate version depending on what calculation is required for periodic boundary conditions.

    if (!cutoff)
        calculateBlockIxnImpl<NoCutoff>(data, blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == NoPeriodic)
        calculateBlockIxnImpl<NoPeriodic>(data, blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicPerAtom)
        calculateBlockIxnImpl<PeriodicPerAtom>(data, blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicPerInteraction)
        calculateBlockIxnImpl<PeriodicPerInteraction>(data, blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicTriclinic)
        calculateBlockIxnImpl<PeriodicTriclinic>(data, blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
}

template<typename FVEC, int BLOCK_SIZE>
template <int PERIODIC_TYPE>
void CpuCustomNonbondedForceFvec<FVEC, BLOCK_SIZE>::calculateBlockIxnImpl(ThreadData& data, int blockIndex, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize, const fvec4& blockCenter) {
    // Load the positions and parameters of the atoms in the block.

    const int32_t* blockAtom = &neighborList->getSortedAtoms()[BLOCK_SIZE*blockIndex];
    fvec4 blockAtomPosq[BLOCK_SIZE];
    FVEC blockAtomForceX(0.0f), blockAtomForceY(0.0f), blockAtomForceZ(0.0f);
    FVEC blockAtomX, blockAtomY, blockAtomZ, blockAtomCharge;
    int numParams = paramNames.size();
    int numComputed = computedValueNames.size();
    for (int i = 0; i < BLOCK_SIZE; i++) {
        blockAtomPosq[i] = fvec4(posq+4*blockAtom[i]);
        if (PERIODIC_TYPE == PeriodicPerAtom)
            blockAtomPosq[i] -= floor((blockAtomPosq[i]-blockCenter)*invBoxSize+0.5f)*boxSize;
        for (int j = 0; j < numParams; j++)
            data.vecParticle1Params[j*BLOCK_SIZE+i] = atomParameters[blockAtom[i]][j];
        for (int j = 0; j < numComputed; j++)
            data.vecParticle1Values[j*BLOCK_SIZE+i] = atomComputedValues[j][blockAtom[i]];
    }
    transpose(blockAtomPosq, blockAtomX, blockAtomY, blockAtomZ, blockAtomCharge);
    const float invSwitchingInterval = 1/(cutoffDistance-switchingDistance);
    const FVEC cutoffDistanceSquared = cutoffDistance * cutoffDistance;

    // Loop over neighbors for this block.

    CpuNeighborList::NeighborIterator neighbors = neighborList->getNeighborIterator(blockIndex);
    FVEC partialEnergy = {};
    while (neighbors.next()) {
        // Load the next neighbor.

        int atom = neighbors.getNeighbor();
        for (int j = 0; j < numParams; j++)
            for (int k = 0; k < BLOCK_SIZE; k++)
                data.vecParticle2Params[j*BLOCK_SIZE+k] = atomParameters[atom][j];
        for (int j = 0; j < numComputed; j++)
            for (int k = 0; k < BLOCK_SIZE; k++)
                data.vecParticle2Values[j*BLOCK_SIZE+k] = atomComputedValues[j][atom];

        // Compute the distances to the block atoms.

        FVEC dx, dy, dz, r2;
        fvec4 atomPos(posq+4*atom);
        if (PERIODIC_TYPE == PeriodicPerAtom)
            atomPos -= floor((atomPos-blockCenter)*invBoxSize+0.5f)*boxSize;
        getDeltaR<PERIODIC_TYPE>(atomPos, blockAtomX, blockAtomY, blockAtomZ, dx, dy, dz, r2, boxSize, invBoxSize);
        auto include = FVEC::expandBitsToMask(~neighbors.getExclusions());
        if (PERIODIC_TYPE != NoCutoff)
            include = blendZero(r2 < cutoffDistanceSquared, include);
        if (!any(include))
            continue; // No interactions to compute.

        // Compute the interactions.

        const auto inverseR = rsqrt(r2);
        const auto r = r2*inverseR;
        r.store(data.rvec.data());
        FVEC dEdR(data.forceVecExpression.evaluate());
        FVEC energy;
        if (includeEnergy || useSwitch)
            energy = FVEC(data.energyVecExpression.evaluate());
        if (useSwitch) {
            const auto t = blendZero((r-switchingDistance)*invSwitchingInterval, r>switchingDistance);
            const auto switchValue = 1+t*t*t*(-10.0f+t*(15.0f-t*6.0f));
            const auto switchDeriv = t*t*(-30.0f+t*(60.0f-t*30.0f))*invSwitchingInterval;
            dEdR = switchValue*dEdR + energy*switchDeriv;
            energy *= switchValue;
        }
        dEdR *= inverseR;

        // Accumulate forces and energies.

        if (includeEnergy) {
            energy = blendZero(energy, include);
            partialEnergy += energy;
        }
        dEdR = blendZero(dEdR, include);
        const auto fx = dx*dEdR;
        const auto fy = dy*dEdR;
        const auto fz = dz*dEdR;
        blockAtomForceX -= fx;
        blockAtomForceY -= fy;
        blockAtomForceZ -= fz;
        float* const atomForce = forces+4*atom;
        const fvec4 newAtomForce = fvec4(atomForce) + reduceToVec3(fx, fy, fz);
        newAtomForce.store(atomForce);
    }
    if (includeEnergy)
        totalEnergy += reduceAdd(partialEnergy);

    // Record the forces on the block atoms.

    fvec4 f[BLOCK_SIZE];
    transpose(blockAtomForceX, blockAtomForceY, blockAtomForceZ, 0.0f, f);
    for (int j = 0; j < BLOCK_SIZE; j++)
        (fvec4(forces+4*blockAtom[j])+f[j]).store(forces+4*blockAtom[j]);
}

template<typename FVEC, int BLOCK_SIZE>
template <int PERIODIC_TYPE>
void CpuCustomNonbondedForceFvec<FVEC, BLOCK_SIZE>::getDeltaR(const fvec4& posI, const FVEC& x, const FVEC& y, const FVEC& z, FVEC& dx, FVEC& dy, FVEC& dz, FVEC& r2, const fvec4& boxSize, const fvec4& invBoxSize) const {
    dx = x-posI[0];
    dy = y-posI[1];
    dz = z-posI[2];
    if (PERIODIC_TYPE == PeriodicTriclinic) {
        const auto scale3 = floor(dz*recipBoxSize[2]+0.5f);
        dx -= scale3*periodicBoxVectors[2][0];
        dy -= scale3*periodicBoxVectors[2][1];
        dz -= scale3*periodicBoxVectors[2][2];
        const auto scale2 = floor(dy*recipBoxSize[1]+0.5f);
        dx -= scale2*periodicBoxVectors[1][0];
        dy -= scale2*periodicBoxVectors[1][1];
        const auto scale1 = floor(dx*recipBoxSize[0]+0.5f);
        dx -= scale1*periodicBoxVectors[0][0];
    }
    else if (PERIODIC_TYPE == PeriodicPerInteraction) {
        dx -= round(dx*invBoxSize[0])*boxSize[0];
        dy -= round(dy*invBoxSize[1])*boxSize[1];
        dz -= round(dz*invBoxSize[2])*boxSize[2];
    }
    r2 = dx*dx + dy*dy + dz*dz;
}

} // namespace OpenMM

#endif // OPENMM_CPU_CUSTOM_NONBONDED_FORCE_FVEC_H__
