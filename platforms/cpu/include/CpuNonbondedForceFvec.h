
/* Portions copyright (c) 2006-2015 Stanford University and Simbios.
 * Contributors: Daniel Towner
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

#ifndef OPENMM_CPU_NONBONDED_FORCE_FVEC_H__
#define OPENMM_CPU_NONBONDED_FORCE_FVEC_H__

#include "CpuNonbondedForce.h"
#include "openmm/internal/vectorize.h"

#include "SimTKOpenMMUtilities.h"

#include <algorithm>
#include <vector>

namespace OpenMM {

enum BlockType {EWALD, NON_EWALD}; // :TODO: Better name for non-ewald.
enum PeriodicType {NoCutoff, NoPeriodic, PeriodicPerAtom, PeriodicPerInteraction, PeriodicTriclinic};

/**
 * Generic SIMD implementation of CpuNonbondedForce. The templating allows the same
 * basic code to be reused for any sort of SIMD type, including SSE, AVX, AVX2, or
 * AVX-512.
 */
template<typename FVEC>
class CpuNonbondedForceFvec : public CpuNonbondedForce {
public:
    /**
     * Store how many elements are contained in each block of atoms.
     */
    static constexpr int blockSize = sizeof(FVEC) / sizeof(float);
    /**
     * Constructor.
     */
    CpuNonbondedForceFvec(const CpuNeighborList& neighbors);

protected:
    /**---------------------------------------------------------------------------------------
      Calculate all the interactions for one atom block. These are part of the virtual function interface
      and consequently have names which explicitly call Ewald variant or not.
      They internally call into the generic handler function below.
      @param blockIndex       the index of the atom block
      @param forces           force array (forces added)
      @param totalEnergy      total energy
      --------------------------------------------------------------------------------------- 
      @{
      */
    void calculateBlockIxn(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize);
    void calculateBlockEwaldIxn(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize);
    /** @} */

    /**---------------------------------------------------------------------------------------
      Calculate all the interactions for one atom block. Identical to function prototypes above but
      with an extra template parameter to choose whether to use Ewald processing or not.
      --------------------------------------------------------------------------------------- */
    template<BlockType BLOCK_TYPE>
    void calculateBlockIxnHandler(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
    * Templatized implementation of calculateBlockIxn. It can handle both Ewald and non-ewald interactions
    * through a template parameter since the code is so similar for the two cases. Note also that the
    * floating-point SIMD type is also templated to allow any suitable type to be used.
    */
    template <int PERIODIC_TYPE, BlockType BLOCK_TYPE>
    void calculateBlockIxnImpl(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize, const fvec4& blockCenter);

    /**
     * Compute the displacement and squared distance between a collection of points, optionally using
     * periodic boundary conditions.
     */
    template <int PERIODIC_TYPE>
    void getDeltaR(const fvec4& posI, const FVEC& x, const FVEC& y, const FVEC& z, FVEC& dx, FVEC& dy, FVEC& dz, FVEC& r2, const fvec4& boxSize, const fvec4& invBoxSize) const;

    /**
     * Compute an approximation of a function using a table lookup.
     **/
      FVEC approximateFunctionFromTable(const std::vector<float>& table, FVEC x, FVEC inverse) const;

};

template <typename FVEC>
CpuNonbondedForceFvec<FVEC>::CpuNonbondedForceFvec(const CpuNeighborList& neighbors) : CpuNonbondedForce(neighbors) {
}

/**
 * Use a table lookup to approximate a function specific function.
 */
template<typename FVEC>
FVEC
CpuNonbondedForceFvec<FVEC>::approximateFunctionFromTable(const std::vector<float>& table,
                                                          const FVEC x, const FVEC inverse) const {
    // Compute the set of 8 index positions from which to gather the table data.
    const auto x1 = x * inverse;
    const auto index = min(floor(x1), float(NUM_TABLE_POINTS));

    FVEC s1, s2;
    gatherVecPair(table.data(), index, s1, s2);

    const auto coeff2 = x1 - FVEC(index);
    return (s1 - coeff2 * s1) + coeff2 * s2;
}

template<typename FVEC>
void CpuNonbondedForceFvec<FVEC>::calculateBlockIxn(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    calculateBlockIxnHandler<BlockType::NON_EWALD>(blockIndex, forces, totalEnergy, boxSize, invBoxSize);
}

template<typename FVEC>
void CpuNonbondedForceFvec<FVEC>::calculateBlockEwaldIxn(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    calculateBlockIxnHandler<BlockType::EWALD>(blockIndex, forces, totalEnergy, boxSize, invBoxSize);
}

template<typename FVEC>
template<BlockType BLOCK_TYPE>
void CpuNonbondedForceFvec<FVEC>::calculateBlockIxnHandler(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Determine whether we need to apply periodic boundary conditions.

    PeriodicType periodicType;
    fvec4 blockCenter;
    if (!periodic) {
        periodicType = NoPeriodic;
        blockCenter = 0.0f;
    }
    else {
        using std::min;
        using std::max;

        const int32_t* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
        float minx, maxx, miny, maxy, minz, maxz;
        minx = maxx = posq[4*blockAtom[0]];
        miny = maxy = posq[4*blockAtom[0]+1];
        minz = maxz = posq[4*blockAtom[0]+2];
        for (int i = 1; i < blockSize; i++) {
            minx = min(minx, posq[4*blockAtom[i]]);
            maxx = max(maxx, posq[4*blockAtom[i]]);
            miny = min(miny, posq[4*blockAtom[i]+1]);
            maxy = max(maxy, posq[4*blockAtom[i]+1]);
            minz = min(minz, posq[4*blockAtom[i]+2]);
            maxz = max(maxz, posq[4*blockAtom[i]+2]);
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
        calculateBlockIxnImpl<NoCutoff, BLOCK_TYPE>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == NoPeriodic)
        calculateBlockIxnImpl<NoPeriodic, BLOCK_TYPE>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicPerAtom)
        calculateBlockIxnImpl<PeriodicPerAtom, BLOCK_TYPE>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicPerInteraction)
        calculateBlockIxnImpl<PeriodicPerInteraction, BLOCK_TYPE>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicTriclinic)
        calculateBlockIxnImpl<PeriodicTriclinic, BLOCK_TYPE>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
}

template<typename FVEC>
template <int PERIODIC_TYPE, BlockType BLOCK_TYPE>
void CpuNonbondedForceFvec<FVEC>::calculateBlockIxnImpl(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize, const fvec4& blockCenter) {
    // Load the positions and parameters of the atoms in the block.

    const int32_t* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
    fvec4 blockAtomPosq[blockSize];
    FVEC blockAtomForceX(0.0f), blockAtomForceY(0.0f), blockAtomForceZ(0.0f);
    FVEC blockAtomX, blockAtomY, blockAtomZ, blockAtomCharge;
    for (int i = 0; i < blockSize; i++) {
        blockAtomPosq[i] = fvec4(posq+4*blockAtom[i]);
        if (PERIODIC_TYPE == PeriodicPerAtom)
            blockAtomPosq[i] -= floor((blockAtomPosq[i]-blockCenter)*invBoxSize+0.5f)*boxSize; // :TODO: Apply one to blockAtom?
    }

    transpose(blockAtomPosq, blockAtomX, blockAtomY, blockAtomZ, blockAtomCharge);
    blockAtomCharge *= ONE_4PI_EPS0;

    // Not the most efficient way to do this, but it works across all types we care about, and this isn't where
    // the cycles are spent anyway.
    FVEC blockAtomSigma = {};
    FVEC blockAtomEpsilon = {};
    for (int i = 0; i < blockSize; ++i) {
        ((float*)&blockAtomSigma)[i] = atomParameters[blockAtom[i]].first;
        ((float*)&blockAtomEpsilon)[i] = atomParameters[blockAtom[i]].second;
    }

    // Ewald needs C6 data gathered from a table. Unused variable for non-ewald.
    const FVEC C6s = (BLOCK_TYPE == BlockType::EWALD) ? FVEC(C6params, blockAtom) : FVEC();

    const float invSwitchingInterval = 1/(cutoffDistance-switchingDistance);
    const FVEC cutoffDistanceSquared = cutoffDistance * cutoffDistance;

    // Loop over neighbors for this block.
    CpuNeighborList::NeighborIterator neighbors = neighborList->getNeighborIterator(blockIndex);
    FVEC partialEnergy = {};
    while (neighbors.next()) {
        // Load the next neighbor.

        int atom = neighbors.getNeighbor();

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
        FVEC energy, dEdR;
        float atomEpsilon = atomParameters[atom].second;
        if (atomEpsilon != 0.0f) {
            const auto sig = blockAtomSigma+atomParameters[atom].first;
            const auto sig2 = (inverseR*sig)*(inverseR*sig);
            const auto sig6 = sig2*sig2*sig2;
            const auto eps = blockAtomEpsilon*atomEpsilon;
            const auto epsSig6 = eps*sig6;
            dEdR = epsSig6*(12.0f*sig6 - 6.0f);
            energy = epsSig6*(sig6-1.0f);
            if (useSwitch) {
                const auto t = blendZero((r-switchingDistance)*invSwitchingInterval, r>switchingDistance);
                const auto switchValue = 1+t*t*t*(-10.0f+t*(15.0f-t*6.0f));
                const auto switchDeriv = t*t*(-30.0f+t*(60.0f-t*30.0f))*invSwitchingInterval;
                dEdR = switchValue*dEdR - energy*switchDeriv*r;
                energy *= switchValue;
            }
            if (BLOCK_TYPE == BlockType::EWALD && ljpme) {
                const auto C6ij = C6s*C6params[atom];
                const auto inverseR2 = inverseR*inverseR;
                const auto mysig2 = sig*sig;
                const auto mysig6 = mysig2*mysig2*mysig2;
                const auto emult = C6ij*inverseR2*inverseR2*inverseR2*approximateFunctionFromTable(exptermsTable, r, FVEC(exptermsDXInv));
                const auto potentialShift = eps*(1.0f-mysig6*inverseRcut6)*mysig6*inverseRcut6 - C6ij*inverseRcut6Expterm;
                dEdR += 6.0f*C6ij*inverseR2*inverseR2*inverseR2*approximateFunctionFromTable(dExptermsTable, r, FVEC(exptermsDXInv));
                energy += emult + potentialShift;
            }

        }
        else {
            energy = 0.0f;
            dEdR = 0.0f;
        }
        const auto chargeProd = blockAtomCharge*posq[4*atom+3];
        if (BLOCK_TYPE == BlockType::EWALD) {
            dEdR += chargeProd*inverseR*approximateFunctionFromTable(ewaldScaleTable, r, FVEC(ewaldDXInv));
        }
        else {
            if (cutoff)
                dEdR += chargeProd*(inverseR-2.0f*krf*r2);
            else
                dEdR += chargeProd*inverseR;
        }
        dEdR *= inverseR*inverseR;

        // Accumulate energies.
        if (totalEnergy) {
            if (BLOCK_TYPE == BlockType::EWALD)
                energy += chargeProd*inverseR*approximateFunctionFromTable(erfcTable, alphaEwald*r, FVEC(erfcDXInv));
            else {  // Non-ewald.
                if (cutoff)
                    energy += chargeProd*(inverseR+krf*r2-crf);
                else
                    energy += chargeProd*inverseR;
            }
            energy = blendZero(energy, include);

            partialEnergy += energy;
        }

        // Accumulate forces.
        dEdR = blendZero(dEdR, include);
        const auto fx = dx*dEdR;
        const auto fy = dy*dEdR;
        const auto fz = dz*dEdR;
        blockAtomForceX += fx;
        blockAtomForceY += fy;
        blockAtomForceZ += fz;

        float* const atomForce = forces+4*atom;
        const fvec4 newAtomForce = fvec4(atomForce) - reduceToVec3(fx, fy, fz);
        newAtomForce.store(atomForce);
    }
    
    if (totalEnergy)
        *totalEnergy += reduceAdd(partialEnergy);

    // Record the forces on the block atoms.
    fvec4 f[blockSize];
    transpose(blockAtomForceX, blockAtomForceY, blockAtomForceZ, 0.0f, f);
    for (int j = 0; j < blockSize; j++)
        (fvec4(forces+4*blockAtom[j])+f[j]).store(forces+4*blockAtom[j]);
}

template<typename FVEC>
template <int PERIODIC_TYPE>
void CpuNonbondedForceFvec<FVEC>::getDeltaR(const fvec4& posI, const FVEC& x, const FVEC& y, const FVEC& z, FVEC& dx, FVEC& dy, FVEC& dz, FVEC& r2, const fvec4& boxSize, const fvec4& invBoxSize) const {
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

#endif // OPENMM_CPU_NONBONDED_FORCE_FVEC_H__
