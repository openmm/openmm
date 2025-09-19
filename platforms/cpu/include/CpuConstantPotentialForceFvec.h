#ifndef OPENMM_CPUCONSTANTPOTENTIALFORCEFVEC_H_
#define OPENMM_CPUCONSTANTPOTENTIALFORCEFVEC_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Evan Pretti                                                       *
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

#include "CpuConstantPotentialForce.h"
#include "SimTKOpenMMRealType.h"

namespace OpenMM {

enum PeriodicType {NoPeriodic, PeriodicPerAtom, PeriodicPerInteraction, PeriodicTriclinic};

template<typename FVEC, typename IVEC>
class CpuConstantPotentialForceFvec : public CpuConstantPotentialForce {
public:
    CpuConstantPotentialForceFvec();

protected:
    /**
     * Computes the direct space contribution to energies and forces for one
     * block of particles.
     * 
     * @param blockIndex  the index of the block
     * @param forces      output forces
     * @param energy      output energy
     */
    void getEnergyForcesBlock(int blockIndex, float* forces, double* energy);
    /**
     * Computes the direct space contribution to charge derivatives for one
     * block of particles.
     * 
     * @param blockIndex   the index of the block
     * @param derivatives  output derivatives
     */
    void getDerivativesBlock(int blockIndex, float* derivatives);
    /**
     * Compute an approximation of a function using a table lookup, where an
     * offset into the table can be provided for each value
     * 
     * @param table         start of tabulated values
     * @param x             values at which to evaluate the tabulated function
     * @param tableOffsets  offsets into the tabulated value array for each value
     */
    FVEC tableLookup(float const * table, FVEC x, IVEC tableOffsets);

private:
    static constexpr int blockSize = sizeof(FVEC) / sizeof(float);

    template<PeriodicType PERIODIC_TYPE>
    void getEnergyForcesBlockImpl(int blockIndex, float* forces, double* energy, fvec4& blockCenter);

    template<PeriodicType PERIODIC_TYPE>
    void getDerivativesBlockImpl(int blockIndex, float* derivatives, fvec4& blockCenter);

    void getBlockPeriodicType(int blockIndex, PeriodicType& periodicType, fvec4& blockCenter);

    template<PeriodicType PERIODIC_TYPE>
    void getDeltaR(const fvec4& posI, const FVEC& x, const FVEC& y, const FVEC& z, FVEC& dx, FVEC& dy, FVEC& dz, FVEC& r2) const;
};

template<typename FVEC, typename IVEC>
CpuConstantPotentialForceFvec<FVEC, IVEC>::CpuConstantPotentialForceFvec() : CpuConstantPotentialForce() {
}

template<typename FVEC, typename IVEC>
void CpuConstantPotentialForceFvec<FVEC, IVEC>::getEnergyForcesBlock(int blockIndex, float* forces, double* energy) {
    PeriodicType periodicType;
    fvec4 blockCenter;
    getBlockPeriodicType(blockIndex, periodicType, blockCenter);

    if (periodicType == NoPeriodic) {
        getEnergyForcesBlockImpl<NoPeriodic>(blockIndex, forces, energy, blockCenter);
    }
    else if (periodicType == PeriodicPerAtom) {
        getEnergyForcesBlockImpl<PeriodicPerAtom>(blockIndex, forces, energy, blockCenter);
    }
    else if (periodicType == PeriodicPerInteraction) {
        getEnergyForcesBlockImpl<PeriodicPerInteraction>(blockIndex, forces, energy, blockCenter);
    }
    else if (periodicType == PeriodicTriclinic) {
        getEnergyForcesBlockImpl<PeriodicTriclinic>(blockIndex, forces, energy, blockCenter);
    }
}

template<typename FVEC, typename IVEC>
void CpuConstantPotentialForceFvec<FVEC, IVEC>::getDerivativesBlock(int blockIndex, float* derivatives) {
    PeriodicType periodicType;
    fvec4 blockCenter;
    getBlockPeriodicType(blockIndex, periodicType, blockCenter);

    if (periodicType == NoPeriodic) {
        getDerivativesBlockImpl<NoPeriodic>(blockIndex, derivatives, blockCenter);
    }
    else if (periodicType == PeriodicPerAtom) {
        getDerivativesBlockImpl<PeriodicPerAtom>(blockIndex, derivatives, blockCenter);
    }
    else if (periodicType == PeriodicPerInteraction) {
        getDerivativesBlockImpl<PeriodicPerInteraction>(blockIndex, derivatives, blockCenter);
    }
    else if (periodicType == PeriodicTriclinic) {
        getDerivativesBlockImpl<PeriodicTriclinic>(blockIndex, derivatives, blockCenter);
    }
}

template<typename FVEC, typename IVEC>
FVEC CpuConstantPotentialForceFvec<FVEC, IVEC>::tableLookup(float const * table, FVEC x, IVEC tableOffsets) {
    FVEC offset = x * FVEC(tableScale);
    FVEC index = min(floor(offset), (float)NUM_TABLE_POINTS);

    FVEC s1, s2;
    gatherVecPair(table, index + tableOffsets, s1, s2);
    FVEC coeff = offset - index;
    return (s1 - coeff * s1) + coeff * s2;
}

template<typename FVEC, typename IVEC>
template<PeriodicType PERIODIC_TYPE>
void CpuConstantPotentialForceFvec<FVEC, IVEC>::getEnergyForcesBlockImpl(int blockIndex, float* forces, double* energy, fvec4& blockCenter) {
    const int32_t* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
    fvec4 blockAtomPosq[blockSize];
    FVEC blockAtomForceX(0.0f), blockAtomForceY(0.0f), blockAtomForceZ(0.0f);
    FVEC blockAtomX, blockAtomY, blockAtomZ, blockAtomCharge;
    for (int k = 0; k < blockSize; k++) {
        blockAtomPosq[k] = fvec4(posq + 4 * blockAtom[k]);
        if (PERIODIC_TYPE == PeriodicPerAtom) {
            blockAtomPosq[k] -= floor((blockAtomPosq[k] - blockCenter) * recipBoxSize + 0.5f) * boxSize;
        }
    }
    transpose(blockAtomPosq, blockAtomX, blockAtomY, blockAtomZ, blockAtomCharge);
    blockAtomCharge *= ONE_4PI_EPS0;

    FVEC cutoffSquared(nonbondedCutoff * nonbondedCutoff);

    // Get offsets into tables for block atoms.
    int tableOffsetsArray[blockSize];
    for (int k = 0; k < blockSize; k++) {
        tableOffsetsArray[k] = (sysElec[blockAtom[k]] + 1) * (numElectrodes + 1) * (NUM_TABLE_POINTS + 4);
    }
    IVEC tableOffsets(tableOffsetsArray);

    CpuNeighborList::NeighborIterator neighbors = neighborList->getNeighborIterator(blockIndex);
    FVEC neighborsEnergy(0.0f);
    while (neighbors.next()) {
        int atom = neighbors.getNeighbor();

        FVEC dx, dy, dz, r2;
        fvec4 atomPos(posq + 4 * atom);
        if (PERIODIC_TYPE == PeriodicPerAtom) {
            atomPos -= floor((atomPos - blockCenter) * recipBoxSize + 0.5f) * boxSize;
        }
        getDeltaR<PERIODIC_TYPE>(atomPos, blockAtomX, blockAtomY, blockAtomZ, dx, dy, dz, r2);

        auto include = FVEC::expandBitsToMask(~neighbors.getExclusions());
        include = blendZero(r2 < cutoffSquared, include);
        if (!any(include)) {
            continue;
        }

        // Get offset into tables for single atom.
        int jOffset = (sysElec[atom] + 1) * (NUM_TABLE_POINTS + 4);

        FVEC inverseR = rsqrt(r2);
        FVEC r = r2 * inverseR;

        FVEC qqFactor = blockAtomCharge * posq[4 * atom + 3];
        FVEC qqInverseR = qqFactor * inverseR;
        FVEC forceFactor = qqInverseR * inverseR * inverseR * tableLookup(&forceLookupTable[jOffset], r, tableOffsets);

        // Accumulate forces for block atoms.
        forceFactor = blendZero(forceFactor, include);
        FVEC forceX = dx * forceFactor;
        FVEC forceY = dy * forceFactor;
        FVEC forceZ = dz * forceFactor;
        blockAtomForceX += forceX;
        blockAtomForceY += forceY;
        blockAtomForceZ += forceZ;

        // Apply forces to single atom.
        float* atomForce = forces + 4 * atom;
        fvec4 newAtomForce = fvec4(atomForce) - reduceToVec3(forceX, forceY, forceZ);
        newAtomForce.store(atomForce);

        if (energy) {
            neighborsEnergy += blendZero(qqInverseR * tableLookup(&energyLookupTable[jOffset], r, tableOffsets), include);
        }
    }

    if (energy) {
        *energy += reduceAdd(neighborsEnergy);
    }

    // Apply forces to block atoms.
    fvec4 f[blockSize];
    transpose(blockAtomForceX, blockAtomForceY, blockAtomForceZ, 0.0f, f);
    for (int k = 0; k < blockSize; k++) {
        (fvec4(forces + 4 * blockAtom[k]) + f[k]).store(forces + 4 * blockAtom[k]);
    }
}

template<typename FVEC, typename IVEC>
template<PeriodicType PERIODIC_TYPE>
void CpuConstantPotentialForceFvec<FVEC, IVEC>::getDerivativesBlockImpl(int blockIndex, float* derivatives, fvec4& blockCenter) {
    const int32_t* blockAtom = &neighborList->getSortedAtoms()[blockSize*blockIndex];
    fvec4 blockAtomPosq[blockSize];
    FVEC blockAtomX, blockAtomY, blockAtomZ, blockAtomCharge;
    FVEC blockAtomDerivatives(0.0f);
    for (int k = 0; k < blockSize; k++) {
        blockAtomPosq[k] = fvec4(posq + 4 * blockAtom[k]);
        if (PERIODIC_TYPE == PeriodicPerAtom) {
            blockAtomPosq[k] -= floor((blockAtomPosq[k] - blockCenter) * recipBoxSize + 0.5f) * boxSize;
        }
    }
    transpose(blockAtomPosq, blockAtomX, blockAtomY, blockAtomZ, blockAtomCharge);
    FVEC cutoffSquared(nonbondedCutoff * nonbondedCutoff);

    // Get offsets into tables for block atoms.
    int iiArray[blockSize];
    int tableOffsetsArray[blockSize];
    for (int k = 0; k < blockSize; k++) {
        iiArray[k] = sysToElec[blockAtom[k]];
        tableOffsetsArray[k] = (sysElec[blockAtom[k]] + 1) * (numElectrodes + 1) * (NUM_TABLE_POINTS + 4);
    }
    IVEC ii(iiArray);
    IVEC tableOffsets(tableOffsetsArray);

    // Check if any block atoms are electrode atoms.
    IVEC iiInclude = ii != IVEC(-1);
    bool iiIncludeAny = any(iiInclude);

    CpuNeighborList::NeighborIterator neighbors = neighborList->getNeighborIterator(blockIndex);
    while (neighbors.next()) {
        int atom = neighbors.getNeighbor();

        FVEC dx, dy, dz, r2;
        fvec4 atomPos(posq + 4 * atom);
        if (PERIODIC_TYPE == PeriodicPerAtom) {
            atomPos -= floor((atomPos - blockCenter) * recipBoxSize + 0.5f) * boxSize;
        }
        getDeltaR<PERIODIC_TYPE>(atomPos, blockAtomX, blockAtomY, blockAtomZ, dx, dy, dz, r2);

        auto include = FVEC::expandBitsToMask(~neighbors.getExclusions());
        include = blendZero(r2 < cutoffSquared, include);
        if (!any(include)) {
            continue;
        }

        // Derivatives must be evaluated if any block atoms are electrode atoms
        // OR the single atom is an electrode atom.
        int jj = sysToElec[atom];
        bool jjInclude = jj != -1;
        if (!(iiIncludeAny || jjInclude)) {
            continue;
        }

        // Get offset into tables for single atom.
        int jOffset = (sysElec[atom] + 1) * (NUM_TABLE_POINTS + 4);

        FVEC inverseR = rsqrt(r2);
        FVEC r = r2 * inverseR;
        FVEC dq = blendZero(ONE_4PI_EPS0 * inverseR * tableLookup(&energyLookupTable[jOffset], r, tableOffsets), include);

        // Accumulate derivatives for block atoms.
        blockAtomDerivatives += blendZero(posq[4 * atom + 3] * dq, iiInclude);

        // Apply derivatives to single atom.
        if (jjInclude) {
            derivatives[4 * jj] += reduceAdd(blockAtomCharge * dq);
        }
    }

    // Apply derivatives to block atoms.
    float blockAtomDerivativesArray[blockSize];
    blockAtomDerivatives.store(blockAtomDerivativesArray);
    for (int k = 0; k < blockSize; k++) {
        derivatives[4 * iiArray[k]] += blockAtomDerivativesArray[k];
    }
}

template<typename FVEC, typename IVEC>
void CpuConstantPotentialForceFvec<FVEC, IVEC>::getBlockPeriodicType(int blockIndex, PeriodicType& periodicType, fvec4& blockCenter) {
    const int32_t* blockAtom = &neighborList->getSortedAtoms()[blockSize * blockIndex];
    float minx, maxx, miny, maxy, minz, maxz;
    minx = maxx = posq[4 * blockAtom[0]];
    miny = maxy = posq[4 * blockAtom[0] + 1];
    minz = maxz = posq[4 * blockAtom[0] + 2];
    for (int i = 1; i < blockSize; i++) {
        minx = std::min(minx, posq[4 * blockAtom[i]]);
        maxx = std::max(maxx, posq[4 * blockAtom[i]]);
        miny = std::min(miny, posq[4 * blockAtom[i] + 1]);
        maxy = std::max(maxy, posq[4 * blockAtom[i] + 1]);
        minz = std::min(minz, posq[4 * blockAtom[i] + 2]);
        maxz = std::max(maxz, posq[4 * blockAtom[i] + 2]);
    }
    blockCenter = fvec4(0.5f * (minx + maxx), 0.5f * (miny + maxy), 0.5f * (minz + maxz), 0.0f);
    if (!(minx < nonbondedCutoff || miny < nonbondedCutoff || minz < nonbondedCutoff ||
            maxx > boxSize[0] - nonbondedCutoff || maxy > boxSize[1] - nonbondedCutoff || maxz > boxSize[2] - nonbondedCutoff)) {
        periodicType = NoPeriodic;
    }
    else if (triclinic) {
        periodicType = PeriodicTriclinic;
    }
    else if (0.5f * (boxSize[0] - (maxx - minx)) >= nonbondedCutoff &&
             0.5f * (boxSize[1] - (maxy - miny)) >= nonbondedCutoff &&
             0.5f * (boxSize[2] - (maxz - minz)) >= nonbondedCutoff) {
        periodicType = PeriodicPerAtom;
    }
    else {
        periodicType = PeriodicPerInteraction;
    }
}

template<typename FVEC, typename IVEC>
template<PeriodicType PERIODIC_TYPE>
void CpuConstantPotentialForceFvec<FVEC, IVEC>::getDeltaR(const fvec4& posI, const FVEC& x, const FVEC& y, const FVEC& z, FVEC& dx, FVEC& dy, FVEC& dz, FVEC& r2) const {
    dx = x - posI[0];
    dy = y - posI[1];
    dz = z - posI[2];
    if (PERIODIC_TYPE == PeriodicTriclinic) {
        const auto scale3 = floor(dz * recipBoxSize[2] + 0.5f);
        dx -= scale3 * boxVectorsVec4[2][0];
        dy -= scale3 * boxVectorsVec4[2][1];
        dz -= scale3 * boxVectorsVec4[2][2];
        const auto scale2 = floor(dy * recipBoxSize[1] + 0.5f);
        dx -= scale2 * boxVectorsVec4[1][0];
        dy -= scale2 * boxVectorsVec4[1][1];
        const auto scale1 = floor(dx * recipBoxSize[0] + 0.5f);
        dx -= scale1 * boxVectorsVec4[0][0];
    }
    else if (PERIODIC_TYPE == PeriodicPerInteraction) {
        dx -= round(dx * recipBoxSize[0]) * boxSize[0];
        dy -= round(dy * recipBoxSize[1]) * boxSize[1];
        dz -= round(dz * recipBoxSize[2]) * boxSize[2];
    }
    r2 = dx * dx + dy * dy + dz * dz;
}

} // namespace OpenMM

#endif // OPENMM_CPUCONSTANTPOTENTIALFORCEFVEC_H_
