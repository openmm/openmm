
/* Portions copyright (c) 2006-2015 Stanford University and Simbios.
 * Contributors: Pande Group
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

#include "SimTKOpenMMUtilities.h"
#include "CpuNonbondedForceVec4.h"
#include <algorithm>
#include <iostream>

using namespace std;
using namespace OpenMM;

/**
 * Factory method to create a CpuNonbondedForceVec4.
 */
CpuNonbondedForce* createCpuNonbondedForceVec4() {
    return new CpuNonbondedForceVec4();
}

/**---------------------------------------------------------------------------------------

   CpuNonbondedForceVec4 constructor

   --------------------------------------------------------------------------------------- */

CpuNonbondedForceVec4::CpuNonbondedForceVec4() {
}

enum PeriodicType {NoPeriodic, PeriodicPerAtom, PeriodicPerInteraction, PeriodicTriclinic};

void CpuNonbondedForceVec4::calculateBlockIxn(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Determine whether we need to apply periodic boundary conditions.
    
    PeriodicType periodicType;
    fvec4 blockCenter;
    if (!periodic) {
        periodicType = NoPeriodic;
        blockCenter = 0.0f;
    }
    else {
        const int* blockAtom = &neighborList->getSortedAtoms()[4*blockIndex];
        float minx, maxx, miny, maxy, minz, maxz;
        minx = maxx = posq[4*blockAtom[0]];
        miny = maxy = posq[4*blockAtom[0]+1];
        minz = maxz = posq[4*blockAtom[0]+2];
        for (int i = 1; i < 4; i++) {
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
    
    if (periodicType == NoPeriodic)
        calculateBlockIxnImpl<NoPeriodic>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicPerAtom)
        calculateBlockIxnImpl<PeriodicPerAtom>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicPerInteraction)
        calculateBlockIxnImpl<PeriodicPerInteraction>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicTriclinic)
        calculateBlockIxnImpl<PeriodicTriclinic>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
}

template <int PERIODIC_TYPE>
void CpuNonbondedForceVec4::calculateBlockIxnImpl(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize, const fvec4& blockCenter) {
    // Load the positions and parameters of the atoms in the block.
    
    const int* blockAtom = &neighborList->getSortedAtoms()[4*blockIndex];
    fvec4 blockAtomPosq[4];
    fvec4 blockAtomForceX(0.0f), blockAtomForceY(0.0f), blockAtomForceZ(0.0f);
    for (int i = 0; i < 4; i++) {
        blockAtomPosq[i] = fvec4(posq+4*blockAtom[i]);
        if (PERIODIC_TYPE == PeriodicPerAtom)
            blockAtomPosq[i] -= floor((blockAtomPosq[i]-blockCenter)*invBoxSize+0.5f)*boxSize;
    }
    fvec4 blockAtomX = fvec4(blockAtomPosq[0][0], blockAtomPosq[1][0], blockAtomPosq[2][0], blockAtomPosq[3][0]);
    fvec4 blockAtomY = fvec4(blockAtomPosq[0][1], blockAtomPosq[1][1], blockAtomPosq[2][1], blockAtomPosq[3][1]);
    fvec4 blockAtomZ = fvec4(blockAtomPosq[0][2], blockAtomPosq[1][2], blockAtomPosq[2][2], blockAtomPosq[3][2]);
    fvec4 blockAtomCharge = fvec4(ONE_4PI_EPS0)*fvec4(blockAtomPosq[0][3], blockAtomPosq[1][3], blockAtomPosq[2][3], blockAtomPosq[3][3]);
    fvec4 blockAtomSigma(atomParameters[blockAtom[0]].first, atomParameters[blockAtom[1]].first, atomParameters[blockAtom[2]].first, atomParameters[blockAtom[3]].first);
    fvec4 blockAtomEpsilon(atomParameters[blockAtom[0]].second, atomParameters[blockAtom[1]].second, atomParameters[blockAtom[2]].second, atomParameters[blockAtom[3]].second);
    const bool needPeriodic = (PERIODIC_TYPE == PeriodicPerInteraction || PERIODIC_TYPE == PeriodicTriclinic);
    const float invSwitchingInterval = 1/(cutoffDistance-switchingDistance);
    
    // Loop over neighbors for this block.
    
    const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
    const vector<char>& exclusions = neighborList->getBlockExclusions(blockIndex);
    for (int i = 0; i < (int) neighbors.size(); i++) {
        // Load the next neighbor.
        
        int atom = neighbors[i];
        
        // Compute the distances to the block atoms.
        
        fvec4 dx, dy, dz, r2;
        fvec4 atomPos(posq+4*atom);
        if (PERIODIC_TYPE == PeriodicPerAtom)
            atomPos -= floor((atomPos-blockCenter)*invBoxSize+0.5f)*boxSize;
        getDeltaR<PERIODIC_TYPE>(atomPos, blockAtomX, blockAtomY, blockAtomZ, dx, dy, dz, r2, needPeriodic, boxSize, invBoxSize);
        ivec4 include;
        char excl = exclusions[i];
        if (excl == 0)
            include = -1;
        else
            include = ivec4(excl&1 ? 0 : -1, excl&2 ? 0 : -1, excl&4 ? 0 : -1, excl&8 ? 0 : -1);
        include = include & (r2 < cutoffDistance*cutoffDistance);
        if (!any(include))
            continue; // No interactions to compute.
        
        // Compute the interactions.
        
        fvec4 inverseR = rsqrt(r2);
        fvec4 energy, dEdR;
        float atomEpsilon = atomParameters[atom].second;
        if (atomEpsilon != 0.0f) {
            fvec4 sig = blockAtomSigma+atomParameters[atom].first;
            fvec4 sig2 = inverseR*sig;
            sig2 *= sig2;
            fvec4 sig6 = sig2*sig2*sig2;
            fvec4 epsSig6 = blockAtomEpsilon*atomEpsilon*sig6;
            dEdR = epsSig6*(12.0f*sig6 - 6.0f);
            energy = epsSig6*(sig6-1.0f);
            if (useSwitch) {
                fvec4 r = r2*inverseR;
                fvec4 t = blend(0.0f, (r-switchingDistance)*invSwitchingInterval, r>switchingDistance);
                fvec4 switchValue = 1+t*t*t*(-10.0f+t*(15.0f-t*6.0f));
                fvec4 switchDeriv = t*t*(-30.0f+t*(60.0f-t*30.0f))*invSwitchingInterval;
                dEdR = switchValue*dEdR - energy*switchDeriv*r;
                energy *= switchValue;
            }
        }
        else {
            energy = 0.0f;
            dEdR = 0.0f;
        }
        fvec4 chargeProd = blockAtomCharge*posq[4*atom+3];
        if (cutoff)
            dEdR += chargeProd*(inverseR-2.0f*krf*r2);
        else
            dEdR += chargeProd*inverseR;
        dEdR *= inverseR*inverseR;

        // Accumulate energies.

        fvec4 one(1.0f);
        if (totalEnergy) {
            if (cutoff)
                energy += chargeProd*(inverseR+krf*r2-crf);
            else
                energy += chargeProd*inverseR;
            energy = blend(0.0f, energy, include);
            *totalEnergy += dot4(energy, one);
        }

        // Accumulate forces.

        dEdR = blend(0.0f, dEdR, include);
        fvec4 fx = dx*dEdR;
        fvec4 fy = dy*dEdR;
        fvec4 fz = dz*dEdR;
        blockAtomForceX += fx;
        blockAtomForceY += fy;
        blockAtomForceZ += fz;
        float* atomForce = forces+4*atom;
        atomForce[0] -= dot4(fx, one);
        atomForce[1] -= dot4(fy, one);
        atomForce[2] -= dot4(fz, one);
    }
    
    // Record the forces on the block atoms.

    fvec4 f[4] = {blockAtomForceX, blockAtomForceY, blockAtomForceZ, 0.0f};
    transpose(f[0], f[1], f[2], f[3]);
    for (int j = 0; j < 4; j++)
        (fvec4(forces+4*blockAtom[j])+f[j]).store(forces+4*blockAtom[j]);
  }

void CpuNonbondedForceVec4::calculateBlockEwaldIxn(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Determine whether we need to apply periodic boundary conditions.
    PeriodicType periodicType;
    fvec4 blockCenter;
    if (!periodic) {
        periodicType = NoPeriodic;
        blockCenter = 0.0f;
    }
    else {
        const int* blockAtom = &neighborList->getSortedAtoms()[4*blockIndex];
        float minx, maxx, miny, maxy, minz, maxz;
        minx = maxx = posq[4*blockAtom[0]];
        miny = maxy = posq[4*blockAtom[0]+1];
        minz = maxz = posq[4*blockAtom[0]+2];
        for (int i = 1; i < 4; i++) {
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
    
    if (periodicType == NoPeriodic)
        calculateBlockEwaldIxnImpl<NoPeriodic>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicPerAtom)
        calculateBlockEwaldIxnImpl<PeriodicPerAtom>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicPerInteraction)
        calculateBlockEwaldIxnImpl<PeriodicPerInteraction>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
    else if (periodicType == PeriodicTriclinic)
        calculateBlockEwaldIxnImpl<PeriodicTriclinic>(blockIndex, forces, totalEnergy, boxSize, invBoxSize, blockCenter);
}

template <int PERIODIC_TYPE>
void CpuNonbondedForceVec4::calculateBlockEwaldIxnImpl(int blockIndex, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize, const fvec4& blockCenter) {
    // Load the positions and parameters of the atoms in the block.
    const int* blockAtom = &neighborList->getSortedAtoms()[4*blockIndex];
    fvec4 blockAtomPosq[4];
    fvec4 blockAtomForceX(0.0f), blockAtomForceY(0.0f), blockAtomForceZ(0.0f);
    for (int i = 0; i < 4; i++) {
        blockAtomPosq[i] = fvec4(posq+4*blockAtom[i]);
        if (PERIODIC_TYPE == PeriodicPerAtom)
            blockAtomPosq[i] -= floor((blockAtomPosq[i]-blockCenter)*invBoxSize+0.5f)*boxSize;
    }
    fvec4 blockAtomX = fvec4(blockAtomPosq[0][0], blockAtomPosq[1][0], blockAtomPosq[2][0], blockAtomPosq[3][0]);
    fvec4 blockAtomY = fvec4(blockAtomPosq[0][1], blockAtomPosq[1][1], blockAtomPosq[2][1], blockAtomPosq[3][1]);
    fvec4 blockAtomZ = fvec4(blockAtomPosq[0][2], blockAtomPosq[1][2], blockAtomPosq[2][2], blockAtomPosq[3][2]);
    fvec4 blockAtomCharge = fvec4(ONE_4PI_EPS0)*fvec4(blockAtomPosq[0][3], blockAtomPosq[1][3], blockAtomPosq[2][3], blockAtomPosq[3][3]);
    fvec4 blockAtomSigma(atomParameters[blockAtom[0]].first, atomParameters[blockAtom[1]].first, atomParameters[blockAtom[2]].first, atomParameters[blockAtom[3]].first);
    fvec4 blockAtomEpsilon(atomParameters[blockAtom[0]].second, atomParameters[blockAtom[1]].second, atomParameters[blockAtom[2]].second, atomParameters[blockAtom[3]].second);
    fvec4 C6s(C6params[blockAtom[0]], C6params[blockAtom[1]], C6params[blockAtom[2]], C6params[blockAtom[3]]);
    const bool needPeriodic = (PERIODIC_TYPE == PeriodicPerInteraction || PERIODIC_TYPE == PeriodicTriclinic);
    const float invSwitchingInterval = 1/(cutoffDistance-switchingDistance);

    // Loop over neighbors for this block.
    
    const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
    const vector<char>& exclusions = neighborList->getBlockExclusions(blockIndex);
    for (int i = 0; i < (int) neighbors.size(); i++) {
        // Load the next neighbor.
        
        int atom = neighbors[i];
        
        // Compute the distances to the block atoms.
        
        fvec4 dx, dy, dz, r2;
        fvec4 atomPos(posq+4*atom);
        if (PERIODIC_TYPE == PeriodicPerAtom)
            atomPos -= floor((atomPos-blockCenter)*invBoxSize+0.5f)*boxSize;
        getDeltaR<PERIODIC_TYPE>(atomPos, blockAtomX, blockAtomY, blockAtomZ, dx, dy, dz, r2, needPeriodic, boxSize, invBoxSize);
        ivec4 include;
        char excl = exclusions[i];
        if (excl == 0)
            include = -1;
        else
            include = ivec4(excl&1 ? 0 : -1, excl&2 ? 0 : -1, excl&4 ? 0 : -1, excl&8 ? 0 : -1);
        include = include & (r2 < cutoffDistance*cutoffDistance);
        if (!any(include))
            continue; // No interactions to compute.
        
        // Compute the interactions.
        
        fvec4 inverseR = rsqrt(r2);
        fvec4 r = r2*inverseR;
        fvec4 energy, dEdR;
        float atomEpsilon = atomParameters[atom].second;
        if (atomEpsilon != 0.0f) {
            fvec4 sig = blockAtomSigma+atomParameters[atom].first;
            fvec4 sig2 = inverseR*sig;
            sig2 *= sig2;
            fvec4 sig6 = sig2*sig2*sig2;
            fvec4 eps = blockAtomEpsilon*atomEpsilon;
            fvec4 epsSig6 = eps*sig6;
            dEdR = epsSig6*(12.0f*sig6 - 6.0f);
            energy = epsSig6*(sig6-1.0f);
            if (useSwitch) {
                fvec4 t = blend(0.0f, (r-switchingDistance)*invSwitchingInterval, r>switchingDistance);
                fvec4 switchValue = 1+t*t*t*(-10.0f+t*(15.0f-t*6.0f));
                fvec4 switchDeriv = t*t*(-30.0f+t*(60.0f-t*30.0f))*invSwitchingInterval;
                dEdR = switchValue*dEdR - energy*switchDeriv*r;
                energy *= switchValue;
            }

            if (ljpme) {
                fvec4 C6ij = C6s*C6params[atom];
                fvec4 inverseR2 = inverseR*inverseR;
                fvec4 mysig2 = sig*sig;
                fvec4 mysig6 = mysig2*mysig2*mysig2;
                fvec4 emult = C6ij*inverseR2*inverseR2*inverseR2*exptermsApprox(r);
                fvec4 potentialShift = eps*(1.0f-mysig6*inverseRcut6)*mysig6*inverseRcut6 - C6ij*inverseRcut6Expterm;
                dEdR += 6.0f*C6ij*inverseR2*inverseR2*inverseR2*dExptermsApprox(r);
                energy += emult + potentialShift;
            }
        }
        else {
            energy = 0.0f;
            dEdR = 0.0f;
        }
        fvec4 chargeProd = blockAtomCharge*posq[4*atom+3];
        dEdR += chargeProd*inverseR*ewaldScaleFunction(r);
        dEdR *= inverseR*inverseR;        

        // Accumulate energies.

        fvec4 one(1.0f);
        if (totalEnergy) {
            energy += chargeProd*inverseR*erfcApprox(alphaEwald*r);
            energy = blend(0.0f, energy, include);
            *totalEnergy += dot4(energy, one);
        }

        // Accumulate forces.

        dEdR = blend(0.0f, dEdR, include);
        fvec4 fx = dx*dEdR;
        fvec4 fy = dy*dEdR;
        fvec4 fz = dz*dEdR;
        blockAtomForceX += fx;
        blockAtomForceY += fy;
        blockAtomForceZ += fz;
        float* atomForce = forces+4*atom;
        atomForce[0] -= dot4(fx, one);
        atomForce[1] -= dot4(fy, one);
        atomForce[2] -= dot4(fz, one);
    }
    
    // Record the forces on the block atoms.

    fvec4 f[4] = {blockAtomForceX, blockAtomForceY, blockAtomForceZ, 0.0f};
    transpose(f[0], f[1], f[2], f[3]);
    for (int j = 0; j < 4; j++)
        (fvec4(forces+4*blockAtom[j])+f[j]).store(forces+4*blockAtom[j]);
}

template <int PERIODIC_TYPE>
void CpuNonbondedForceVec4::getDeltaR(const fvec4& posI, const fvec4& x, const fvec4& y, const fvec4& z, fvec4& dx, fvec4& dy, fvec4& dz, fvec4& r2, bool periodic, const fvec4& boxSize, const fvec4& invBoxSize) const {
    dx = x-posI[0];
    dy = y-posI[1];
    dz = z-posI[2];
    if (PERIODIC_TYPE == PeriodicTriclinic) {
        fvec4 scale3 = floor(dz*recipBoxSize[2]+0.5f);
        dx -= scale3*periodicBoxVectors[2][0];
        dy -= scale3*periodicBoxVectors[2][1];
        dz -= scale3*periodicBoxVectors[2][2];
        fvec4 scale2 = floor(dy*recipBoxSize[1]+0.5f);
        dx -= scale2*periodicBoxVectors[1][0];
        dy -= scale2*periodicBoxVectors[1][1];
        fvec4 scale1 = floor(dx*recipBoxSize[0]+0.5f);
        dx -= scale1*periodicBoxVectors[0][0];
    }
    else if (PERIODIC_TYPE == PeriodicPerInteraction) {
        dx -= round(dx*invBoxSize[0])*boxSize[0];
        dy -= round(dy*invBoxSize[1])*boxSize[1];
        dz -= round(dz*invBoxSize[2])*boxSize[2];
    }
    r2 = dx*dx + dy*dy + dz*dz;
}

fvec4 CpuNonbondedForceVec4::erfcApprox(const fvec4& x) {
    fvec4 x1 = x*erfcDXInv;
    ivec4 index = min(floor(x1), NUM_TABLE_POINTS);
    fvec4 coeff2 = x1-index;
    fvec4 coeff1 = 1.0f-coeff2;
    fvec4 t1(&erfcTable[index[0]]);
    fvec4 t2(&erfcTable[index[1]]);
    fvec4 t3(&erfcTable[index[2]]);
    fvec4 t4(&erfcTable[index[3]]);
    transpose(t1, t2, t3, t4);
    return coeff1*t1 + coeff2*t2;
}

fvec4 CpuNonbondedForceVec4::ewaldScaleFunction(const fvec4& x) {
    // Compute the tabulated Ewald scale factor: erfc(alpha*r) + 2*alpha*r*exp(-alpha*alpha*r*r)/sqrt(PI)

    fvec4 x1 = x*ewaldDXInv;
    ivec4 index = min(floor(x1), NUM_TABLE_POINTS);
    fvec4 coeff2 = x1-index;
    fvec4 coeff1 = 1.0f-coeff2;
    fvec4 t1(&ewaldScaleTable[index[0]]);
    fvec4 t2(&ewaldScaleTable[index[1]]);
    fvec4 t3(&ewaldScaleTable[index[2]]);
    fvec4 t4(&ewaldScaleTable[index[3]]);
    transpose(t1, t2, t3, t4);
    return coeff1*t1 + coeff2*t2;
}

fvec4 CpuNonbondedForceVec4::exptermsApprox(const fvec4& r) {
    fvec4 r1 = r*exptermsDXInv;
    ivec4 index = min(floor(r1), NUM_TABLE_POINTS);
    fvec4 coeff2 = r1-index;
    fvec4 coeff1 = 1.0f-coeff2;
    fvec4 t1(&exptermsTable[index[0]]);
    fvec4 t2(&exptermsTable[index[1]]);
    fvec4 t3(&exptermsTable[index[2]]);
    fvec4 t4(&exptermsTable[index[3]]);
    transpose(t1, t2, t3, t4);
    return coeff1*t1 + coeff2*t2;
}

fvec4 CpuNonbondedForceVec4::dExptermsApprox(const fvec4& r) {
    fvec4 r1 = r*exptermsDXInv;
    ivec4 index = min(floor(r1), NUM_TABLE_POINTS);
    fvec4 coeff2 = r1-index;
    fvec4 coeff1 = 1.0f-coeff2;
    fvec4 t1(&dExptermsTable[index[0]]);
    fvec4 t2(&dExptermsTable[index[1]]);
    fvec4 t3(&dExptermsTable[index[2]]);
    fvec4 t4(&dExptermsTable[index[3]]);
    transpose(t1, t2, t3, t4);
    return coeff1*t1 + coeff2*t2;
}

