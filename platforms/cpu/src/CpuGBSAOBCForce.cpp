/* Portions copyright (c) 2006-2018 Stanford University and Simbios.
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

#include "CpuGBSAOBCForce.h"
#include "SimTKOpenMMRealType.h"
#include "openmm/internal/vectorize.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace OpenMM;

const int CpuGBSAOBCForce::NUM_TABLE_POINTS = 4096;
const float CpuGBSAOBCForce::TABLE_MIN = 0.25f;
const float CpuGBSAOBCForce::TABLE_MAX = 1.5f;

CpuGBSAOBCForce::CpuGBSAOBCForce() : cutoff(false), periodic(false) {
    logDX = (TABLE_MAX-TABLE_MIN)/NUM_TABLE_POINTS;
    logDXInv = 1.0f/logDX;
    logTable.resize(NUM_TABLE_POINTS+4);
    for (int i = 0; i < NUM_TABLE_POINTS+4; i++) {
        double x = TABLE_MIN+i*logDX;
        logTable[i] = log(x);
    }
}

void CpuGBSAOBCForce::setUseCutoff(float distance) {
    cutoff = true;
    cutoffDistance = distance;
}

void CpuGBSAOBCForce::setPeriodic(float* periodicBoxSize) {
    periodic = true;
    this->periodicBoxSize[0] = periodicBoxSize[0];
    this->periodicBoxSize[1] = periodicBoxSize[1];
    this->periodicBoxSize[2] = periodicBoxSize[2];
}

void CpuGBSAOBCForce::setSoluteDielectric(float dielectric) {
    soluteDielectric = dielectric;
}

void CpuGBSAOBCForce::setSolventDielectric(float dielectric) {
    solventDielectric = dielectric;
}

void CpuGBSAOBCForce::setSurfaceAreaEnergy(float energy) {
    surfaceAreaFactor = 4*M_PI*energy;
}

const std::vector<std::pair<float, float> >& CpuGBSAOBCForce::getParticleParameters() const {
    return particleParams;
}

void CpuGBSAOBCForce::setParticleParameters(const std::vector<std::pair<float, float> >& params) {
    particleParams = params;
    bornRadii.resize(params.size()+3);
    obcChain.resize(params.size()+3);
    for (int i = bornRadii.size()-3; i < bornRadii.size(); i++) {
        bornRadii[i] = 0;
        obcChain[i] = 0;
    }
}

void CpuGBSAOBCForce::computeForce(const AlignedArray<float>& posq, vector<AlignedArray<float> >& threadForce, double* totalEnergy, ThreadPool& threads) {
    // Record the parameters for the threads.
    
    this->posq = &posq[0];
    this->threadForce = &threadForce;
    includeEnergy = (totalEnergy != NULL);
    int numThreads = threads.getNumThreads();
    threadEnergy.resize(numThreads);
    threadBornForces.resize(numThreads);
    for (int i = 0; i < numThreads; i++)
        threadBornForces[i].resize(particleParams.size()+3);
    
    // Signal the threads to start running and wait for them to finish.
    
    atomicCounter = 0;
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadComputeForce(threads, threadIndex); });
    threads.waitForThreads(); // Compute Born radii
    atomicCounter = 0;
    threads.resumeThreads();
    threads.waitForThreads(); // Compute surface area term
    atomicCounter = 0;
    threads.resumeThreads();
    threads.waitForThreads(); // First loop
    atomicCounter = 0;
    threads.resumeThreads();
    threads.waitForThreads(); // Second loop
    
    // Combine the energies from all the threads.
    
    if (totalEnergy != NULL) {
        double energy = 0;
        for (int i = 0; i < numThreads; i++)
            energy += threadEnergy[i];
        *totalEnergy += (float) energy;
    }
}

void CpuGBSAOBCForce::threadComputeForce(ThreadPool& threads, int threadIndex) {
    int numParticles = particleParams.size();
    int numThreads = threads.getNumThreads();
    const float dielectricOffset = 0.009;
    const float alphaObc = 1.0f;
    const float betaObc = 0.8f;
    const float gammaObc = 4.85f;
    fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
    fvec4 invBoxSize((1/periodicBoxSize[0]), (1/periodicBoxSize[1]), (1/periodicBoxSize[2]), 0);
    int start = (threadIndex*numParticles)/numThreads;
    int end = ((threadIndex+1)*numParticles)/numThreads;

    // Calculate Born radii

    while (true) {
        int blockStart = atomicCounter.fetch_add(4);
        if (blockStart >= numParticles)
            break;
        int numInBlock = min(4, numParticles-blockStart);
        ivec4 blockAtomIndex(blockStart, blockStart+1, blockStart+2, blockStart+3);
        float atomRadius[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        float atomx[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        float atomy[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        float atomz[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        int blockMask[4] = {0, 0, 0, 0};
        for (int i = 0; i < numInBlock; i++) {
            int atomIndex = blockStart+i;
            atomRadius[i] = particleParams[atomIndex].first;
            atomx[i] = posq[4*atomIndex];
            atomy[i] = posq[4*atomIndex+1];
            atomz[i] = posq[4*atomIndex+2];
            blockMask[i] = 0xFFFFFFFF;
        }
        fvec4 offsetRadiusI(atomRadius);
        fvec4 radiusIInverse = 1.0f/offsetRadiusI;
        fvec4 x(atomx);
        fvec4 y(atomy);
        fvec4 z(atomz);
        ivec4 mask(blockMask);
        float sum[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        for (int atomJ = 0; atomJ < numParticles; atomJ++) {
            fvec4 posJ(posq+4*atomJ);
            fvec4 dx, dy, dz, r2;
            getDeltaR(posJ, x, y, z, dx, dy, dz, r2, periodic, boxSize, invBoxSize);
            ivec4 include = mask & (blockAtomIndex != ivec4(atomJ));
            if (cutoff)
                include = include & (r2 < cutoffDistance*cutoffDistance);
            if (!any(include))
                continue;
            fvec4 r = sqrt(r2);
            float scaledRadiusJ = particleParams[atomJ].second;
            float scaledRadiusJ2 = scaledRadiusJ*scaledRadiusJ;
            fvec4 rScaledRadiusJ = r + scaledRadiusJ;
            include = include & (offsetRadiusI < rScaledRadiusJ);
            fvec4 l_ij = 1.0f/max(offsetRadiusI, abs(r-scaledRadiusJ));
            fvec4 u_ij = 1.0f/rScaledRadiusJ;
            fvec4 l_ij2 = l_ij*l_ij;
            fvec4 u_ij2 = u_ij*u_ij;
            fvec4 rInverse = 1.0f/r;
            fvec4 r2Inverse = rInverse*rInverse;
            fvec4 logRatio = fastLog(u_ij/l_ij);
            fvec4 term = l_ij - u_ij + 0.25f*r*(u_ij2 - l_ij2) + (0.5f*rInverse*logRatio) + (0.25f*scaledRadiusJ*scaledRadiusJ*rInverse)*(l_ij2 - u_ij2);
            for (int j = 0; j < 4; j++) {
                if (include[j]) {
                    sum[j] += term[j];
                    if (offsetRadiusI[j] < scaledRadiusJ-r[j])
                        sum[j] += 2.0f*(radiusIInverse[j]-l_ij[j]);
                }
            }
        }
        for (int i = 0; i < numInBlock; i++) {
            int atomIndex = blockStart+i;
            sum[i] *= 0.5f*atomRadius[i];
            float sum2 = sum[i]*sum[i];
            float sum3 = sum[i]*sum2;
            float tanhSum = tanh(alphaObc*sum[i] - betaObc*sum2 + gammaObc*sum3);
            float radiusI = atomRadius[i] + dielectricOffset;
            bornRadii[atomIndex] = 1.0f/(1.0f/atomRadius[i] - tanhSum/radiusI);
            obcChain[atomIndex] = atomRadius[i]*(alphaObc - 2.0f*betaObc*sum[i] + 3.0f*gammaObc*sum2);
            obcChain[atomIndex] = (1.0f - tanhSum*tanhSum)*obcChain[atomIndex]/radiusI;
        }
    }
    threads.syncThreads();

    // Calculate ACE surface area term.

    const float probeRadius = 0.14f;
    double energy = 0.0;
    AlignedArray<float>& bornForces = threadBornForces[threadIndex];
    for (int i = 0; i < numParticles; i++)
        bornForces[i] = 0.0f;
    while (true) {
        int atomI = atomicCounter++;
        if (atomI >= numParticles)
            break;
        if (bornRadii[atomI] > 0) {
            float radiusI = particleParams[atomI].first + dielectricOffset;
            float r = radiusI + probeRadius;
            float ratio6 = powf(radiusI/bornRadii[atomI], 6.0f);
            float saTerm = surfaceAreaFactor*r*r*ratio6;
            energy += saTerm;
            bornForces[atomI] = -6.0f*saTerm/bornRadii[atomI]; 
        }
        else
            bornForces[atomI] = 0.0f;
    }
    threads.syncThreads();
 
    // First loop of Born energy computation.

    float* forces = &(*threadForce)[threadIndex][0];
    float preFactor;
    if (soluteDielectric != 0.0f && solventDielectric != 0.0f)
        preFactor = ONE_4PI_EPS0*((1.0f/solventDielectric) - (1.0f/soluteDielectric));
    else
        preFactor = 0.0f;
    while (true) {
        int blockStart = atomicCounter.fetch_add(4);
        if (blockStart >= numParticles)
            break;
        int numInBlock = min(4, numParticles-blockStart);
        ivec4 blockAtomIndex(blockStart, blockStart+1, blockStart+2, blockStart+3);
        float atomCharge[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        float atomx[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        float atomy[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        float atomz[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        int blockMask[4] = {0, 0, 0, 0};
        fvec4 blockAtomForceX(0.0f), blockAtomForceY(0.0f), blockAtomForceZ(0.0f), blockAtomBornForce(0.0f);
        for (int i = 0; i < numInBlock; i++) {
            int atomIndex = blockStart+i;
            atomx[i] = posq[4*atomIndex];
            atomy[i] = posq[4*atomIndex+1];
            atomz[i] = posq[4*atomIndex+2];
            atomCharge[i] = preFactor*posq[4*atomIndex+3];
            blockMask[i] = 0xFFFFFFFF;
        }
        fvec4 radii(&bornRadii[blockStart]);
        fvec4 x(atomx);
        fvec4 y(atomy);
        fvec4 z(atomz);
        fvec4 partialChargeI(atomCharge);
        ivec4 mask(blockMask);
        for (int atomJ = blockStart; atomJ < numParticles; atomJ++) {
            fvec4 posJ(posq+4*atomJ);
            fvec4 dx, dy, dz, r2;
            getDeltaR(posJ, x, y, z, dx, dy, dz, r2, periodic, boxSize, invBoxSize);
            ivec4 include = mask & (blockAtomIndex <= ivec4(atomJ));
            if (cutoff)
                include = include & (r2 < cutoffDistance*cutoffDistance);
            if (!any(include))
                continue;
            fvec4 r = sqrt(r2);
            fvec4 alpha2_ij = radii*bornRadii[atomJ];
            fvec4 D_ij = r2/(4.0f*alpha2_ij);
            fvec4 expTerm = exp(-D_ij);
            fvec4 denominator2 = r2 + alpha2_ij*expTerm;
            fvec4 denominator = sqrt(denominator2);
            fvec4 Gpol = (partialChargeI*posJ[3])/denominator; 
            fvec4 dGpol_dr = -Gpol*(1.0f - 0.25f*expTerm)/denominator2;  
            fvec4 dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f + D_ij)/denominator2;
            dGpol_dr = blend(0.0f, dGpol_dr, include);
            dGpol_dalpha2_ij = blend(0.0f, dGpol_dalpha2_ij, include);
            fvec4 fx = dx*dGpol_dr;
            fvec4 fy = dy*dGpol_dr;
            fvec4 fz = dz*dGpol_dr;
            blockAtomForceX -= fx;
            blockAtomForceY -= fy;
            blockAtomForceZ -= fz;
            blockAtomBornForce += dGpol_dalpha2_ij*bornRadii[atomJ];
            float* atomForce = forces+4*atomJ;
            fvec4 one(1.0f);
            atomForce[0] += dot4(fx, one);
            atomForce[1] += dot4(fy, one);
            atomForce[2] += dot4(fz, one);
            ivec4 atomJMask = include & (blockAtomIndex != ivec4(atomJ));
            fvec4 termEnergy = blend(0.0f, Gpol, include);
            if (cutoff)
                termEnergy -= blend(0.0f, partialChargeI*posJ[3]/cutoffDistance, atomJMask);
            termEnergy *= blend(0.5f, 1.0f, atomJMask);
            energy += dot4(termEnergy, one);
            bornForces[atomJ] += dot4(blend(0.0f, dGpol_dalpha2_ij, atomJMask), radii);
        }
        fvec4 f[4] = {blockAtomForceX, blockAtomForceY, blockAtomForceZ, 0.0f};
        transpose(f[0], f[1], f[2], f[3]);
        for (int i = 0; i < numInBlock; i++) {
            int atomIndex = blockStart+i;
            (fvec4(forces+4*atomIndex)+f[i]).store(forces+4*atomIndex);
            bornForces[atomIndex] += blockAtomBornForce[i];
        }
    }
    threads.syncThreads();

    // Second loop of Born energy computation.

    while (true) {
        int blockStart = atomicCounter.fetch_add(4);
        if (blockStart >= numParticles)
            break;
        fvec4 bornForce(0.0f);
        for (int i = 0; i < numThreads; i++)
            bornForce += fvec4(&threadBornForces[i][blockStart]);
        fvec4 radii(&bornRadii[blockStart]);
        bornForce *= radii*radii*fvec4(&obcChain[blockStart]);
        int numInBlock = min(4, numParticles-blockStart);
        ivec4 blockAtomIndex(blockStart, blockStart+1, blockStart+2, blockStart+3);
        float atomRadius[4], atomx[4], atomy[4], atomz[4];
        int blockMask[4] = {0, 0, 0, 0};
        fvec4 blockAtomForceX(0.0f), blockAtomForceY(0.0f), blockAtomForceZ(0.0f);
        for (int i = 0; i < numInBlock; i++) {
            int atomIndex = blockStart+i;
            atomRadius[i] = particleParams[atomIndex].first;
            atomx[i] = posq[4*atomIndex];
            atomy[i] = posq[4*atomIndex+1];
            atomz[i] = posq[4*atomIndex+2];
            blockMask[i] = 0xFFFFFFFF;
        }
        for (int i = numInBlock; i < 4; i++) {
            atomx[i] = 0.0f;
            atomy[i] = 0.0f;
            atomz[i] = 0.0f;
        }
        fvec4 offsetRadiusI(atomRadius);
        fvec4 x(atomx);
        fvec4 y(atomy);
        fvec4 z(atomz);
        ivec4 mask(blockMask);
        for (int atomJ = 0; atomJ < numParticles; atomJ++) {
            fvec4 posJ(posq+4*atomJ);
            fvec4 dx, dy, dz, r2;
            getDeltaR(posJ, x, y, z, dx, dy, dz, r2, periodic, boxSize, invBoxSize);
            ivec4 include = mask & (blockAtomIndex != ivec4(atomJ));
            if (cutoff)
                include = include & (r2 < cutoffDistance*cutoffDistance);
            if (!any(include))
                continue;
            fvec4 r = sqrt(r2);
            float scaledRadiusJ = particleParams[atomJ].second;
            float scaledRadiusJ2 = scaledRadiusJ*scaledRadiusJ;
            fvec4 rScaledRadiusJ = r + scaledRadiusJ;
            include = include & (offsetRadiusI < rScaledRadiusJ);
            fvec4 l_ij = 1.0f/max(offsetRadiusI, abs(r-scaledRadiusJ));
            fvec4 u_ij = 1.0f/rScaledRadiusJ;
            fvec4 l_ij2 = l_ij*l_ij;
            fvec4 u_ij2 = u_ij*u_ij;
            fvec4 rInverse = 1.0f/r;
            fvec4 r2Inverse = rInverse*rInverse;
            fvec4 logRatio = fastLog(u_ij/l_ij);
            fvec4 t3 = 0.125f*(1.0f + scaledRadiusJ2*r2Inverse)*(l_ij2 - u_ij2) + 0.25f*logRatio*r2Inverse;
            fvec4 de = bornForce*t3*rInverse;
            de = blend(0.0f, de, include);
            fvec4 fx = dx*de;
            fvec4 fy = dy*de;
            fvec4 fz = dz*de;
            blockAtomForceX += fx;
            blockAtomForceY += fy;
            blockAtomForceZ += fz;
            float* atomForce = forces+4*atomJ;
            fvec4 one(1.0f);
            atomForce[0] -= dot4(fx, one);
            atomForce[1] -= dot4(fy, one);
            atomForce[2] -= dot4(fz, one);
        }
        fvec4 f[4] = {blockAtomForceX, blockAtomForceY, blockAtomForceZ, 0.0f};
        transpose(f[0], f[1], f[2], f[3]);
        for (int i = 0; i < numInBlock; i++) {
            int atomIndex = blockStart+i;
            (fvec4(forces+4*atomIndex)+f[i]).store(forces+4*atomIndex);
        }
    }
    threadEnergy[threadIndex] = energy;
}

void CpuGBSAOBCForce::getDeltaR(const fvec4& posI, const fvec4& x, const fvec4& y, const fvec4& z, fvec4& dx, fvec4& dy, fvec4& dz, fvec4& r2, bool periodic, const fvec4& boxSize, const fvec4& invBoxSize) const {
    dx = x-posI[0];
    dy = y-posI[1];
    dz = z-posI[2];
    if (periodic) {
        dx -= round(dx*invBoxSize[0])*boxSize[0];
        dy -= round(dy*invBoxSize[1])*boxSize[1];
        dz -= round(dz*invBoxSize[2])*boxSize[2];
    }
    r2 = dx*dx + dy*dy + dz*dz;
}

fvec4 CpuGBSAOBCForce::fastLog(const fvec4& x) {
    // Evaluate log(x) using a lookup table for speed.

    fvec4 x1 = (x-TABLE_MIN)*logDXInv;
    ivec4 index = floor(x1);
    if (any((index < 0) | (index >= NUM_TABLE_POINTS)))
        return fvec4(logf(x[0]), logf(x[1]), logf(x[2]), logf(x[3]));
    fvec4 coeff2 = x1-index;
    fvec4 coeff1 = 1.0f-coeff2;
    fvec4 t1(&logTable[index[0]]);
    fvec4 t2(&logTable[index[1]]);
    fvec4 t3(&logTable[index[2]]);
    fvec4 t4(&logTable[index[3]]);
    transpose(t1, t2, t3, t4);
    return coeff1*t1 + coeff2*t2;
}
