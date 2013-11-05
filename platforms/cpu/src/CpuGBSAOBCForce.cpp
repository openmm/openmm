
/* Portions copyright (c) 2006-2013 Stanford University and Simbios.
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
#include <cmath>

using namespace std;
using namespace OpenMM;

class CpuGBSAOBCForce::ComputeTask : public ThreadPool::Task {
public:
    ComputeTask(CpuGBSAOBCForce& owner) : owner(owner) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        owner.threadComputeForce(threads, threadIndex);
    }
    CpuGBSAOBCForce& owner;
};

CpuGBSAOBCForce::CpuGBSAOBCForce() : cutoff(false), periodic(false) {
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

const std::vector<std::pair<float, float> >& CpuGBSAOBCForce::getParticleParameters() const {
    return particleParams;
}

void CpuGBSAOBCForce::setParticleParameters(const std::vector<std::pair<float, float> >& params) {
    particleParams = params;
    bornRadii.resize(params.size());
    obcChain.resize(params.size());
}

void CpuGBSAOBCForce::computeForce(const std::vector<float>& posq, vector<vector<float> >& threadForce, double* totalEnergy, ThreadPool& threads) {
    // Record the parameters for the threads.
    
    this->posq = &posq[0];
    this->threadForce = &threadForce;
    includeEnergy = (totalEnergy != NULL);
    int numThreads = threads.getNumThreads();
    threadEnergy.resize(numThreads);
    threadBornForces.resize(numThreads);
    for (int i = 0; i < numThreads; i++)
        threadBornForces[i].resize(particleParams.size());
    
    // Signal the threads to start running and wait for them to finish.
    
    ComputeTask task(*this);
    threads.execute(task);
    threads.waitForThreads(); // Compute Born radii
    threads.resumeThreads();
    threads.waitForThreads(); // Compute surface area term
    threads.resumeThreads();
    threads.waitForThreads(); // First loop
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

    for (int atomI = start; atomI < end; atomI++) {
        fvec4 posI(posq+4*atomI);
        float offsetRadiusI = particleParams[atomI].first;
        float radiusIInverse = 1.0f/offsetRadiusI;
        float sum = 0.0f;
        for (int atomJ = 0; atomJ < numParticles; atomJ++) {
            if (atomJ != atomI) {
                fvec4 posJ(posq+4*atomJ);
                fvec4 deltaR;
                float r2;
                getDeltaR(posI, posJ, deltaR, r2, periodic, boxSize, invBoxSize);
                if (cutoff && r2 >= cutoffDistance*cutoffDistance)
                    continue;
                float r = sqrtf(r2);
                float scaledRadiusJ = particleParams[atomJ].second;
                float rScaledRadiusJ = r + scaledRadiusJ;
                if (offsetRadiusI < rScaledRadiusJ) {
                    float rInverse = 1.0f/r;
                    float l_ij = 1.0f/(offsetRadiusI > fabs(r - scaledRadiusJ) ? offsetRadiusI : fabs(r - scaledRadiusJ));
                    float u_ij = 1.0f/rScaledRadiusJ;
                    float l_ij2 = l_ij*l_ij;
                    float u_ij2 = u_ij*u_ij;
                    float ratio = log((u_ij/l_ij));
                    float term = l_ij - u_ij + 0.25f*r*(u_ij2 - l_ij2) + (0.5f*rInverse*ratio) + (0.25f*scaledRadiusJ*scaledRadiusJ*rInverse)*(l_ij2 - u_ij2);
                    if (offsetRadiusI < (scaledRadiusJ - r))
                        term += 2.0f*(radiusIInverse - l_ij);
                    sum += term;

                }
            }
        }
        sum *= 0.5f*offsetRadiusI;
        float sum2 = sum*sum;
        float sum3 = sum*sum2;
        float tanhSum = tanh(alphaObc*sum - betaObc*sum2 + gammaObc*sum3);
        float radiusI = offsetRadiusI + dielectricOffset;
        bornRadii[atomI] = 1.0f/(1.0f/offsetRadiusI - tanhSum/radiusI);
        obcChain[atomI] = offsetRadiusI*(alphaObc - 2.0f*betaObc*sum + 3.0f*gammaObc*sum2);
        obcChain[atomI] = (1.0f - tanhSum*tanhSum)*obcChain[atomI]/radiusI;
    }
    threads.syncThreads();

    // Calculate ACE surface area term.

    const float probeRadius = 0.14f;
    const float surfaceAreaFactor = 28.3919551;
    double energy = 0.0;
    vector<float>& bornForces = threadBornForces[threadIndex];
    for (int i = 0; i < numParticles; i++)
        bornForces[i] = 0.0f;
    for (int atomI = start; atomI < end; atomI++) {
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
    for (int atomI = start; atomI < end; atomI++) {
        fvec4 posI(posq+4*atomI);
        fvec4 forceI(0.0f);
        float partialChargeI = preFactor*posI[3];
        for (int atomJ = atomI; atomJ < numParticles; atomJ++) {
            fvec4 posJ(posq+4*atomJ);
            fvec4 deltaR;
            float r2;
            getDeltaR(posI, posJ, deltaR, r2, periodic, boxSize, invBoxSize);
            if (cutoff && r2 >= cutoffDistance*cutoffDistance)
                continue;
            float r = sqrtf(r2);
            float alpha2_ij = bornRadii[atomI]*bornRadii[atomJ];
            float D_ij = r2/(4.0f*alpha2_ij);
            float expTerm = exp(-D_ij);
            float denominator2 = r2 + alpha2_ij*expTerm; 
            float denominator = sqrt(denominator2);
            float Gpol = (partialChargeI*posJ[3])/denominator; 
            float dGpol_dr = -Gpol*(1.0f - 0.25f*expTerm)/denominator2;  
            float dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f + D_ij)/denominator2;
            float termEnergy = Gpol;
            if (atomI != atomJ) {
                if (cutoff)
                    termEnergy -= partialChargeI*posJ[3]/cutoffDistance;
                bornForces[atomJ] += dGpol_dalpha2_ij*bornRadii[atomI];
                fvec4 result = deltaR*dGpol_dr;
                forceI += result;
                (fvec4(forces+4*atomJ)-result).store(forces+4*atomJ);
            }
            else
                termEnergy *= 0.5f;
          energy += termEnergy;
          bornForces[atomI] += dGpol_dalpha2_ij*bornRadii[atomJ];

        }
        (fvec4(forces+4*atomI)+forceI).store(forces+4*atomI);
    }
    threads.syncThreads();

    // Second loop of Born energy computation.

    for (int atomI = start; atomI < end; atomI++) {
        float bornForce = 0;
        for (int i = 0; i < numThreads; i++)
            bornForce += threadBornForces[i][atomI];
        bornForce *= bornRadii[atomI]*bornRadii[atomI]*obcChain[atomI];      
        float offsetRadiusI = particleParams[atomI].first;
        fvec4 posI(posq+4*atomI);
        fvec4 forceI(0.0f);

        for (int atomJ = 0; atomJ < numParticles; atomJ++) {
            if (atomJ != atomI) {
                fvec4 posJ(posq+4*atomJ);
                fvec4 deltaR;
                float r2;
                getDeltaR(posI, posJ, deltaR, r2, periodic, boxSize, invBoxSize);
                if (cutoff && r2 >= cutoffDistance*cutoffDistance)
                    continue;
                float r = sqrtf(r2);
                float scaledRadiusJ = particleParams[atomJ].second;
                float scaledRadiusJ2 = scaledRadiusJ*scaledRadiusJ;
                float rScaledRadiusJ = r + scaledRadiusJ;
                if (offsetRadiusI < rScaledRadiusJ) {
                    float l_ij = 1.0f/(offsetRadiusI > fabs(r - scaledRadiusJ) ? offsetRadiusI : fabs(r - scaledRadiusJ));
                    float u_ij = 1.0f/rScaledRadiusJ;
                    float l_ij2 = l_ij*l_ij;
                    float u_ij2 = u_ij*u_ij;
                    float rInverse = 1.0f/r;
                    float r2Inverse = rInverse*rInverse;
                    float t3 = 0.125f*(1.0f + scaledRadiusJ2*r2Inverse)*(l_ij2 - u_ij2) + 0.25f*log(u_ij/l_ij)*r2Inverse;
                    float de = bornForce*t3*rInverse;
                    fvec4 result = deltaR*de;
                    forceI -= result;
                    (fvec4(forces+4*atomJ)+result).store(forces+4*atomJ);
                }
            }
        }
        (fvec4(forces+4*atomI)+forceI).store(forces+4*atomI);
    }
    threadEnergy[threadIndex] = energy;
}

void CpuGBSAOBCForce::getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, bool periodic, const fvec4& boxSize, const fvec4& invBoxSize) const {
    deltaR = posJ-posI;
    if (periodic) {
        fvec4 base = round(deltaR*invBoxSize)*boxSize;
        deltaR = deltaR-base;
    }
    r2 = dot3(deltaR, deltaR);
}
