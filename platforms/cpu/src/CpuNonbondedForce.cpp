
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

#include <complex>

#include "SimTKOpenMMUtilities.h"
#include "CpuNonbondedForce.h"
#include "ReferenceForce.h"
#include "ReferencePME.h"
#include <algorithm>
#include <iostream>

// In case we're using some primitive version of Visual Studio this will
// make sure that erf() and erfc() are defined.
#include "openmm/internal/MSVC_erfc.h"

using namespace std;
using namespace OpenMM;

const float CpuNonbondedForce::TWO_OVER_SQRT_PI = (float) (2/sqrt(PI_M));
const int CpuNonbondedForce::NUM_TABLE_POINTS = 2048;

/**---------------------------------------------------------------------------------------

   CpuNonbondedForce constructor

   --------------------------------------------------------------------------------------- */

CpuNonbondedForce::CpuNonbondedForce() : cutoff(false), useSwitch(false), periodic(false), ewald(false), pme(false), ljpme(false), tableIsValid(false), expTableIsValid(false),
    cutoffDistance(0.0f), alphaDispersionEwald(0.0f), alphaEwald(0.0f) {
}

CpuNonbondedForce::~CpuNonbondedForce() {
}

/**---------------------------------------------------------------------------------------

   Set the force to use a cutoff.

   @param distance            the cutoff distance
   @param neighbors           the neighbor list to use
   @param solventDielectric   the dielectric constant of the bulk solvent

     --------------------------------------------------------------------------------------- */

void CpuNonbondedForce::setUseCutoff(float distance, const CpuNeighborList& neighbors, float solventDielectric) {
    if (distance != cutoffDistance)
        tableIsValid = false;
    cutoff = true;
    cutoffDistance = distance;
    inverseRcut6 = pow(cutoffDistance, -6);
    neighborList = &neighbors;
    krf = pow(cutoffDistance, -3.0f)*(solventDielectric-1.0)/(2.0*solventDielectric+1.0);
    crf = (1.0/cutoffDistance)*(3.0*solventDielectric)/(2.0*solventDielectric+1.0);
    if(alphaDispersionEwald != 0.0f){
        // We set this here, in case setUseCutoff is called after the dispersion alpha is set.
        double dalphaR = alphaDispersionEwald*cutoffDistance;
        double dar2 = dalphaR * dalphaR;
        double dar4 = dar2*dar2;
        double dar6 = dar4*dar2;
        double expterm = EXP(-dar2);
        inverseRcut6Expterm  = inverseRcut6*(1.0 - expterm * (1.0 + dar2 + 0.5*dar4));
    }

}

/**---------------------------------------------------------------------------------------

   Set the force to use a switching function on the Lennard-Jones interaction.

   @param distance            the switching distance

   --------------------------------------------------------------------------------------- */

void CpuNonbondedForce::setUseSwitchingFunction(float distance) {
    useSwitch = true;
    switchingDistance = distance;
}

/**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param periodicBoxVectors    the vectors defining the periodic box

     --------------------------------------------------------------------------------------- */

void CpuNonbondedForce::setPeriodic(Vec3* periodicBoxVectors) {

    assert(cutoff);
    assert(periodicBoxVectors[0][0] >= 2.0*cutoffDistance);
    assert(periodicBoxVectors[1][1] >= 2.0*cutoffDistance);
    assert(periodicBoxVectors[2][2] >= 2.0*cutoffDistance);
    periodic = true;
    this->periodicBoxVectors[0] = periodicBoxVectors[0];
    this->periodicBoxVectors[1] = periodicBoxVectors[1];
    this->periodicBoxVectors[2] = periodicBoxVectors[2];
    recipBoxSize[0] = (float) (1.0/periodicBoxVectors[0][0]);
    recipBoxSize[1] = (float) (1.0/periodicBoxVectors[1][1]);
    recipBoxSize[2] = (float) (1.0/periodicBoxVectors[2][2]);
    periodicBoxVec4.resize(3);
    periodicBoxVec4[0] = fvec4(periodicBoxVectors[0][0], periodicBoxVectors[0][1], periodicBoxVectors[0][2], 0);
    periodicBoxVec4[1] = fvec4(periodicBoxVectors[1][0], periodicBoxVectors[1][1], periodicBoxVectors[1][2], 0);
    periodicBoxVec4[2] = fvec4(periodicBoxVectors[2][0], periodicBoxVectors[2][1], periodicBoxVectors[2][2], 0);
    triclinic = (periodicBoxVectors[0][1] != 0.0 || periodicBoxVectors[0][2] != 0.0 ||
            periodicBoxVectors[1][0] != 0.0 || periodicBoxVectors[1][2] != 0.0 ||
            periodicBoxVectors[2][0] != 0.0 || periodicBoxVectors[2][1] != 0.0);
}

/**---------------------------------------------------------------------------------------

     Set the force to use Ewald summation.

     @param alpha  the Ewald separation parameter
     @param kmaxx  the largest wave vector in the x direction
     @param kmaxy  the largest wave vector in the y direction
     @param kmaxz  the largest wave vector in the z direction

     --------------------------------------------------------------------------------------- */

void CpuNonbondedForce::setUseEwald(float alpha, int kmaxx, int kmaxy, int kmaxz) {
    if (alpha != alphaEwald)
        tableIsValid = false;
    alphaEwald = alpha;
    numRx = kmaxx;
    numRy = kmaxy;
    numRz = kmaxz;
    ewald = true;
    tabulateEwaldScaleFactor();
}

/**---------------------------------------------------------------------------------------

     Set the force to use Particle-Mesh Ewald (PME) summation.

     @param alpha  the Ewald separation parameter
     @param gridSize the dimensions of the mesh

     --------------------------------------------------------------------------------------- */

void CpuNonbondedForce::setUsePME(float alpha, int meshSize[3]) {
    if (alpha != alphaEwald)
        tableIsValid = false;
    alphaEwald = alpha;
    meshDim[0] = meshSize[0];
    meshDim[1] = meshSize[1];
    meshDim[2] = meshSize[2];
    pme = true;
    tabulateEwaldScaleFactor();
}


/**---------------------------------------------------------------------------------------

     Set the force to use Particle-Mesh Ewald (PME) summation for dispersion.

     @param alpha  the Ewald separation parameter
     @param gridSize the dimensions of the mesh

     --------------------------------------------------------------------------------------- */

void CpuNonbondedForce::setUseLJPME(float alpha, int meshSize[3]) {
    if (alpha != alphaDispersionEwald)
        expTableIsValid = false;
    alphaDispersionEwald = alpha;
    dispersionMeshDim[0] = meshSize[0];
    dispersionMeshDim[1] = meshSize[1];
    dispersionMeshDim[2] = meshSize[2];
    ljpme = true;
    tabulateExpTerms();
    if(cutoffDistance != 0.0f){
        // We set this here, in case setUseLJPME is called after the cutoff is set
        double dalphaR = alphaDispersionEwald*cutoffDistance;
        double dar2 = dalphaR * dalphaR;
        double dar4 = dar2*dar2;
        double dar6 = dar4*dar2;
        double expterm = EXP(-dar2);
        inverseRcut6Expterm  = inverseRcut6*(1.0 - expterm * (1.0 + dar2 + 0.5*dar4));
    }
}


void CpuNonbondedForce::tabulateEwaldScaleFactor() {
    if (tableIsValid)
        return;
    tableIsValid = true;
    ewaldDX = cutoffDistance/NUM_TABLE_POINTS;
    ewaldDXInv = 1.0f/ewaldDX;
    erfcDXInv = 1.0f/(ewaldDX*alphaEwald);
    erfcTable.resize(NUM_TABLE_POINTS+4);
    ewaldScaleTable.resize(NUM_TABLE_POINTS+4);
    for (int i = 0; i < NUM_TABLE_POINTS+4; i++) {
        double r = i*ewaldDX;
        double alphaR = alphaEwald*r;
        erfcTable[i] = erfc(alphaR);
        ewaldScaleTable[i] = erfcTable[i] + TWO_OVER_SQRT_PI*alphaR*exp(-alphaR*alphaR);
    }
}

void CpuNonbondedForce::tabulateExpTerms() {
    if (expTableIsValid)
        return;
    expTableIsValid = true;
    exptermsDX = cutoffDistance/NUM_TABLE_POINTS;
    exptermsDXInv = 1.0f/exptermsDX;
    exptermsTable.resize(NUM_TABLE_POINTS+4);
    dExptermsTable.resize(NUM_TABLE_POINTS+4);
    for (int i = 0; i < NUM_TABLE_POINTS+4; i++) {
        double r = i*ewaldDX;
        double dalphaR = alphaDispersionEwald*r;
        double dar2 = dalphaR * dalphaR;
        double dar4 = dar2*dar2;
        double dar6 = dar4*dar2;
        double expterm = EXP(-dar2);
        exptermsTable[i]  = (1.0 - expterm * (1.0 + dar2 + 0.5*dar4));
        dExptermsTable[i] = (1.0 - expterm * (1.0 + dar2 + 0.5*dar4 + dar6/6.0));
    }
}

void CpuNonbondedForce::calculateReciprocalIxn(int numberOfAtoms, float* posq, const vector<Vec3>& atomCoordinates,
                                               const vector<pair<float, float> >& atomParameters, const vector<float> &C6params, const vector<set<int> >& exclusions,
                                               vector<Vec3>& forces, double* totalEnergy) const {
    typedef std::complex<float> d_complex;

    static const float epsilon     =  1.0;

    int kmax                       = (ewald ? max(numRx, max(numRy,numRz)) : 0);
    float factorEwald              = -1 / (4*alphaEwald*alphaEwald);
    float TWO_PI                   = 2.0 * PI_M;
    float recipCoeff               = (float)(ONE_4PI_EPS0*4*PI_M/(periodicBoxVectors[0][0] * periodicBoxVectors[1][1] * periodicBoxVectors[2][2]) /epsilon);

    if (pme) {
        pme_t pmedata;
        pme_init(&pmedata, alphaEwald, numberOfAtoms, meshDim, 5, 1);
        vector<double> charges(numberOfAtoms);
        for (int i = 0; i < numberOfAtoms; i++)
            charges[i] = posq[4*i+3];
        double recipEnergy = 0.0;
        pme_exec(pmedata, atomCoordinates, forces, charges, periodicBoxVectors, &recipEnergy);
        if (totalEnergy)
            *totalEnergy += recipEnergy;
        pme_destroy(pmedata);

        if (ljpme) {
            // Dispersion reciprocal space terms
            pme_init(&pmedata,alphaDispersionEwald,numberOfAtoms,dispersionMeshDim,5,1);

            std::vector<Vec3> dpmeforces;
            for (int i = 0; i < numberOfAtoms; i++){
                charges[i] = C6params[i];
                dpmeforces.push_back(Vec3());
            }
            double recipDispersionEnergy = 0.0;
            pme_exec_dpme(pmedata,atomCoordinates,dpmeforces,charges,periodicBoxVectors,&recipDispersionEnergy);
            for (int i = 0; i < numberOfAtoms; i++){
                forces[i][0] += dpmeforces[i][0];
                forces[i][1] += dpmeforces[i][1];
                forces[i][2] += dpmeforces[i][2];
            }
            if (totalEnergy)
                *totalEnergy += recipDispersionEnergy;

            pme_destroy(pmedata);
        }

    }

    // Ewald method

    else if (ewald) {

        // setup reciprocal box

        float recipBoxSize[3] = {(float) (TWO_PI/periodicBoxVectors[0][0]), (float) (TWO_PI/periodicBoxVectors[1][1]), (float) (TWO_PI/periodicBoxVectors[2][2])};


        // setup K-vectors

#define EIR(x, y, z) eir[(x)*numberOfAtoms*3+(y)*3+z]
        vector<d_complex> eir(kmax*numberOfAtoms*3);
        vector<d_complex> tab_xy(numberOfAtoms);
        vector<d_complex> tab_qxyz(numberOfAtoms);

        for (int i = 0; (i < numberOfAtoms); i++) {
            float* pos = posq+4*i;
            for (int m = 0; (m < 3); m++)
                EIR(0, i, m) = d_complex(1,0);

            for (int m=0; (m<3); m++)
                EIR(1, i, m) = d_complex(cos(pos[m]*recipBoxSize[m]),
                                         sin(pos[m]*recipBoxSize[m]));

            for (int j=2; (j<kmax); j++)
                for (int m=0; (m<3); m++)
                    EIR(j, i, m) = EIR(j-1, i, m) * EIR(1, i, m);
        }

        // calculate reciprocal space energy and forces

        int lowry = 0;
        int lowrz = 1;

        for (int rx = 0; rx < numRx; rx++) {
            float kx = rx * recipBoxSize[0];
            for (int ry = lowry; ry < numRy; ry++) {
                float ky = ry * recipBoxSize[1];
                if (ry >= 0) {
                    for (int n = 0; n < numberOfAtoms; n++)
                        tab_xy[n] = EIR(rx, n, 0) * EIR(ry, n, 1);
                }
                else {
                    for (int n = 0; n < numberOfAtoms; n++)
                        tab_xy[n]= EIR(rx, n, 0) * conj (EIR(-ry, n, 1));
                }
                for (int rz = lowrz; rz < numRz; rz++) {
                    if (rz >= 0) {
                        for (int n = 0; n < numberOfAtoms; n++)
                            tab_qxyz[n] = posq[4*n+3] * (tab_xy[n] * EIR(rz, n, 2));
                    }
                    else {
                        for (int n = 0; n < numberOfAtoms; n++)
                            tab_qxyz[n] = posq[4*n+3] * (tab_xy[n] * conj(EIR(-rz, n, 2)));
                    }
                    float cs = 0.0f;
                    float ss = 0.0f;

                    for (int n = 0; n < numberOfAtoms; n++) {
                        cs += tab_qxyz[n].real();
                        ss += tab_qxyz[n].imag();
                    }

                    float kz = rz * recipBoxSize[2];
                    float k2 = kx * kx + ky * ky + kz * kz;
                    float ak = exp(k2*factorEwald) / k2;

                    for (int n = 0; n < numberOfAtoms; n++) {
                        float force = ak * (cs * tab_qxyz[n].imag() - ss * tab_qxyz[n].real());
                        forces[n][0] += 2 * recipCoeff * force * kx;
                        forces[n][1] += 2 * recipCoeff * force * ky;
                        forces[n][2] += 2 * recipCoeff * force * kz;
                    }

                    if (totalEnergy)
                        *totalEnergy += recipCoeff * ak * (cs * cs + ss * ss);

                    lowrz = 1 - numRz;
                }
                lowry = 1 - numRy;
            }
        }
    }
}


void CpuNonbondedForce::calculateDirectIxn(int numberOfAtoms, float* posq, const vector<Vec3>& atomCoordinates, const vector<pair<float, float> >& atomParameters,
                                           const vector<float>& C6params, const vector<set<int> >& exclusions, vector<AlignedArray<float> >& threadForce, double* totalEnergy, ThreadPool& threads) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->posq = posq;
    this->atomCoordinates = &atomCoordinates[0];
    this->atomParameters = &atomParameters[0];
    this->C6params = &C6params[0];
    this->exclusions = &exclusions[0];
    this->threadForce = &threadForce;
    includeEnergy = (totalEnergy != NULL);
    threadEnergy.resize(threads.getNumThreads());
    atomicCounter = 0;
    
    // Signal the threads to start running and wait for them to finish.
    
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadComputeDirect(threads, threadIndex); });
    threads.waitForThreads();
    
    // Signal the threads to subtract the exclusions.
    
    if (ewald || pme) {
        atomicCounter = 0;
        threads.resumeThreads();
        threads.waitForThreads();
    }
    
    // Combine the energies from all the threads.
    
    if (totalEnergy != NULL) {
        double directEnergy = 0;
        int numThreads = threads.getNumThreads();
        for (int i = 0; i < numThreads; i++)
            directEnergy += threadEnergy[i];
        *totalEnergy += directEnergy;
    }
}

void CpuNonbondedForce::threadComputeDirect(ThreadPool& threads, int threadIndex) {
    // Compute this thread's subset of interactions.

    int numThreads = threads.getNumThreads();
    threadEnergy[threadIndex] = 0;
    double* energyPtr = (includeEnergy ? &threadEnergy[threadIndex] : NULL);
    float* forces = &(*threadForce)[threadIndex][0];
    fvec4 boxSize(periodicBoxVectors[0][0], periodicBoxVectors[1][1], periodicBoxVectors[2][2], 0);
    fvec4 invBoxSize(recipBoxSize[0], recipBoxSize[1], recipBoxSize[2], 0);
    if (ewald || pme || ljpme) {
        // Compute the interactions from the neighbor list.
        while (true) {
            int nextBlock = atomicCounter++;
            if (nextBlock >= neighborList->getNumBlocks())
                break;
            calculateBlockEwaldIxn(nextBlock, forces, energyPtr, boxSize, invBoxSize);
        }

        // Now subtract off the exclusions, since they were implicitly included in the reciprocal space sum.

        threads.syncThreads();
        const int groupSize = max(1, numberOfAtoms/(10*numThreads));
        while (true) {
            int start = atomicCounter.fetch_add(groupSize);
            if (start >= numberOfAtoms)
                break;
            int end = min(start+groupSize, numberOfAtoms);
            for (int i = start; i < end; i++) {
                fvec4 posI((float) atomCoordinates[i][0], (float) atomCoordinates[i][1], (float) atomCoordinates[i][2], 0.0f);
                float scaledChargeI = (float) (ONE_4PI_EPS0*posq[4*i+3]);
                for (int excluded : exclusions[i]) {
                    if (excluded > i) {
                        int j = excluded;
                        fvec4 deltaR;
                        fvec4 posJ((float) atomCoordinates[j][0], (float) atomCoordinates[j][1], (float) atomCoordinates[j][2], 0.0f);
                        float r2;
                        getDeltaR(posJ, posI, deltaR, r2, false, boxSize, invBoxSize);
                        float r = sqrtf(r2);
                        float alphaR = alphaEwald*r;
                        float erfAlphaR = erf(alphaR);
                        if (erfAlphaR > 1e-6f) {
                            float inverseR = 1/r;
                            float chargeProdOverR = scaledChargeI*posq[4*j+3]*inverseR;
                            float dEdR = chargeProdOverR*inverseR*inverseR;
                            dEdR = dEdR * (erfAlphaR-TWO_OVER_SQRT_PI*alphaR*(float)exp(-alphaR*alphaR));
                            fvec4 result = deltaR*dEdR;
                            (fvec4(forces+4*i)-result).store(forces+4*i);
                            (fvec4(forces+4*j)+result).store(forces+4*j);
                            if (includeEnergy)
                                threadEnergy[threadIndex] -= chargeProdOverR*erfAlphaR;
                        }
                        else if (includeEnergy)
                            threadEnergy[threadIndex] -= alphaEwald*TWO_OVER_SQRT_PI*scaledChargeI*posq[4*j+3];
                        if (ljpme) {
                            float C6ij = C6params[i]*C6params[j];
                            float inverseR2 = 1.0f/r2;
                            float emult = C6ij*inverseR2*inverseR2*inverseR2*exptermsApprox(r);
                            if(includeEnergy)
                                threadEnergy[threadIndex] += emult;
                            float dEdR = -6.0f*C6ij*inverseR2*inverseR2*inverseR2*inverseR2*dExptermsApprox(r);
                            fvec4 result = deltaR*dEdR;
                            (fvec4(forces+4*i)-result).store(forces+4*i);
                            (fvec4(forces+4*j)+result).store(forces+4*j);
                        }
                    }
                }
            }
        }
    }
    else if (cutoff) {
        // Compute the interactions from the neighbor list.

        while (true) {
            int nextBlock = atomicCounter++;
            if (nextBlock >= neighborList->getNumBlocks())
                break;
            calculateBlockIxn(nextBlock, forces, energyPtr, boxSize, invBoxSize);
        }
    }
    else {
        // Loop over all atom pairs

        while (true) {
            int i = atomicCounter++;
            if (i >= numberOfAtoms)
                break;
            for (int j = i+1; j < numberOfAtoms; j++)
                if (exclusions[j].find(i) == exclusions[j].end())
                    calculateOneIxn(i, j, forces, energyPtr, boxSize, invBoxSize);
        }
    }
}

void CpuNonbondedForce::calculateOneIxn(int ii, int jj, float* forces, double* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    // get deltaR, R2, and R between 2 atoms

    fvec4 deltaR;
    fvec4 posI(posq+4*ii);
    fvec4 posJ(posq+4*jj);
    float r2;
    getDeltaR(posJ, posI, deltaR, r2, periodic, boxSize, invBoxSize);
    if (cutoff && r2 >= cutoffDistance*cutoffDistance)
        return;
    float r = sqrtf(r2);
    float inverseR = 1/r;
    float switchValue = 1, switchDeriv = 0;
    if (useSwitch && r > switchingDistance) {
        float t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
        switchValue = 1+t*t*t*(-10+t*(15-t*6));
        switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
    }
    float sig       = atomParameters[ii].first + atomParameters[jj].first;
    float sig2      = inverseR*sig;
    sig2     *= sig2;
    float sig6      = sig2*sig2*sig2;

    float eps       = atomParameters[ii].second*atomParameters[jj].second;
    float dEdR      = switchValue*eps*(12.0f*sig6 - 6.0f)*sig6;
    float chargeProd = ONE_4PI_EPS0*posq[4*ii+3]*posq[4*jj+3];
    if (cutoff)
        dEdR += (float) (chargeProd*(inverseR-2.0f*krf*r2));
    else
        dEdR += (float) (chargeProd*inverseR);
    dEdR *= inverseR*inverseR;
    float energy = eps*(sig6-1.0f)*sig6;
    if (useSwitch) {
        dEdR -= energy*switchDeriv*inverseR;
        energy *= switchValue;
    }

    // accumulate energies

    if (totalEnergy) {
        if (cutoff)
            energy += (float) (chargeProd*(inverseR+krf*r2-crf));
        else
            energy += (float) (chargeProd*inverseR);
        *totalEnergy += energy;
    }

    // accumulate forces

    fvec4 result = deltaR*dEdR;
    (fvec4(forces+4*ii)+result).store(forces+4*ii);
    (fvec4(forces+4*jj)-result).store(forces+4*jj);
}

void CpuNonbondedForce::getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, bool periodic, const fvec4& boxSize, const fvec4& invBoxSize) const {
    deltaR = posJ-posI;
    if (periodic) {
        if (triclinic) {
            deltaR -= periodicBoxVec4[2]*floorf(deltaR[2]*recipBoxSize[2]+0.5f);
            deltaR -= periodicBoxVec4[1]*floorf(deltaR[1]*recipBoxSize[1]+0.5f);
            deltaR -= periodicBoxVec4[0]*floorf(deltaR[0]*recipBoxSize[0]+0.5f);
        }
        else {
            fvec4 base = round(deltaR*invBoxSize)*boxSize;
            deltaR = deltaR-base;
        }
    }
    r2 = dot3(deltaR, deltaR);
}

float CpuNonbondedForce::erfcApprox(float x) {
    float x1 = x*erfcDXInv;
    int index = min((int) floor(x1), NUM_TABLE_POINTS);
    float coeff2 = x1-index;
    float coeff1 = 1.0f-coeff2;
    return coeff1*erfcTable[index] + coeff2*erfcTable[index+1];
}

float CpuNonbondedForce::exptermsApprox(float x) {
    float x1 = x*exptermsDXInv;
    int index = min((int) floor(x1), NUM_TABLE_POINTS);
    float coeff2 = x1-index;
    float coeff1 = 1.0f-coeff2;
    return coeff1*exptermsTable[index] + coeff2*exptermsTable[index+1];
}

float CpuNonbondedForce::dExptermsApprox(float x) {
    float x1 = x*exptermsDXInv;
    int index = min((int) floor(x1), NUM_TABLE_POINTS);
    float coeff2 = x1-index;
    float coeff1 = 1.0f-coeff2;
    return coeff1*dExptermsTable[index] + coeff2*dExptermsTable[index+1];
}
