
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

#include <string.h>
#include <complex>

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMUtilities.h"
#include "CpuNonbondedForce.h"
#include "ReferenceForce.h"
#include "ReferencePME.h"
#include "openmm/internal/hardware.h"

// In case we're using some primitive version of Visual Studio this will
// make sure that erf() and erfc() are defined.
#include "openmm/internal/MSVC_erfc.h"

using namespace std;

float CpuNonbondedForce::TWO_OVER_SQRT_PI = (float) (2/sqrt(PI_M));


class CpuNonbondedForce::ThreadData {
public:
    ThreadData(int index, CpuNonbondedForce& owner) : index(index), owner(owner) {
    }
    int index;
    CpuNonbondedForce& owner;
    vector<float> threadForce;
    double threadEnergy;
};

static void* threadBody(void* args) {
    CpuNonbondedForce::ThreadData& data = *reinterpret_cast<CpuNonbondedForce::ThreadData*>(args);
    data.owner.runThread(data.index, data.threadForce, data.threadEnergy);
    delete &data;
    return 0;
}

/**---------------------------------------------------------------------------------------

   CpuNonbondedForce constructor

   --------------------------------------------------------------------------------------- */

CpuNonbondedForce::CpuNonbondedForce() : cutoff(false), useSwitch(false), periodic(false), ewald(false), pme(false) {
    isDeleted = false;
    numThreads = getNumProcessors();
    pthread_cond_init(&startCondition, NULL);
    pthread_cond_init(&endCondition, NULL);
    pthread_mutex_init(&lock, NULL);
    thread.resize(numThreads);
    pthread_mutex_lock(&lock);
    waitCount = 0;
    for (int i = 0; i < numThreads; i++) {
        ThreadData* data = new ThreadData(i, *this);
        threadData.push_back(data);
        pthread_create(&thread[i], NULL, threadBody, data);
    }
    while (waitCount < numThreads)
        pthread_cond_wait(&endCondition, &lock);
    pthread_mutex_unlock(&lock);
}

/**---------------------------------------------------------------------------------------

   CpuNonbondedForce destructor

   --------------------------------------------------------------------------------------- */

CpuNonbondedForce::~CpuNonbondedForce(){
    isDeleted = true;
    pthread_mutex_lock(&lock);
    pthread_cond_broadcast(&startCondition);
    pthread_mutex_unlock(&lock);
    for (int i = 0; i < (int) thread.size(); i++)
        pthread_join(thread[i], NULL);
    pthread_mutex_destroy(&lock);
    pthread_cond_destroy(&startCondition);
    pthread_cond_destroy(&endCondition);
}

  /**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance
     @param neighbors           the neighbor list to use
     @param solventDielectric   the dielectric constant of the bulk solvent

     --------------------------------------------------------------------------------------- */

  void CpuNonbondedForce::setUseCutoff(float distance, const vector<pair<int, int> >& neighbors, float solventDielectric) {

    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;
    krf = pow(cutoffDistance, -3.0f)*(solventDielectric-1.0)/(2.0*solventDielectric+1.0);
    crf = (1.0/cutoffDistance)*(3.0*solventDielectric)/(2.0*solventDielectric+1.0);
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

     @param boxSize             the X, Y, and Z widths of the periodic box

     --------------------------------------------------------------------------------------- */

  void CpuNonbondedForce::setPeriodic(float* periodicBoxSize) {

    assert(cutoff);
    assert(periodicBoxSize[0] >= 2*cutoffDistance);
    assert(periodicBoxSize[1] >= 2*cutoffDistance);
    assert(periodicBoxSize[2] >= 2*cutoffDistance);
    periodic = true;
    this->periodicBoxSize[0] = periodicBoxSize[0];
    this->periodicBoxSize[1] = periodicBoxSize[1];
    this->periodicBoxSize[2] = periodicBoxSize[2];
    boxSize = _mm_set_ps(0, periodicBoxSize[2], periodicBoxSize[1], periodicBoxSize[0]);
    invBoxSize = _mm_set_ps(0, (1/periodicBoxSize[2]), (1/periodicBoxSize[1]), (1/periodicBoxSize[0]));
    half = _mm_set1_ps(0.5);
  }

  /**---------------------------------------------------------------------------------------

     Set the force to use Ewald summation.

     @param alpha  the Ewald separation parameter
     @param kmaxx  the largest wave vector in the x direction
     @param kmaxy  the largest wave vector in the y direction
     @param kmaxz  the largest wave vector in the z direction

     --------------------------------------------------------------------------------------- */

  void CpuNonbondedForce::setUseEwald(float alpha, int kmaxx, int kmaxy, int kmaxz) {
      alphaEwald = alpha;
      numRx = kmaxx;
      numRy = kmaxy;
      numRz = kmaxz;
      ewald = true;
  }

  /**---------------------------------------------------------------------------------------

     Set the force to use Particle-Mesh Ewald (PME) summation.

     @param alpha  the Ewald separation parameter
     @param gridSize the dimensions of the mesh

     --------------------------------------------------------------------------------------- */

  void CpuNonbondedForce::setUsePME(float alpha, int meshSize[3]) {
      alphaEwald = alpha;
      meshDim[0] = meshSize[0];
      meshDim[1] = meshSize[1];
      meshDim[2] = meshSize[2];
      pme = true;
  }

void CpuNonbondedForce::calculateReciprocalIxn(int numberOfAtoms, float* posq, vector<OpenMM::RealVec>& atomCoordinates,
                                             const vector<pair<float, float> >& atomParameters, const vector<set<int> >& exclusions,
                                             vector<OpenMM::RealVec>& forces, float* totalEnergy) const {
    typedef std::complex<float> d_complex;

    static const float epsilon     =  1.0;

    int kmax                            = (ewald ? std::max(numRx, std::max(numRy,numRz)) : 0);
    float factorEwald              = -1 / (4*alphaEwald*alphaEwald);
    float TWO_PI                   = 2.0 * PI_M;
    float recipCoeff               = (float)(ONE_4PI_EPS0*4*PI_M/(periodicBoxSize[0] * periodicBoxSize[1] * periodicBoxSize[2]) /epsilon);

    if (pme) {
        pme_t pmedata;
        RealOpenMM virial[3][3];
        pme_init(&pmedata, alphaEwald, numberOfAtoms, meshDim, 5, 1);
        vector<RealOpenMM> charges(numberOfAtoms);
        for (int i = 0; i < numberOfAtoms; i++)
            charges[i] = posq[4*i+3];
        RealOpenMM boxSize[3] = {periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2]};
        RealOpenMM recipEnergy = 0.0;
        pme_exec(pmedata, atomCoordinates, forces, charges, boxSize, &recipEnergy, virial);
        if (totalEnergy)
            *totalEnergy += recipEnergy;
        pme_destroy(pmedata);
    }

    // Ewald method

    else if (ewald) {

        // setup reciprocal box

        float recipBoxSize[3] = { TWO_PI / periodicBoxSize[0], TWO_PI / periodicBoxSize[1], TWO_PI / periodicBoxSize[2]};


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


void CpuNonbondedForce::calculateDirectIxn(int numberOfAtoms, float* posq, const vector<pair<float, float> >& atomParameters,
                const vector<set<int> >& exclusions, float* forces, float* totalEnergy) {
    // Record the parameters for the threads.
    
    this->posq = posq;
    this->atomParameters = atomParameters;
    this->exclusions = exclusions;
    includeEnergy = (totalEnergy != NULL);
    
    // Signal the threads to start running and wait for them to finish.
    
    pthread_mutex_lock(&lock);
    waitCount = 0;
    pthread_cond_broadcast(&startCondition);
    while (waitCount < numThreads)
        pthread_cond_wait(&endCondition, &lock);
    pthread_mutex_unlock(&lock);
    
    // Combine the results from all the threads.
    
    double directEnergy = 0;
    for (int i = 0; i < numThreads; i++)
        directEnergy += threadData[i]->threadEnergy;
    for (int i = 0; i < numberOfAtoms; i++) {
        __m128 f = _mm_loadu_ps(forces+4*i);
        for (int j = 0; j < numThreads; j++)
            f = _mm_add_ps(f, _mm_loadu_ps(&threadData[j]->threadForce[4*i]));
        _mm_storeu_ps(forces+4*i, f);
    }

    if (ewald || pme) {
        // Now subtract off the exclusions, since they were implicitly included in the reciprocal space sum.

        for (int i = 0; i < numberOfAtoms; i++)
            for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter) {
                if (*iter > i) {
                   int ii = i;
                   int jj = *iter;
                   __m128 deltaR;
                   __m128 posI = _mm_loadu_ps(posq+4*ii);
                   __m128 posJ = _mm_loadu_ps(posq+4*jj);
                   float r2;
                   getDeltaR(posJ, posI, deltaR, r2, false);
                   float r         = sqrtf(r2);
                   float inverseR  = 1/r;
                   float alphaR    = alphaEwald * r;
                   float erfAlphaR = 1.0f-erfcApprox(alphaR);
                   if (erfAlphaR > 1e-6) {
                       float chargeProd = ONE_4PI_EPS0*posq[4*ii+3]*posq[4*jj+3];
                       float dEdR      = (float) (chargeProd * inverseR * inverseR * inverseR);
                             dEdR      = (float) (dEdR * (erfAlphaR - TWO_OVER_SQRT_PI * alphaR * exp (- alphaR * alphaR)));
                       __m128 result = _mm_mul_ps(deltaR, _mm_set1_ps(dEdR));
                       _mm_storeu_ps(forces+4*ii, _mm_sub_ps(_mm_loadu_ps(forces+4*ii), result));
                       _mm_storeu_ps(forces+4*jj, _mm_add_ps(_mm_loadu_ps(forces+4*jj), result));
                       if (includeEnergy)
                           directEnergy -= chargeProd*inverseR*erfAlphaR;
                   }
                }
            }
    }
    if (totalEnergy != NULL)
        *totalEnergy += (float) directEnergy;
}


void CpuNonbondedForce::runThread(int index, vector<float>& threadForce, double& threadEnergy) {
    while (true) {
        // Wait for the signal to start running.
        
        pthread_mutex_lock(&lock);
        waitCount++;
        pthread_cond_signal(&endCondition);
        pthread_cond_wait(&startCondition, &lock);
        pthread_mutex_unlock(&lock);
        if (isDeleted)
            break;
        
        // Compute this thread's subset of interactions.
        
        threadEnergy = 0;
        double* energyPtr = (includeEnergy ? &threadEnergy : NULL);
        int numberOfAtoms = atomParameters.size();
        threadForce.resize(4*numberOfAtoms, 0.0f);
        for (int i = 0; i < 4*numberOfAtoms; i++)
            threadForce[i] = 0.0f;
        if (ewald || pme) {
            // Compute the interactions from the neighbor list.

            for (int i = index; i < (int) neighborList->size(); i += numThreads) {
                pair<int, int> pair = (*neighborList)[i];
                calculateOneEwaldIxn(pair.first, pair.second, &threadForce[0], energyPtr);
            }
        }
        else if (cutoff) {
            // Compute the interactions from the neighbor list.

            for (int i = index; i < (int) neighborList->size(); i += numThreads) {
                pair<int, int> pair = (*neighborList)[i];
                calculateOneIxn(pair.first, pair.second, &threadForce[0], energyPtr);
            }
        }
        else {
            // Loop over all atom pairs

            for (int i = index; i < numberOfAtoms; i += numThreads){
                for (int j = i+1; j < numberOfAtoms; j++)
                    if (exclusions[j].find(i) == exclusions[j].end())
                        calculateOneIxn(i, j, &threadForce[0], energyPtr);
            }
        }
    }
}

void CpuNonbondedForce::calculateOneIxn(int ii, int jj, float* forces, double* totalEnergy) {
    // get deltaR, R2, and R between 2 atoms

    __m128 deltaR;
    __m128 posI = _mm_loadu_ps(posq+4*ii);
    __m128 posJ = _mm_loadu_ps(posq+4*jj);
    float r2;
    getDeltaR(posJ, posI, deltaR, r2, periodic);
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
    dEdR     *= inverseR*inverseR;
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

    __m128 result = _mm_mul_ps(deltaR, _mm_set1_ps(dEdR));
    _mm_storeu_ps(forces+4*ii, _mm_add_ps(_mm_loadu_ps(forces+4*ii), result));
    _mm_storeu_ps(forces+4*jj, _mm_sub_ps(_mm_loadu_ps(forces+4*jj), result));
  }

void CpuNonbondedForce::calculateOneEwaldIxn(int ii, int jj, float* forces, double* totalEnergy) {
    __m128 deltaR;
    __m128 posI = _mm_loadu_ps(posq+4*ii);
    __m128 posJ = _mm_loadu_ps(posq+4*jj);
    float r2;
    getDeltaR(posJ, posI, deltaR, r2, true);
    if (r2 >= cutoffDistance*cutoffDistance)
        return;
    float r         = sqrtf(r2);
    float inverseR  = 1/r;
    float switchValue = 1, switchDeriv = 0;
    if (useSwitch && r > switchingDistance) {
        float t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
        switchValue = 1+t*t*t*(-10+t*(15-t*6));
        switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
    }
    float alphaR    = alphaEwald * r;
    float erfcAlphaR = erfcApprox(alphaR);
    float chargeProd = ONE_4PI_EPS0*posq[4*ii+3]*posq[4*jj+3];
    float dEdR      = (float) (chargeProd * inverseR * inverseR * inverseR);
          dEdR      = (float) (dEdR * (erfcAlphaR + TWO_OVER_SQRT_PI * alphaR * exp(-alphaR*alphaR)));

    float sig       = atomParameters[ii].first +  atomParameters[jj].first;
    float sig2      = inverseR*sig;
          sig2     *= sig2;
    float sig6      = sig2*sig2*sig2;
    float eps       = atomParameters[ii].second*atomParameters[jj].second;
          dEdR     += switchValue*eps*(12.0f*sig6 - 6.0f)*sig6*inverseR*inverseR;
    float energy = eps*(sig6-1.0f)*sig6;
    if (useSwitch) {
        dEdR -= energy*switchDeriv*inverseR;
        energy *= switchValue;
    }

    // accumulate forces

    __m128 result = _mm_mul_ps(deltaR, _mm_set1_ps(dEdR));
    _mm_storeu_ps(forces+4*ii, _mm_add_ps(_mm_loadu_ps(forces+4*ii), result));
    _mm_storeu_ps(forces+4*jj, _mm_sub_ps(_mm_loadu_ps(forces+4*jj), result));

    // accumulate energies

    if (totalEnergy) {
        energy += (float) (chargeProd*inverseR*erfcAlphaR);
        *totalEnergy += energy;
    }
}

void CpuNonbondedForce::getDeltaR(const __m128& posI, const __m128& posJ, __m128& deltaR, float& r2, bool periodic) const {
    deltaR = _mm_sub_ps(posJ, posI);
    if (periodic) {
        __m128 base = _mm_mul_ps(_mm_floor_ps(_mm_add_ps(_mm_mul_ps(deltaR, invBoxSize), half)), boxSize);
        deltaR = _mm_sub_ps(deltaR, base);
    }
    r2 = _mm_cvtss_f32(_mm_dp_ps(deltaR, deltaR, 0x71));
}

float CpuNonbondedForce::erfcApprox(float x) {
    // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
    // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
    // error of 3e-7.

    float t = 1.0f+(0.0705230784f+(0.0422820123f+(0.0092705272f+(0.0001520143f+(0.0002765672f+0.0000430638f*x)*x)*x)*x)*x)*x;
    t *= t;
    t *= t;
    t *= t;
    return 1.0f/(t*t);
}
