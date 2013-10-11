
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
#include <sstream>
#include <complex>

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "CpuNonbondedForce.h"
#include "ReferenceForce.h"
#include "ReferencePME.h"

// In case we're using some primitive version of Visual Studio this will
// make sure that erf() and erfc() are defined.
#include "openmm/internal/MSVC_erfc.h"

using namespace std;

/**---------------------------------------------------------------------------------------

   CpuNonbondedForce constructor

   --------------------------------------------------------------------------------------- */

CpuNonbondedForce::CpuNonbondedForce() : cutoff(false), useSwitch(false), periodic(false), ewald(false), pme(false) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuNonbondedForce::CpuNonbondedForce";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   CpuNonbondedForce destructor

   --------------------------------------------------------------------------------------- */

CpuNonbondedForce::~CpuNonbondedForce(){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuNonbondedForce::~CpuNonbondedForce";

   // ---------------------------------------------------------------------------------------

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

/**---------------------------------------------------------------------------------------

   Calculate Ewald ixn

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomParameters   atom parameters                             atomParameters[atomIndex][paramterIndex]
   @param exclusions       atom exclusion indices
                           exclusions[atomIndex] contains the list of exclusions for that atom
   @param fixedParameters  non atom parameters (not currently used)
   @param forces           force array (forces added)
   @param totalEnergy      total energy
   @param includeDirect      true if direct space interactions should be included
   @param includeReciprocal  true if reciprocal space interactions should be included

   --------------------------------------------------------------------------------------- */

void CpuNonbondedForce::calculateEwaldIxn(int numberOfAtoms, float* atomCoordinates,
                                             float** atomParameters, vector<set<int> >& exclusions,
                                             float* fixedParameters, float* forces,
                                             float* totalEnergy, bool includeDirect, bool includeReciprocal) const {
    typedef std::complex<float> d_complex;

    static const float epsilon     =  1.0;
    static const float one         =  1.0;
    static const float six         =  6.0;
    static const float twelve      = 12.0;

    int kmax                            = (ewald ? std::max(numRx, std::max(numRy,numRz)) : 0);
    float  factorEwald             = -1 / (4*alphaEwald*alphaEwald);
    float SQRT_PI                  = sqrt(PI_M);
    float TWO_PI                   = 2.0 * PI_M;
    float recipCoeff               = (float)(ONE_4PI_EPS0*4*PI_M/(periodicBoxSize[0] * periodicBoxSize[1] * periodicBoxSize[2]) /epsilon);

    float totalSelfEwaldEnergy     = 0.0;
    float realSpaceEwaldEnergy     = 0.0;
    float recipEnergy              = 0.0;
    float vdwEnergy                = 0.0;

// **************************************************************************************
// SELF ENERGY
// **************************************************************************************

    if (includeReciprocal) {
        for (int atomID = 0; atomID < numberOfAtoms; atomID++){
            float selfEwaldEnergy       = (float) (ONE_4PI_EPS0*atomCoordinates[4*atomID+3]*atomCoordinates[4*atomID+3] * alphaEwald/SQRT_PI);
            totalSelfEwaldEnergy            -= selfEwaldEnergy;
        }
    }

    if (totalEnergy){
        *totalEnergy += totalSelfEwaldEnergy;
    }

// **************************************************************************************
// RECIPROCAL SPACE EWALD ENERGY AND FORCES
// **************************************************************************************
    // PME

  if (pme && includeReciprocal) {
    pme_t          pmedata; /* abstract handle for PME data */
    float virial[3][3];

    pme_init(&pmedata,alphaEwald,numberOfAtoms,meshDim,5,1);

//    pme_exec(pmedata,atomCoordinates,forces,atomParameters,periodicBoxSize,&recipEnergy,virial);

    if (totalEnergy)
       *totalEnergy += recipEnergy;

        pme_destroy(pmedata);
  }

    // Ewald method

  else if (ewald && includeReciprocal) {

    // setup reciprocal box

    float recipBoxSize[3] = { TWO_PI / periodicBoxSize[0], TWO_PI / periodicBoxSize[1], TWO_PI / periodicBoxSize[2]};


    // setup K-vectors

  #define EIR(x, y, z) eir[(x)*numberOfAtoms*3+(y)*3+z]
  vector<d_complex> eir(kmax*numberOfAtoms*3);
  vector<d_complex> tab_xy(numberOfAtoms);
  vector<d_complex> tab_qxyz(numberOfAtoms);

  if (kmax < 1) {
      std::stringstream message;
      message << " kmax < 1 , Aborting" << std::endl;
      SimTKOpenMMLog::printError(message);
  }

  for (int i = 0; (i < numberOfAtoms); i++) {
    float* pos = atomCoordinates+4*i;
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
             tab_qxyz[n] = atomParameters[n][QIndex] * (tab_xy[n] * EIR(rz, n, 2));
          }

          else {
            for (int n = 0; n < numberOfAtoms; n++)
              tab_qxyz[n] = atomParameters[n][QIndex] * (tab_xy[n] * conj(EIR(-rz, n, 2)));
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
            float* f = forces+4*n;
            f[0] += 2 * recipCoeff * force * kx;
            f[1] += 2 * recipCoeff * force * ky;
            f[2] += 2 * recipCoeff * force * kz;
          }

          if (totalEnergy)
             *totalEnergy += recipCoeff * ak * (cs * cs + ss * ss);

          lowrz = 1 - numRz;
        }
        lowry = 1 - numRy;
      }
    }
  }

// **************************************************************************************
// SHORT-RANGE ENERGY AND FORCES
// **************************************************************************************

    if (!includeDirect)
        return;
    double totalVdwEnergy            = 0.0f;
    double totalRealSpaceEwaldEnergy = 0.0f;

    for (int i = 0; i < (int) neighborList->size(); i++) {
       pair<int, int> pair = (*neighborList)[i];
       int ii = pair.first;
       int jj = pair.second;

       __m128 deltaR;
       __m128 posI = _mm_loadu_ps(atomCoordinates+4*ii);
       __m128 posJ = _mm_loadu_ps(atomCoordinates+4*jj);
       float r2;
       getDeltaR(posJ, posI, deltaR, r2, true);
       float r         = sqrtf(r2);
       float inverseR  = 1/r;
       float switchValue = 1, switchDeriv = 0;
       if (useSwitch && r > switchingDistance) {
           float t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
           switchValue = 1+t*t*t*(-10+t*(15-t*6));
           switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
       }
       float alphaR    = alphaEwald * r;


       float chargeProd = ONE_4PI_EPS0*atomCoordinates[4*ii+3]*atomCoordinates[4*jj+3];
       float dEdR      = (float) (chargeProd * inverseR * inverseR * inverseR);
             dEdR      = (float) (dEdR * (erfc(alphaR) + 2 * alphaR * exp (- alphaR * alphaR) / SQRT_PI));

       float sig       = atomParameters[ii][SigIndex] +  atomParameters[jj][SigIndex];
       float sig2      = inverseR*sig;
                  sig2     *= sig2;
       float sig6      = sig2*sig2*sig2;
       float eps       = atomParameters[ii][EpsIndex]*atomParameters[jj][EpsIndex];
                  dEdR     += switchValue*eps*(twelve*sig6 - six)*sig6*inverseR*inverseR;
       vdwEnergy = eps*(sig6-one)*sig6;
       if (useSwitch) {
           dEdR -= vdwEnergy*switchDeriv*inverseR;
           vdwEnergy *= switchValue;
       }

       // accumulate forces

       __m128 result = _mm_mul_ps(deltaR, _mm_set1_ps(dEdR));
       _mm_storeu_ps(forces+4*ii, _mm_add_ps(_mm_loadu_ps(forces+4*ii), result));
       _mm_storeu_ps(forces+4*jj, _mm_sub_ps(_mm_loadu_ps(forces+4*jj), result));

       // accumulate energies

       realSpaceEwaldEnergy        = (float) (chargeProd*inverseR*erfc(alphaR));

       totalVdwEnergy             += vdwEnergy;
       totalRealSpaceEwaldEnergy  += realSpaceEwaldEnergy;

    }

    if (totalEnergy)
        *totalEnergy += totalRealSpaceEwaldEnergy + totalVdwEnergy;

    // Now subtract off the exclusions, since they were implicitly included in the reciprocal space sum.

    float totalExclusionEnergy = 0.0f;
    for (int i = 0; i < numberOfAtoms; i++)
        for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter) {
            if (*iter > i) {
               int ii = i;
               int jj = *iter;

               __m128 deltaR;
               __m128 posI = _mm_loadu_ps(atomCoordinates+4*ii);
               __m128 posJ = _mm_loadu_ps(atomCoordinates+4*jj);
               float r2;
               getDeltaR(posJ, posI, deltaR, r2, false);
               float r         = sqrtf(r2);
               float inverseR  = 1/r;
               float alphaR    = alphaEwald * r;
               if (erf(alphaR) > 1e-6) {
                   float chargeProd = ONE_4PI_EPS0*atomCoordinates[4*ii+3]*atomCoordinates[4*jj+3];
                   float dEdR      = (float) (chargeProd * inverseR * inverseR * inverseR);
                         dEdR      = (float) (dEdR * (erf(alphaR) - 2 * alphaR * exp (- alphaR * alphaR) / SQRT_PI));

                   // accumulate forces

                   __m128 result = _mm_mul_ps(deltaR, _mm_set1_ps(dEdR));
                   _mm_storeu_ps(forces+4*ii, _mm_add_ps(_mm_loadu_ps(forces+4*ii), result));
                   _mm_storeu_ps(forces+4*jj, _mm_sub_ps(_mm_loadu_ps(forces+4*jj), result));

                   // accumulate energies

                   realSpaceEwaldEnergy = (float) (chargeProd*inverseR*erf(alphaR));

                   totalExclusionEnergy += realSpaceEwaldEnergy;
               }
            }
        }

    if (totalEnergy)
        *totalEnergy -= totalExclusionEnergy;
}


/**---------------------------------------------------------------------------------------

   Calculate LJ Coulomb pair ixn

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomParameters   atom parameters                             atomParameters[atomIndex][paramterIndex]
   @param exclusions       atom exclusion indices
                           exclusions[atomIndex] contains the list of exclusions for that atom
   @param fixedParameters  non atom parameters (not currently used)
   @param forces           force array (forces added)
   @param totalEnergy      total energy
   @param includeDirect      true if direct space interactions should be included
   @param includeReciprocal  true if reciprocal space interactions should be included

   --------------------------------------------------------------------------------------- */

void CpuNonbondedForce::calculatePairIxn(int numberOfAtoms, float* atomCoordinates,
                                             float** atomParameters, vector<set<int> >& exclusions,
                                             float* fixedParameters, float* forces,
                                             float* totalEnergy, bool includeDirect, bool includeReciprocal) const {

   if (ewald || pme) {
       calculateEwaldIxn(numberOfAtoms, atomCoordinates, atomParameters, exclusions, fixedParameters, forces,
               totalEnergy, includeDirect, includeReciprocal);
       return;
   }
   if (!includeDirect)
       return;
   double directEnergy = 0;
   double* energyPtr = (totalEnergy == NULL ? NULL : &directEnergy);
   if (cutoff) {
       for (int i = 0; i < (int) neighborList->size(); i++) {
           pair<int, int> pair = (*neighborList)[i];
           calculateOneIxn(pair.first, pair.second, atomCoordinates, atomParameters, forces, energyPtr);
       }
   }
   else {
       for (int ii = 0; ii < numberOfAtoms; ii++){
          // loop over atom pairs

          for (int jj = ii+1; jj < numberOfAtoms; jj++)
              if (exclusions[jj].find(ii) == exclusions[jj].end())
                  calculateOneIxn(ii, jj, atomCoordinates, atomParameters, forces, energyPtr);
       }
   }
   if (totalEnergy != NULL)
       *totalEnergy += (float) directEnergy;
}

  /**---------------------------------------------------------------------------------------

     Calculate LJ Coulomb pair ixn between two atoms

     @param ii               the index of the first atom
     @param jj               the index of the second atom
     @param atomCoordinates  atom coordinates
     @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
     @param forces           force array (forces added)
     @param totalEnergy      total energy

     --------------------------------------------------------------------------------------- */

void CpuNonbondedForce::calculateOneIxn(int ii, int jj, float* atomCoordinates,
                        float** atomParameters, float* forces,
                        double* totalEnergy) const {
    // get deltaR, R2, and R between 2 atoms

    __m128 deltaR;
    __m128 posI = _mm_loadu_ps(atomCoordinates+4*ii);
    __m128 posJ = _mm_loadu_ps(atomCoordinates+4*jj);
    float r2;
    getDeltaR(posJ, posI, deltaR, r2, periodic);
    float r = sqrtf(r2);
    float inverseR = 1/r;
    float switchValue = 1, switchDeriv = 0;
    if (useSwitch && r > switchingDistance) {
        float t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
        switchValue = 1+t*t*t*(-10+t*(15-t*6));
        switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
    }
    float sig       = atomParameters[ii][SigIndex] +  atomParameters[jj][SigIndex];
    float sig2      = inverseR*sig;
          sig2     *= sig2;
    float sig6      = sig2*sig2*sig2;

    float eps       = atomParameters[ii][EpsIndex]*atomParameters[jj][EpsIndex];
    float dEdR      = switchValue*eps*(12.0f*sig6 - 6.0f)*sig6;
    float chargeProd = ONE_4PI_EPS0*atomCoordinates[4*ii+3]*atomCoordinates[4*jj+3];
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

void CpuNonbondedForce::getDeltaR(const __m128& posI, const __m128& posJ, __m128& deltaR, float& r2, bool periodic) const {
    deltaR = _mm_sub_ps(posJ, posI);
    if (periodic) {
        __m128 base = _mm_mul_ps(_mm_floor_ps(_mm_add_ps(_mm_mul_ps(deltaR, invBoxSize), half)), boxSize);
        deltaR = _mm_sub_ps(deltaR, base);
    }
    r2 = _mm_cvtss_f32(_mm_dp_ps(deltaR, deltaR, 0x71));
}
