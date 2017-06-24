
/* Portions copyright (c) 2006 Stanford University and Simbios.
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

#include "AmoebaReferenceForce.h"
#include "AmoebaReferenceInPlaneAngleForce.h"
#include "ReferenceForce.h"
#include "SimTKOpenMMRealType.h"

using std::vector;
using namespace OpenMM;

void AmoebaReferenceInPlaneAngleForce::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

/**---------------------------------------------------------------------------------------

   Get dEdT and energy prefactor given cosine of angle :: the calculation for different
   the force types is identical 

   @param  cosine               cosine of angle
   @param  idealAngle           ideal angle
   @param  angleK               angle k (quadratic prefactor)
   @param  angleCubic           cubic prefactor
   @param  angleQuartic         quartiic prefactor
   @param  anglePentic          pentic prefactor
   @param  angleSextic          sextic prefactor
   @param  dEdR                 dEdR

   @return energy

   --------------------------------------------------------------------------------------- */

double AmoebaReferenceInPlaneAngleForce::getPrefactorsGivenAngleCosine(double cosine,
                                                                       double idealAngle,     double angleK,
                                                                       double angleCubic,     double angleQuartic,
                                                                       double anglePentic,    double angleSextic,
                                                                       double* dEdR) const {
   double angle;
   if (cosine >= 1.0) {
      angle = 0.0;
   } else if (cosine <= -1.0) {
      angle = RADIAN*PI_M;
   } else {
      angle = RADIAN*ACOS(cosine);
   }   
   double deltaIdeal         = angle - idealAngle;
   double deltaIdeal2        = deltaIdeal*deltaIdeal;
   double deltaIdeal3        = deltaIdeal*deltaIdeal2;
   double deltaIdeal4        = deltaIdeal2*deltaIdeal2;

   *dEdR                     = (2.0 + 3.0*angleCubic*deltaIdeal +
                                4.0*angleQuartic*deltaIdeal2     + 
                                5.0*anglePentic*deltaIdeal3      +
                                6.0*angleSextic*deltaIdeal4);

   *dEdR                    *= RADIAN*angleK*deltaIdeal;

   double energy             = 1.0f + angleCubic*deltaIdeal + angleQuartic*deltaIdeal2 +
                               anglePentic*deltaIdeal3 + angleSextic*deltaIdeal4;
   energy                   *= angleK*deltaIdeal2;

   return energy;

}

/**---------------------------------------------------------------------------------------

   Calculate Amoeba angle ixn (force and energy)

   @param positionAtomA           Cartesian coordinates of atom A
   @param positionAtomB           Cartesian coordinates of atom B
   @param positionAtomC           Cartesian coordinates of atom C
   @param positionAtomD           Cartesian coordinates of atom D
   @param angleLength             angle
   @param angleK                  quadratic angle force
   @param angleCubic              cubic angle force parameter
   @param angleQuartic            quartic angle force parameter
   @param anglePentic             pentic angle force parameter
   @param angleSextic             sextic angle force parameter
   @param forces                  force vector

   @return energy

   --------------------------------------------------------------------------------------- */

double AmoebaReferenceInPlaneAngleForce::calculateAngleIxn(const Vec3& positionAtomA, const Vec3& positionAtomB,
                                                           const Vec3& positionAtomC, const Vec3& positionAtomD,
                                                           double angle,                      double angleK,
                                                           double angleCubic,                 double angleQuartic,
                                                           double anglePentic,                double angleSextic,
                                                           Vec3* forces) const {
    // T   = AD x CD
    // P   = B + T*delta
    // AP  = A - P
    // CP  = A - P
    // M   = CP x AP

    enum { AD, BD, CD, T, AP, P, CP, M, APxM, CPxM, ADxBD, BDxCD, TxCD, ADxT, dBxAD, CDxdB, LastDeltaAtomIndex };

    Vec3 deltaR[LastDeltaAtomIndex];
    if (usePeriodic) {
        deltaR[AD] = ReferenceForce::getDeltaRPeriodic(positionAtomD, positionAtomA, boxVectors);
        deltaR[BD] = ReferenceForce::getDeltaRPeriodic(positionAtomD, positionAtomB, boxVectors);
        deltaR[CD] = ReferenceForce::getDeltaRPeriodic(positionAtomD, positionAtomC, boxVectors);
    }
    else {
        deltaR[AD] = ReferenceForce::getDeltaR(positionAtomD, positionAtomA);
        deltaR[BD] = ReferenceForce::getDeltaR(positionAtomD, positionAtomB);
        deltaR[CD] = ReferenceForce::getDeltaR(positionAtomD, positionAtomC);
    }
 
    deltaR[T] = deltaR[AD].cross(deltaR[CD]);

    double rT2 = deltaR[T].dot(deltaR[T]);
    double delta = -deltaR[T].dot(deltaR[BD])/rT2;

    deltaR[P]  = positionAtomB + deltaR[T]*delta;
    if (usePeriodic) {
        deltaR[AP] = ReferenceForce::getDeltaRPeriodic(deltaR[P], positionAtomA, boxVectors);
        deltaR[CP] = ReferenceForce::getDeltaRPeriodic(deltaR[P], positionAtomC, boxVectors);
    }
    else {
        deltaR[AP] = ReferenceForce::getDeltaR(deltaR[P], positionAtomA);
        deltaR[CP] = ReferenceForce::getDeltaR(deltaR[P], positionAtomC);
    }
 
    double rAp2 = deltaR[AP].dot(deltaR[AP]);
    double rCp2 = deltaR[CP].dot(deltaR[CP]);
    if (rAp2 <= 0 && rCp2 <= 0) {
       return 0;
    }
 
    deltaR[M] = deltaR[CP].cross(deltaR[AP]);
 
    double rm = sqrt(deltaR[M].dot(deltaR[M]));
    if (rm < 1.0e-06) {
       rm = 1.0e-06;
    }
 
    double dot     = deltaR[AP].dot(deltaR[CP]);
    double cosine  = dot/SQRT(rAp2*rCp2);
 
    double dEdR;
    double energy = getPrefactorsGivenAngleCosine(cosine, angle, angleK, angleCubic, angleQuartic,
                                                      anglePentic, angleSextic, &dEdR);

    double termA   = -dEdR/(rAp2*rm);
    double termC   =  dEdR/(rCp2*rm);
 
    deltaR[APxM] = deltaR[AP].cross(deltaR[M]);
    deltaR[CPxM] = deltaR[CP].cross(deltaR[M]);
 
    // forces will be gathered here
 
    enum { dA, dB, dC, dD, LastDIndex };
    Vec3 forceTerm[LastDIndex];
 
    forceTerm[dA] = deltaR[APxM]*termA;
    forceTerm[dC] = deltaR[CPxM]*termC;
    forceTerm[dB] = -(forceTerm[dA] + forceTerm[dC]);

    double pTrT2  = forceTerm[dB].dot(deltaR[T]);
           pTrT2 /= rT2;
 
    deltaR[CDxdB] = deltaR[CD].cross(forceTerm[dB]);
    deltaR[dBxAD] = forceTerm[dB].cross(deltaR[AD]);
 
    if (FABS(pTrT2) > 1.0e-08) {
 
        double delta2 = delta*2;
 
        deltaR[BDxCD] = forceTerm[dB].cross(deltaR[CD]);
        deltaR[TxCD] = forceTerm[T].cross(deltaR[CD]);
        deltaR[ADxBD] = forceTerm[AD].cross(deltaR[BD]);
        deltaR[ADxT] = forceTerm[AD].cross(deltaR[T]);
        Vec3 term = deltaR[BDxCD] + deltaR[TxCD]*delta2;
        forceTerm[dA] += deltaR[CDxdB]*delta + term*pTrT2;
        term = deltaR[ADxBD] + deltaR[ADxT]*delta2;
        forceTerm[dC] += deltaR[dBxAD]*delta + term*pTrT2;
        forceTerm[dD] = -(forceTerm[dA] + forceTerm[dB] + forceTerm[dC]);
    } else {
        forceTerm[dA] += deltaR[CDxdB]*delta;
        forceTerm[dC] += deltaR[dBxAD]*delta;
        forceTerm[dD]  = -(forceTerm[dA] + forceTerm[dB] + forceTerm[dC]);
    }
 
    // accumulate forces
 
    for (int jj = 0; jj < 4; jj++)
        forces[jj] = forceTerm[jj];
 
    return energy;

}

double AmoebaReferenceInPlaneAngleForce::calculateForceAndEnergy(int numAngles, vector<Vec3>& posData,
                                                                 const std::vector<int>&  particle1,
                                                                 const std::vector<int>&  particle2,
                                                                 const std::vector<int>&  particle3,
                                                                 const std::vector<int>&  particle4,
                                                                 const std::vector<double>&  angle,
                                                                 const std::vector<double>&  kQuadratic,
                                                                 double angleCubic,
                                                                 double angleQuartic,
                                                                 double anglePentic,
                                                                 double angleSextic,
                                                                 vector<Vec3>& forceData) const {
    double energy      = 0.0; 
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numAngles); ii++) {
        int particle1Index = particle1[ii];
        int particle2Index = particle2[ii];
        int particle3Index = particle3[ii];
        int particle4Index = particle4[ii];
        double idealAngle  = angle[ii];
        double angleK      = kQuadratic[ii];
        Vec3 forces[4];
        energy            += calculateAngleIxn(posData[particle1Index], posData[particle2Index], posData[particle3Index], posData[particle4Index],
                                               idealAngle, angleK, angleCubic, angleQuartic, anglePentic, angleSextic, forces);

        // accumulate forces
     
        for (int jj = 0; jj < 3; jj++) {
            forceData[particle1Index][jj] -= forces[0][jj];
            forceData[particle2Index][jj] -= forces[1][jj];
            forceData[particle3Index][jj] -= forces[2][jj];
            forceData[particle4Index][jj] -= forces[3][jj];
        }
 
    }   
    return energy;
}

