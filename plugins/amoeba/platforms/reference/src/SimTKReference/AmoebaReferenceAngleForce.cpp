
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
#include "AmoebaReferenceAngleForce.h"

using std::vector;
using OpenMM::RealVec;

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

RealOpenMM AmoebaReferenceAngleForce::getPrefactorsGivenAngleCosine( RealOpenMM cosine,
                                                                             RealOpenMM idealAngle,     RealOpenMM angleK,
                                                                             RealOpenMM angleCubic,     RealOpenMM angleQuartic,
                                                                             RealOpenMM anglePentic,    RealOpenMM angleSextic,
                                                                             RealOpenMM* dEdR ) const {

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM zero          = 0.0;
   static const RealOpenMM one           = 1.0;
   static const RealOpenMM two           = 2.0;
   static const RealOpenMM three         = 3.0;
   static const RealOpenMM four          = 4.0;
   static const RealOpenMM five          = 5.0;
   static const RealOpenMM six           = 6.0;

   // static const std::string methodName = "AmoebaReferenceAngleForce::getPrefactorsGivenAngleCosine";

   // ---------------------------------------------------------------------------------------

   RealOpenMM angle;
   if( cosine >= one ){
      angle = zero;
   } else if( cosine <= -one ){
      angle = RADIAN*PI_M;
   } else {
      angle = RADIAN*ACOS(cosine);
   }   
   RealOpenMM deltaIdeal         = angle - idealAngle;
   RealOpenMM deltaIdeal2        = deltaIdeal*deltaIdeal;
   RealOpenMM deltaIdeal3        = deltaIdeal*deltaIdeal2;
   RealOpenMM deltaIdeal4        = deltaIdeal2*deltaIdeal2;

   *dEdR                         = ( two + three*angleCubic*deltaIdeal +
                                     four*angleQuartic*deltaIdeal2     + 
                                     five*anglePentic*deltaIdeal3      +
                                     six*angleSextic*deltaIdeal4 );

   *dEdR                        *= RADIAN*angleK*deltaIdeal;

   RealOpenMM energy             = 1.0f + angleCubic*deltaIdeal + angleQuartic*deltaIdeal2 +
                                   anglePentic*deltaIdeal3 + angleSextic*deltaIdeal4;
   energy                       *= angleK*deltaIdeal2;

   return energy;

}

/**---------------------------------------------------------------------------------------

   Calculate Amoeba angle ixn (force and energy)

   @param positionAtomA           Cartesian coordinates of atom A
   @param positionAtomB           Cartesian coordinates of atom B
   @param positionAtomC           Cartesian coordinates of atom C
   @param angleLength             angle
   @param angleK                  quadratic angle force
   @param angleCubic              cubic angle force parameter
   @param angleQuartic            quartic angle force parameter
   @param anglePentic             pentic angle force parameter
   @param angleSextic             sextic angle force parameter
   @param forces                  force vector

   @return energy

   --------------------------------------------------------------------------------------- */

RealOpenMM AmoebaReferenceAngleForce::calculateAngleIxn( const RealVec& positionAtomA, const RealVec& positionAtomB,
                                                                 const RealVec& positionAtomC,
                                                                 RealOpenMM angle,          RealOpenMM angleK,
                                                                 RealOpenMM angleCubic,     RealOpenMM angleQuartic,
                                                                 RealOpenMM anglePentic,    RealOpenMM angleSextic,
                                                                 RealVec* forces ) const {

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "AmoebaReferenceAngleForce::calculateAngleIxn";

   static const RealOpenMM zero          = 0.0;
   static const RealOpenMM one           = 1.0;
   static const RealOpenMM onePt5        = 1.5;
   static const RealOpenMM two           = 2.0;

   // ---------------------------------------------------------------------------------------

    std::vector<RealOpenMM> deltaR[2];
    AmoebaReferenceForce::loadDeltaR( positionAtomA, positionAtomB, deltaR[0] );
    RealOpenMM rAB2      = AmoebaReferenceForce::getNormSquared3( deltaR[0] );
    RealOpenMM rAB       = SQRT( rAB2 );
 
    AmoebaReferenceForce::loadDeltaR( positionAtomC, positionAtomB, deltaR[1] );
    RealOpenMM rCB2      = AmoebaReferenceForce::getNormSquared3( deltaR[1] );
    RealOpenMM rCB       = SQRT( rCB2 );
 
    if( rAB <= zero || rCB <= zero ){
       return zero;
    }
 
    std::vector<RealOpenMM> pVector(3);
    AmoebaReferenceForce::getCrossProduct( deltaR[0], deltaR[1], pVector );
    RealOpenMM rp = AmoebaReferenceForce::getNorm3( pVector );
    if( rp < 1.0e-06 ){
       rp = 1.0e-06;
    }
    RealOpenMM dot    = AmoebaReferenceForce::getDotProduct3( deltaR[0], deltaR[1] );
    RealOpenMM cosine = dot/(rAB*rCB);
 
    RealOpenMM dEdR;
    RealOpenMM energy = getPrefactorsGivenAngleCosine( cosine, angle, angleK, angleCubic, angleQuartic,
                                                       anglePentic, angleSextic, &dEdR );
 
    RealOpenMM termA  =  dEdR/(rAB2*rp);
    RealOpenMM termC  = -dEdR/(rCB2*rp);
 
    std::vector<RealOpenMM> deltaCrossP[3];
    deltaCrossP[0].resize(3); 
    deltaCrossP[1].resize(3); 
    deltaCrossP[2].resize(3); 
    AmoebaReferenceForce::getCrossProduct( deltaR[0], pVector, deltaCrossP[0] );
    AmoebaReferenceForce::getCrossProduct( deltaR[1], pVector, deltaCrossP[2] );
    for( unsigned int ii = 0; ii < 3; ii++ ){
       deltaCrossP[0][ii] *= termA;
       deltaCrossP[2][ii] *= termC;
       deltaCrossP[1][ii]  = -1.0f*(deltaCrossP[0][ii] + deltaCrossP[2][ii]);
    }
 
    // accumulate forces
 
    for( int jj = 0; jj < 3; jj++ ){
        forces[jj][0] = deltaCrossP[jj][0];
        forces[jj][1] = deltaCrossP[jj][1];
        forces[jj][2] = deltaCrossP[jj][2];
    }
 
    return energy;
}

RealOpenMM AmoebaReferenceAngleForce::calculateForceAndEnergy( int numAngles, vector<RealVec>& posData,
                                                                       const std::vector<int>&  particle1,
                                                                       const std::vector<int>&  particle2,
                                                                       const std::vector<int>&  particle3,
                                                                       const std::vector<RealOpenMM>&  angle,
                                                                       const std::vector<RealOpenMM>&  kQuadratic,
                                                                       RealOpenMM angleCubic,
                                                                       RealOpenMM angleQuartic,
                                                                       RealOpenMM anglePentic,
                                                                       RealOpenMM angleSextic,
                                                                       vector<RealVec>& forceData) const {
    RealOpenMM energy      = 0.0; 
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numAngles); ii++) {
        int particle1Index      = particle1[ii];
        int particle2Index      = particle2[ii];
        int particle3Index      = particle3[ii];
        RealOpenMM idealAngle   = angle[ii];
        RealOpenMM angleK       = kQuadratic[ii];
        RealVec forces[3];
        energy                 += calculateAngleIxn( posData[particle1Index], posData[particle2Index], posData[particle3Index],
                                                     idealAngle, angleK, angleCubic, angleQuartic, anglePentic, angleSextic, forces );

        for( unsigned int jj = 0; jj < 3; jj++ ){
            forceData[particle1Index][jj] += forces[0][jj];
            forceData[particle2Index][jj] += forces[1][jj];
            forceData[particle3Index][jj] += forces[2][jj];
        }
    }   
    return energy;
}
