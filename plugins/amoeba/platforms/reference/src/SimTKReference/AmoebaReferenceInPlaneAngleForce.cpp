
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

RealOpenMM AmoebaReferenceInPlaneAngleForce::getPrefactorsGivenAngleCosine( RealOpenMM cosine,
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

   // static const std::string methodName = "AmoebaReferenceInPlaneAngleForce::getPrefactorsGivenAngleCosine";

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

RealOpenMM AmoebaReferenceInPlaneAngleForce::calculateAngleIxn( const RealVec& positionAtomA, const RealVec& positionAtomB,
                                                                        const RealVec& positionAtomC, const RealVec& positionAtomD,
                                                                        RealOpenMM angle,                      RealOpenMM angleK,
                                                                        RealOpenMM angleCubic,                 RealOpenMM angleQuartic,
                                                                        RealOpenMM anglePentic,                RealOpenMM angleSextic,
                                                                        RealVec* forces ) const {

   // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "AmoebaReferenceInPlaneAngleForce::calculateAngleIxn";
 
    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM one           = 1.0;
    static const RealOpenMM two           = 2.0;

   // ---------------------------------------------------------------------------------------

    // T   = AD x CD
    // P   = B + T*delta
    // AP  = A - P
    // CP  = A - P
    // M   = CP x AP

    enum { AD, BD, CD, T, AP, P, CP, M, APxM, CPxM, ADxBD, BDxCD, TxCD, ADxT, dBxAD, CDxdB, LastDeltaAtomIndex };

    std::vector<RealOpenMM> deltaR[LastDeltaAtomIndex];
    for( int ii = 0; ii < LastDeltaAtomIndex; ii++ ){
        deltaR[ii].resize(3);
    }
    AmoebaReferenceForce::loadDeltaR( positionAtomD, positionAtomA, deltaR[AD] );
    AmoebaReferenceForce::loadDeltaR( positionAtomD, positionAtomB, deltaR[BD] );
    AmoebaReferenceForce::loadDeltaR( positionAtomD, positionAtomC, deltaR[CD] );
 
    AmoebaReferenceForce::getCrossProduct( deltaR[AD], deltaR[CD], deltaR[T] );

    RealOpenMM rT2     = AmoebaReferenceForce::getNormSquared3( deltaR[T] );
    RealOpenMM delta   = AmoebaReferenceForce::getDotProduct3( deltaR[T], deltaR[BD] )/rT2;
         delta        *= -one;
 
    for( int ii = 0; ii < 3; ii++ ){
       deltaR[P][ii]  = positionAtomB[ii] + deltaR[T][ii]*delta;
       deltaR[AP][ii] = positionAtomA[ii] - deltaR[P][ii];
       deltaR[CP][ii] = positionAtomC[ii] - deltaR[P][ii];
    }
 
    RealOpenMM rAp2 = AmoebaReferenceForce::getNormSquared3( deltaR[AP] );
    RealOpenMM rCp2 = AmoebaReferenceForce::getNormSquared3( deltaR[CP] );
    if( rAp2 <= zero && rCp2 <= zero ){
       return zero;
    }
 
    AmoebaReferenceForce::getCrossProduct( deltaR[CP], deltaR[AP], deltaR[M] );
 
    RealOpenMM rm = AmoebaReferenceForce::getNorm3( deltaR[M] );
    if( rm < 1.0e-06 ){
       rm = 1.0e-06;
    }
 
    RealOpenMM dot     = AmoebaReferenceForce::getDotProduct3( deltaR[AP], deltaR[CP] );
    RealOpenMM cosine  = dot/SQRT( rAp2*rCp2 );
 
    RealOpenMM dEdR;
    RealOpenMM energy = getPrefactorsGivenAngleCosine( cosine, angle, angleK, angleCubic, angleQuartic,
                                                       anglePentic, angleSextic, &dEdR );

    RealOpenMM termA   = -dEdR/(rAp2*rm);
    RealOpenMM termC   =  dEdR/(rCp2*rm);
 
    AmoebaReferenceForce::getCrossProduct( deltaR[AP], deltaR[M], deltaR[APxM] );
    AmoebaReferenceForce::getCrossProduct( deltaR[CP], deltaR[M], deltaR[CPxM] );
 
    // forces will be gathered here
 
    enum { dA, dB, dC, dD, LastDIndex };
    std::vector<RealOpenMM> forceTerm[LastDIndex];
    for( int ii = 0; ii < LastDIndex; ii++ ){
        forceTerm[ii].resize(3);
    }
 
    for( int ii = 0; ii < 3; ii++ ){
       forceTerm[dA][ii] = deltaR[APxM][ii]*termA;
       forceTerm[dC][ii] = deltaR[CPxM][ii]*termC;
       forceTerm[dB][ii] = -one*( forceTerm[dA][ii] + forceTerm[dC][ii] );
    }
 
    RealOpenMM pTrT2  = AmoebaReferenceForce::getDotProduct3( forceTerm[dB], deltaR[T] );
               pTrT2 /= rT2;
 
    AmoebaReferenceForce::getCrossProduct( deltaR[CD], forceTerm[dB], deltaR[CDxdB] );
    AmoebaReferenceForce::getCrossProduct( forceTerm[dB], deltaR[AD], deltaR[dBxAD] );
 
    if( FABS( pTrT2 ) > 1.0e-08 ){
 
       RealOpenMM delta2 = delta*two;
 
       AmoebaReferenceForce::getCrossProduct( deltaR[BD], deltaR[CD], deltaR[BDxCD] );
       AmoebaReferenceForce::getCrossProduct( deltaR[T],  deltaR[CD], deltaR[TxCD]  );
       AmoebaReferenceForce::getCrossProduct( deltaR[AD], deltaR[BD], deltaR[ADxBD] );
       AmoebaReferenceForce::getCrossProduct( deltaR[AD], deltaR[T],  deltaR[ADxT]  );
       for( int ii = 0; ii < 3; ii++ ){
    
          RealOpenMM term     = deltaR[BDxCD][ii] + delta2*deltaR[TxCD][ii];
          forceTerm[dA][ii]  += delta*deltaR[CDxdB][ii] + term*pTrT2;
    
               term           = deltaR[ADxBD][ii] + delta2*deltaR[ADxT][ii];
          forceTerm[dC][ii]  += delta*deltaR[dBxAD][ii] + term*pTrT2;
    
          forceTerm[dD][ii]  = -( forceTerm[dA][ii] + forceTerm[dB][ii] + forceTerm[dC][ii] );
       }
    } else {
       for( int ii = 0; ii < 3; ii++ ){
    
          forceTerm[dA][ii] += delta*deltaR[CDxdB][ii];
          forceTerm[dC][ii] += delta*deltaR[dBxAD][ii];
 
          forceTerm[dD][ii]  = -( forceTerm[dA][ii] + forceTerm[dB][ii] + forceTerm[dC][ii] );
       }
    }
 
    // accumulate forces
 
    for( int jj = 0; jj < 4; jj++ ){
        forces[jj][0] = forceTerm[jj][0];
        forces[jj][1] = forceTerm[jj][1];
        forces[jj][2] = forceTerm[jj][2];
    }
 
    return energy;

}

RealOpenMM AmoebaReferenceInPlaneAngleForce::calculateForceAndEnergy( int numAngles, vector<RealVec>& posData,
                                                                       const std::vector<int>&  particle1,
                                                                       const std::vector<int>&  particle2,
                                                                       const std::vector<int>&  particle3,
                                                                       const std::vector<int>&  particle4,
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
        int particle4Index      = particle4[ii];
        RealOpenMM idealAngle   = angle[ii];
        RealOpenMM angleK       = kQuadratic[ii];
        RealVec forces[4];
        energy                 += calculateAngleIxn( posData[particle1Index], posData[particle2Index], posData[particle3Index], posData[particle4Index],
                                                     idealAngle, angleK, angleCubic, angleQuartic, anglePentic, angleSextic, forces );

        // accumulate forces
     
        for( int jj = 0; jj < 3; jj++ ){
            forceData[particle1Index][jj] -= forces[0][jj];
            forceData[particle2Index][jj] -= forces[1][jj];
            forceData[particle3Index][jj] -= forces[2][jj];
            forceData[particle4Index][jj] -= forces[3][jj];
        }
 
    }   
    return energy;
}

