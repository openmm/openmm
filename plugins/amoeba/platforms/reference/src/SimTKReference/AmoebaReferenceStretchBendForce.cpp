
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
#include "AmoebaReferenceStretchBendForce.h"
#include <vector>

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   Calculate Amoeba stretch bend angle ixn (force and energy)

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

RealOpenMM AmoebaReferenceStretchBendForce::calculateStretchBendIxn( const RealVec& positionAtomA, const RealVec& positionAtomB,
                                                                     const RealVec& positionAtomC,
                                                                     RealOpenMM lengthAB,      RealOpenMM lengthCB,
                                                                     RealOpenMM idealAngle,    RealOpenMM kParameter,
                                                                     RealVec* forces ) const {

   // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "AmoebaReferenceStretchBendForce::calculateStretchBendIxn";
 
    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM one           = 1.0;

   // ---------------------------------------------------------------------------------------

   enum { A, B, C, LastAtomIndex };
   enum { AB, CB, CBxAB, ABxP, CBxP, LastDeltaIndex };

   // ---------------------------------------------------------------------------------------

   // get deltaR between various combinations of the 3 atoms
   // and various intermediate terms

    std::vector<RealOpenMM> deltaR[LastDeltaIndex];
    for( unsigned int ii = 0; ii < LastDeltaIndex; ii++ ){
        deltaR[ii].resize(3);
    }
    AmoebaReferenceForce::loadDeltaR( positionAtomB, positionAtomA, deltaR[AB] );
    AmoebaReferenceForce::loadDeltaR( positionAtomB, positionAtomC, deltaR[CB] );
    RealOpenMM  rAB2 = AmoebaReferenceForce::getNormSquared3( deltaR[AB] );
    RealOpenMM  rAB  = SQRT( rAB2 );
    RealOpenMM  rCB2 = AmoebaReferenceForce::getNormSquared3( deltaR[CB] );
    RealOpenMM  rCB  = SQRT( rCB2 );

    AmoebaReferenceForce::getCrossProduct( deltaR[CB], deltaR[AB], deltaR[CBxAB] );
    RealOpenMM  rP   = AmoebaReferenceForce::getNorm3( deltaR[CBxAB] );
    if( rP <= zero ){
       return zero;
    }
    RealOpenMM dot    = AmoebaReferenceForce::getDotProduct3( deltaR[CB], deltaR[AB] );
    RealOpenMM cosine = dot/(rAB*rCB);
 
    RealOpenMM angle;
    if( cosine >= one ){
       angle = zero;
    } else if( cosine <= -one ){
       angle = PI_M;
    } else {
       angle = RADIAN*ACOS(cosine);
    }
 
    RealOpenMM termA = -RADIAN/(rAB2*rP);
    RealOpenMM termC =  RADIAN/(rCB2*rP);
 
    // P = CBxAB
 
    AmoebaReferenceForce::getCrossProduct( deltaR[AB], deltaR[CBxAB], deltaR[ABxP] );
    AmoebaReferenceForce::getCrossProduct( deltaR[CB], deltaR[CBxAB], deltaR[CBxP] );
    for( int ii = 0; ii < 3; ii++ ){
       deltaR[ABxP][ii] *= termA;
       deltaR[CBxP][ii] *= termC;
    }
 
    RealOpenMM dr    = rAB - lengthAB + rCB - lengthCB;
    termA            = one/rAB;
    termC            = one/rCB;
 
    // ---------------------------------------------------------------------------------------
 
    // forces
 
    // calculate forces for atoms a, b, c
    // the force for b is then -( a + c)
 
    std::vector<RealOpenMM> subForce[LastAtomIndex];
    for( int ii = 0; ii < LastAtomIndex; ii++ ){
        subForce[ii].resize(3);
    }
    RealOpenMM dt = angle - idealAngle*RADIAN;
    for( int jj = 0; jj < 3; jj++ ){
       subForce[A][jj] = kParameter*(dt*termA*deltaR[AB][jj] + dr*deltaR[ABxP][jj] );
       subForce[C][jj] = kParameter*(dt*termC*deltaR[CB][jj] + dr*deltaR[CBxP][jj] );
       subForce[B][jj] = -( subForce[A][jj] + subForce[C][jj] );
    }
 
    // add in forces
 
    for( int jj = 0; jj < LastAtomIndex; jj++ ){
        forces[jj][0] = subForce[jj][0];
        forces[jj][1] = subForce[jj][1];
        forces[jj][2] = subForce[jj][2];
    }
 
    // ---------------------------------------------------------------------------------------
 
    return (kParameter*dt*dr);
}

RealOpenMM AmoebaReferenceStretchBendForce::calculateForceAndEnergy( int numStretchBends, vector<RealVec>& posData,
                                                                       const std::vector<int>&  particle1,
                                                                       const std::vector<int>&  particle2,
                                                                       const std::vector<int>&  particle3,
                                                                       const std::vector<RealOpenMM>& lengthABParameters,
                                                                       const std::vector<RealOpenMM>& lengthCBParameters,
                                                                       const std::vector<RealOpenMM>&  angle,
                                                                       const std::vector<RealOpenMM>&  kQuadratic,
                                                                       vector<RealVec>& forceData) const {
    RealOpenMM energy      = 0.0; 
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numStretchBends); ii++) {
        int particle1Index      = particle1[ii];
        int particle2Index      = particle2[ii];
        int particle3Index      = particle3[ii];
        RealOpenMM abLength     = lengthABParameters[ii];
        RealOpenMM cbLength     = lengthCBParameters[ii];
        RealOpenMM idealAngle   = angle[ii];
        RealOpenMM angleK       = kQuadratic[ii];
        RealVec forces[3];
        energy                 += calculateStretchBendIxn( posData[particle1Index], posData[particle2Index], posData[particle3Index],
                                                           abLength, cbLength, idealAngle, angleK, forces );
        // accumulate forces
    
        for( int jj = 0; jj < 3; jj++ ){
            forceData[particle1Index][jj] -= forces[0][jj];
            forceData[particle2Index][jj] -= forces[1][jj];
            forceData[particle3Index][jj] -= forces[2][jj];
        }

    }   
    return energy;
}

