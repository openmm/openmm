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
#include "AmoebaReferenceTorsionForce.h"

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   Calculate Amoeba torsion ixn (force and energy)

   @param positionAtomA           Cartesian coordinates of atom A
   @param positionAtomB           Cartesian coordinates of atom B
   @param positionAtomC           Cartesian coordinates of atom C
   @param positionAtomD           Cartesian coordinates of atom D
   @param torsion1                vector of torsion params for first index (amplitude, phase, fold)
   @param torsion2                vector of torsion params for second index (amplitude, phase, fold)
   @param torsion3                vector of torsion params for third index (amplitude, phase, fold)
   @param forces                  force vector

   @return energy

   --------------------------------------------------------------------------------------- */

RealOpenMM AmoebaReferenceTorsionForce::calculateTorsionIxn( const RealVec& positionAtomA, const RealVec& positionAtomB,
                                                             const RealVec& positionAtomC, const RealVec& positionAtomD,
                                                             const std::vector<RealOpenMM>& torsionParameters1,
                                                             const std::vector<RealOpenMM>& torsionParameters2,
                                                             const std::vector<RealOpenMM>& torsionParameters3,
                                                             RealVec* forces ) const {

   // ---------------------------------------------------------------------------------------

   //static const std::string methodName = "AmoebaReferenceTorsionForce::calculateHarmonicForce";

   static const RealOpenMM zero          = 0.0;
   static const RealOpenMM one           = 1.0;
   static const RealOpenMM two           = 2.0;
   static const RealOpenMM three         = 3.0;

   // ---------------------------------------------------------------------------------------

    enum { BA, CB, DC, CA, DB, LastDeltaIndex };
    std::vector<RealOpenMM> deltaR[LastDeltaIndex];
    AmoebaReferenceForce::loadDeltaR( positionAtomA, positionAtomB, deltaR[BA] );
    AmoebaReferenceForce::loadDeltaR( positionAtomB, positionAtomC, deltaR[CB] );
    AmoebaReferenceForce::loadDeltaR( positionAtomC, positionAtomD, deltaR[DC] );
    AmoebaReferenceForce::loadDeltaR( positionAtomA, positionAtomC, deltaR[CA] );
    AmoebaReferenceForce::loadDeltaR( positionAtomB, positionAtomD, deltaR[DB] );

    enum { Xt, Xu, Xtu, LastXtIndex };
    std::vector<RealOpenMM> crossProducts[LastXtIndex];
    for( unsigned int ii = 0; ii < LastXtIndex; ii++ ){
        crossProducts[ii].resize(3);
    }
    AmoebaReferenceForce::getCrossProduct( deltaR[BA], deltaR[CB], crossProducts[Xt] );
    AmoebaReferenceForce::getCrossProduct( deltaR[CB], deltaR[DC], crossProducts[Xu] );
    AmoebaReferenceForce::getCrossProduct( crossProducts[Xt], crossProducts[Xu], crossProducts[Xtu] );
 
    RealOpenMM rT2   = AmoebaReferenceForce::getNormSquared3( crossProducts[Xt] );
    RealOpenMM rU2   = AmoebaReferenceForce::getNormSquared3( crossProducts[Xu] );
    RealOpenMM rTrU  = SQRT( rT2*rU2 );
    if( rTrU <= zero ){
       return zero;
    }
 
    RealOpenMM rCB   = AmoebaReferenceForce::getNorm3( deltaR[CB] );
 
    // ---------------------------------------------------------------------------------------
 
    // cos(w), cos(2w), cos(3w), ... 
    // sin(w), sin(2w), sin(3w), ... 
  
    RealOpenMM cosine           = AmoebaReferenceForce::getDotProduct3( crossProducts[Xt], crossProducts[Xu] );
               cosine          /= rTrU;
 
    RealOpenMM sine             = AmoebaReferenceForce::getDotProduct3( deltaR[CB], crossProducts[Xtu] );
               sine            /= (rCB*rTrU);

    RealOpenMM cosine2          = cosine*cosine - sine*sine;
    RealOpenMM sine2            = two*cosine*sine;
    RealOpenMM cosine3          = cosine*cosine2 - sine*sine2;
    RealOpenMM sine3            = cosine*sine2 + sine*cosine2;

/*
    RealOpenMM cosine4          = cosine*cosine3 - sine*sine3;
    RealOpenMM sine4            = cosine*sine3 + sine*cosine3;
    RealOpenMM cosine5          = cosine*cosine4 - sine*sine4;
    RealOpenMM sine5            = cosine*sine4 + sine*cosine4;
    RealOpenMM cosine6          = cosine*cosine5 - sine*sine5;
    RealOpenMM sine6            = cosine*sine5 + sine*cosine5;
*/
 
    // ---------------------------------------------------------------------------------------
 
    // dEdPhi prefactor
  
    RealOpenMM dEdPhi  = torsionParameters1[0]*      (cosine* sin( torsionParameters1[1] ) - sine* cos( torsionParameters1[1] ) );
    dEdPhi            += torsionParameters2[0]*two*  (cosine2*sin( torsionParameters2[1] ) - sine2*cos( torsionParameters2[1] ) );
    dEdPhi            += torsionParameters3[0]*three*(cosine3*sin( torsionParameters3[1] ) - sine3*cos( torsionParameters3[1] ) );

    // ---------------------------------------------------------------------------------------
  
    // dEdtu[0]      = dEdT
    // dEdtu[1]      = dEdU
  
    // tempVector[0] == dEdA: dEdT x CB
    // tempVector[1] == dEdB: (CA x dEdT) + (dEdU x DC)
    // tempVector[2] == dEdC: (dEdT x BA) + (DB x dEdU)
    // tempVector[3] == dEdD: (dEdU x CB)
   
    std::vector<RealOpenMM> dEdtu[2];
    std::vector<RealOpenMM> tempVector[6];
    for( unsigned int ii = 0; ii < 2; ii++ ){
        dEdtu[ii].resize(3);
    }
    for( unsigned int ii = 0; ii < 6; ii++ ){
        tempVector[ii].resize(3);
    }
 
    // dEdT & dEdU
  
    AmoebaReferenceForce::getCrossProduct( crossProducts[Xt], deltaR[CB], tempVector[0] );
    AmoebaReferenceForce::getCrossProduct( crossProducts[Xu], deltaR[CB], tempVector[1] );
    RealOpenMM norm[2] = { dEdPhi/(rT2*rCB ), -dEdPhi/(rU2*rCB ) };

    dEdtu[0][0] = norm[0]*tempVector[0][0];
    dEdtu[0][1] = norm[0]*tempVector[0][1];
    dEdtu[0][2] = norm[0]*tempVector[0][2];

    dEdtu[1][0] = norm[1]*tempVector[1][0];
    dEdtu[1][1] = norm[1]*tempVector[1][1];
    dEdtu[1][2] = norm[1]*tempVector[1][2];
 
    // dEdA
   
    AmoebaReferenceForce::getCrossProduct( dEdtu[0], deltaR[CB], tempVector[0] );
 
    // dEdB
  
    AmoebaReferenceForce::getCrossProduct( deltaR[CA], dEdtu[0], tempVector[4] );
    AmoebaReferenceForce::getCrossProduct( dEdtu[1], deltaR[DC], tempVector[1] );
 
    // dEdC
 
    AmoebaReferenceForce::getCrossProduct( dEdtu[0], deltaR[BA], tempVector[5] );
    AmoebaReferenceForce::getCrossProduct( deltaR[DB], dEdtu[1], tempVector[2] );

    tempVector[1][0] += tempVector[4][0];
    tempVector[2][0] += tempVector[5][0];
 
    tempVector[1][1] += tempVector[4][1];
    tempVector[2][1] += tempVector[5][1];
 
    tempVector[1][2] += tempVector[4][2];
    tempVector[2][2] += tempVector[5][2];
 
    // dEdD
  
    AmoebaReferenceForce::getCrossProduct( dEdtu[1], deltaR[CB], tempVector[3] );
 
    // ---------------------------------------------------------------------------------------
     
    // forces
  
    forces[0][0] -= tempVector[0][0];
    forces[0][1] -= tempVector[0][1];
    forces[0][2] -= tempVector[0][2];
 
    forces[1][0] -= tempVector[1][0];
    forces[1][1] -= tempVector[1][1];
    forces[1][2] -= tempVector[1][2];
 
    forces[2][0] -= tempVector[2][0];
    forces[2][1] -= tempVector[2][1];
    forces[2][2] -= tempVector[2][2];
 
    forces[3][0] -= tempVector[3][0];
    forces[3][1] -= tempVector[3][1];
    forces[3][2] -= tempVector[3][2];
 
    // ---------------------------------------------------------------------------------------
 
    // calculate energy
  
    RealOpenMM energy  = torsionParameters1[0]*( one + cosine* cos( torsionParameters1[1] ) + sine *sin( torsionParameters1[1] ) );
               energy += torsionParameters2[0]*( one + cosine2*cos( torsionParameters2[1] ) + sine2*sin( torsionParameters2[1] ) );
               energy += torsionParameters3[0]*( one + cosine3*cos( torsionParameters3[1] ) + sine3*sin( torsionParameters3[1] ) );

    return energy;
}

RealOpenMM AmoebaReferenceTorsionForce::calculateForceAndEnergy( int numTorsions, vector<RealVec>& posData,
                                                                 const std::vector<int>&  particle1,
                                                                 const std::vector<int>&  particle2,
                                                                 const std::vector<int>&  particle3,
                                                                 const std::vector<int>&  particle4,
                                                                 const std::vector< std::vector<RealOpenMM> >& torsionParameters1,
                                                                 const std::vector< std::vector<RealOpenMM> >& torsionParameters2,
                                                                 const std::vector< std::vector<RealOpenMM> >& torsionParameters3,
                                                                 vector<RealVec>& forceData ) const {
    RealOpenMM energy      = 0.0; 
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numTorsions); ii++) {
        int particle1Index      = particle1[ii];
        int particle2Index      = particle2[ii];
        int particle3Index      = particle3[ii];
        int particle4Index      = particle4[ii];
        RealVec forces[4];
        forces[0]               = forceData[particle1Index];
        forces[1]               = forceData[particle2Index];
        forces[2]               = forceData[particle3Index];
        forces[3]               = forceData[particle4Index];
        energy                 += calculateTorsionIxn( posData[particle1Index], posData[particle2Index], posData[particle3Index], posData[particle4Index],
                                                       torsionParameters1[ii], torsionParameters2[ii], torsionParameters3[ii], forces );
    }   
    return energy;
}


