
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
#include "AmoebaReferenceOutOfPlaneBendForce.h"

using std::vector;
using namespace OpenMM;

void AmoebaReferenceOutOfPlaneBendForce::setPeriodic(OpenMM::RealVec* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

/**---------------------------------------------------------------------------------------

   Calculate Amoeba Out-Of-Plane-Bend  ixn (force and energy)

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

RealOpenMM AmoebaReferenceOutOfPlaneBendForce::calculateOutOfPlaneBendIxn(const RealVec& positionAtomA, const RealVec& positionAtomB,
                                                                          const RealVec& positionAtomC, const RealVec& positionAtomD,
                                                                          RealOpenMM angleK,
                                                                          RealOpenMM angleCubic,                 RealOpenMM angleQuartic,
                                                                          RealOpenMM anglePentic,                RealOpenMM angleSextic,
                                                                          RealVec* forces) const {

   // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "AmoebaReferenceOutOfPlaneBendForce::calculateOutOfPlaneBendIxn";
 
    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM one           = 1.0;
    static const RealOpenMM two           = 2.0;
    static const RealOpenMM three         = 3.0;
    static const RealOpenMM four          = 4.0;
    static const RealOpenMM five          = 5.0;
    static const RealOpenMM six           = 6.0;

    enum { A, B, C, D, LastAtomIndex };
    enum { AB, CB, DB, AD, CD, LastDeltaIndex };
 
    // ---------------------------------------------------------------------------------------
 
    // get deltaR between various combinations of the 4 atoms
    // and various intermediate terms
 
    std::vector<RealOpenMM> deltaR[LastDeltaIndex];
    for (int ii = 0; ii < LastDeltaIndex; ii++) {
        deltaR[ii].resize(3);
    }
    if (usePeriodic) {
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomB, positionAtomA, deltaR[AB], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomB, positionAtomC, deltaR[CB], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomB, positionAtomD, deltaR[DB], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomD, positionAtomA, deltaR[AD], boxVectors);
        AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomD, positionAtomC, deltaR[CD], boxVectors);
    }
    else {
        AmoebaReferenceForce::loadDeltaR(positionAtomB, positionAtomA, deltaR[AB]);
        AmoebaReferenceForce::loadDeltaR(positionAtomB, positionAtomC, deltaR[CB]);
        AmoebaReferenceForce::loadDeltaR(positionAtomB, positionAtomD, deltaR[DB]);
        AmoebaReferenceForce::loadDeltaR(positionAtomD, positionAtomA, deltaR[AD]);
        AmoebaReferenceForce::loadDeltaR(positionAtomD, positionAtomC, deltaR[CD]);
    }

    RealOpenMM rDB2  = AmoebaReferenceForce::getNormSquared3(deltaR[DB]);
    RealOpenMM rAD2  = AmoebaReferenceForce::getNormSquared3(deltaR[AD]);
    RealOpenMM rCD2  = AmoebaReferenceForce::getNormSquared3(deltaR[CD]);
 
    std::vector<RealOpenMM> tempVector(3);
    AmoebaReferenceForce::getCrossProduct(deltaR[CB], deltaR[DB], tempVector);
    RealOpenMM   eE  = AmoebaReferenceForce::getDotProduct3(deltaR[AB], tempVector);
    RealOpenMM  dot  = AmoebaReferenceForce::getDotProduct3(deltaR[AD],  deltaR[CD]);
    RealOpenMM   cc  = rAD2*rCD2 - dot*dot;
 
    if (rDB2 <= zero || cc == zero) {
       return zero;
    }
    RealOpenMM bkk2   = rDB2 - eE*eE/cc;
    RealOpenMM cosine = SQRT(bkk2/rDB2);
    RealOpenMM angle;
    if (cosine >= one) {
       angle = zero;
    } else if (cosine <= -one) {
       angle = PI_M;
    } else {
       angle = RADIAN*ACOS(cosine);
    }
 
    // chain rule
 
    RealOpenMM dt    = angle;
    RealOpenMM dt2   = dt*dt;
    RealOpenMM dt3   = dt2*dt;
    RealOpenMM dt4   = dt2*dt2;
 
    RealOpenMM dEdDt = two + three*angleCubic*dt + four*angleQuartic*dt2 +
                       five*anglePentic*dt3 + six*angleSextic*dt4;
 
          dEdDt     *= angleK*dt*RADIAN;
 
    RealOpenMM dEdCos  = dEdDt/SQRT(cc*bkk2);
    if (eE > zero) {
       dEdCos *= -one;
    }
 
    RealOpenMM term = eE/cc;
 
    std::vector<RealOpenMM> dccd[LastAtomIndex];
    std::vector<RealOpenMM> deed[LastAtomIndex];
    std::vector<RealOpenMM> subForce[LastAtomIndex];
    for (int ii = 0; ii < LastAtomIndex; ii++) {
        dccd[ii].resize(3);
        deed[ii].resize(3);
        subForce[ii].resize(3);
    }   
    for (int ii = 0; ii < 3; ii++) {
       dccd[A][ii] = (deltaR[AD][ii]*rCD2 - deltaR[CD][ii]*dot)*term;
       dccd[C][ii] = (deltaR[CD][ii]*rAD2 - deltaR[AD][ii]*dot)*term;
       dccd[D][ii] = -one*(dccd[A][ii] + dccd[C][ii]);
    }
 
    AmoebaReferenceForce::getCrossProduct(deltaR[DB], deltaR[CB], deed[A]);
    AmoebaReferenceForce::getCrossProduct(deltaR[AB], deltaR[DB], deed[C]);
    AmoebaReferenceForce::getCrossProduct(deltaR[CB], deltaR[AB], deed[D]);
 
    term        = eE/rDB2;
    deed[D][0] += deltaR[DB][0]*term;
    deed[D][1] += deltaR[DB][1]*term;
    deed[D][2] += deltaR[DB][2]*term;
 
    // ---------------------------------------------------------------------------------------
 
    // forces
 
    // calculate forces for atoms a, c, d
    // the force for b is then -(a+ c + d)
 
 
    for (int jj = 0; jj < LastAtomIndex; jj++) {
 
       // A, C, D
 
       for (int ii = 0; ii < 3; ii++) {
          subForce[jj][ii] = dEdCos*(dccd[jj][ii] + deed[jj][ii]);
       }
 
       if (jj == 0)jj++; // skip B
 
       // now compute B
 
       if (jj == 3) {
          for (int ii = 0; ii < 3; ii++) {
             subForce[1][ii] = -one*(subForce[0][ii] + subForce[2][ii] + subForce[3][ii]);
          }
       }
    }
 
    // add in forces
 
    for (int jj = 0; jj < LastAtomIndex; jj++) {
       for (int ii = 0; ii < 3; ii++) {
          forces[jj][ii] = subForce[jj][ii];
       }
    }
 
    // ---------------------------------------------------------------------------------------
 
    // calculate energy if 'energy' is set
 
    RealOpenMM energy   = one + angleCubic*dt + angleQuartic*dt2 + anglePentic*dt3 + angleSextic*dt4;
    energy             *= angleK*dt2;

    return energy;
}

RealOpenMM AmoebaReferenceOutOfPlaneBendForce::calculateForceAndEnergy(int numOutOfPlaneBends, vector<RealVec>& posData,
                                                                       const std::vector<int>&  particle1,
                                                                       const std::vector<int>&  particle2,
                                                                       const std::vector<int>&  particle3,
                                                                       const std::vector<int>&  particle4,
                                                                       const std::vector<RealOpenMM>&  kQuadratic,
                                                                       RealOpenMM angleCubic,
                                                                       RealOpenMM angleQuartic,
                                                                       RealOpenMM anglePentic,
                                                                       RealOpenMM angleSextic,
                                                                       vector<RealVec>& forceData) const {
    RealOpenMM energy      = 0.0; 
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numOutOfPlaneBends); ii++) {
        int particle1Index      = particle1[ii];
        int particle2Index      = particle2[ii];
        int particle3Index      = particle3[ii];
        int particle4Index      = particle4[ii];
        RealOpenMM kAngle       = kQuadratic[ii];
        RealVec forces[4];
        energy                 += calculateOutOfPlaneBendIxn(posData[particle1Index], posData[particle2Index], posData[particle3Index], posData[particle4Index],
                                                              kAngle, angleCubic, angleQuartic, anglePentic, angleSextic, forces);
        for (int jj = 0; jj < 3; jj++) {
            forceData[particle1Index][jj] -= forces[0][jj];
            forceData[particle2Index][jj] -= forces[1][jj];
            forceData[particle3Index][jj] -= forces[2][jj];
            forceData[particle4Index][jj] -= forces[3][jj];
        }

    }   
    return energy;
}

