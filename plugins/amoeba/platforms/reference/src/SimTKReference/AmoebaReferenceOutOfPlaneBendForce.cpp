
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
#include "SimTKOpenMMRealType.h"

using std::vector;
using namespace OpenMM;

void AmoebaReferenceOutOfPlaneBendForce::setPeriodic(OpenMM::Vec3* vectors) {
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

double AmoebaReferenceOutOfPlaneBendForce::calculateOutOfPlaneBendIxn(const Vec3& positionAtomA, const Vec3& positionAtomB,
                                                                      const Vec3& positionAtomC, const Vec3& positionAtomD,
                                                                      double angleK,
                                                                      double angleCubic,                 double angleQuartic,
                                                                      double anglePentic,                double angleSextic,
                                                                      Vec3* forces) const {

    enum { A, B, C, D, LastAtomIndex };
    enum { AB, CB, DB, AD, CD, LastDeltaIndex };
 
    // get deltaR between various combinations of the 4 atoms
    // and various intermediate terms
 
    std::vector<double> deltaR[LastDeltaIndex];
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

    double rDB2  = AmoebaReferenceForce::getNormSquared3(deltaR[DB]);
    double rAD2  = AmoebaReferenceForce::getNormSquared3(deltaR[AD]);
    double rCD2  = AmoebaReferenceForce::getNormSquared3(deltaR[CD]);
 
    std::vector<double> tempVector(3);
    AmoebaReferenceForce::getCrossProduct(deltaR[CB], deltaR[DB], tempVector);
    double eE = AmoebaReferenceForce::getDotProduct3(deltaR[AB], tempVector);
    double dot = AmoebaReferenceForce::getDotProduct3(deltaR[AD],  deltaR[CD]);
    double cc = rAD2*rCD2 - dot*dot;
 
    if (rDB2 <= 0.0 || cc == 0.0) {
       return 0.0;
    }
    double bkk2   = rDB2 - eE*eE/cc;
    double cosine = sqrt(bkk2/rDB2);
    double angle;
    if (cosine >= 1.0) {
       angle = 0.0;
    } else if (cosine <= -1.0) {
       angle = M_PI;
    } else {
       angle = RADIAN*acos(cosine);
    }
 
    // chain rule
 
    double dt    = angle;
    double dt2   = dt*dt;
    double dt3   = dt2*dt;
    double dt4   = dt2*dt2;
 
    double dEdDt = 2.0 + 3.0*angleCubic*dt + 4.0*angleQuartic*dt2 +
                   5.0*anglePentic*dt3 + 6.0*angleSextic*dt4;
 
          dEdDt     *= angleK*dt*RADIAN;
 
    double dEdCos  = dEdDt/sqrt(cc*bkk2);
    if (eE > 0.0) {
       dEdCos *= -1.0;
    }
 
    double term = eE/cc;
 
    std::vector<double> dccd[LastAtomIndex];
    std::vector<double> deed[LastAtomIndex];
    std::vector<double> subForce[LastAtomIndex];
    for (int ii = 0; ii < LastAtomIndex; ii++) {
        dccd[ii].resize(3);
        deed[ii].resize(3);
        subForce[ii].resize(3);
    }   
    for (int ii = 0; ii < 3; ii++) {
       dccd[A][ii] = (deltaR[AD][ii]*rCD2 - deltaR[CD][ii]*dot)*term;
       dccd[C][ii] = (deltaR[CD][ii]*rAD2 - deltaR[AD][ii]*dot)*term;
       dccd[D][ii] = -1.0*(dccd[A][ii] + dccd[C][ii]);
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
             subForce[1][ii] = -1.0*(subForce[0][ii] + subForce[2][ii] + subForce[3][ii]);
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
 
    double energy = 1.0 + angleCubic*dt + angleQuartic*dt2 + anglePentic*dt3 + angleSextic*dt4;
    energy *= angleK*dt2;

    return energy;
}

double AmoebaReferenceOutOfPlaneBendForce::calculateForceAndEnergy(int numOutOfPlaneBends, vector<Vec3>& posData,
                                                                   const std::vector<int>&  particle1,
                                                                   const std::vector<int>&  particle2,
                                                                   const std::vector<int>&  particle3,
                                                                   const std::vector<int>&  particle4,
                                                                   const std::vector<double>&  kQuadratic,
                                                                   double angleCubic,
                                                                   double angleQuartic,
                                                                   double anglePentic,
                                                                   double angleSextic,
                                                                   vector<Vec3>& forceData) const {
    double energy      = 0.0; 
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numOutOfPlaneBends); ii++) {
        int particle1Index      = particle1[ii];
        int particle2Index      = particle2[ii];
        int particle3Index      = particle3[ii];
        int particle4Index      = particle4[ii];
        double kAngle           = kQuadratic[ii];
        Vec3 forces[4];
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

