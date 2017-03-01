
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

#include "AmoebaReferenceBondForce.h"
#include "AmoebaReferenceForce.h"

using std::vector;
using namespace OpenMM;

void AmoebaReferenceBondForce::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

/**---------------------------------------------------------------------------------------

   Calculate Amoeba bond ixn (force and energy)

   @param positionAtomA           Cartesian coordinates of atom A
   @param positionAtomB           Cartesian coordinates of atom B
   @param bondLength              bond length
   @param bondK                   bond force
   @param bondCubic               cubic bond force parameter
   @param bondQuartic             quartic bond force parameter
   @param forces                  force vector

   @return energy

   --------------------------------------------------------------------------------------- */

double AmoebaReferenceBondForce::calculateBondIxn(const Vec3& positionAtomA, const Vec3& positionAtomB,
                                                  double bondLength, double bondK,
                                                  double bondCubic, double bondQuartic,
                                                  Vec3* forces) const {
   // get deltaR, R2, and R between 2 atoms

   std::vector<double> deltaR;
   if (usePeriodic)
       AmoebaReferenceForce::loadDeltaRPeriodic(positionAtomA, positionAtomB, deltaR, boxVectors);
   else
       AmoebaReferenceForce::loadDeltaR(positionAtomA, positionAtomB, deltaR);
   double r = AmoebaReferenceForce::getNorm3(deltaR);

   // deltaIdeal = r - r_0

   double deltaIdeal      = r - bondLength;
   double deltaIdeal2     = deltaIdeal*deltaIdeal;

   double dEdR            = (1.0 + 1.5*bondCubic*deltaIdeal + 2.0*bondQuartic*deltaIdeal2);
   dEdR                  *= 2.0*bondK*deltaIdeal;
   dEdR                   = r > 0.0 ? (dEdR/r) : 0.0;

   forces[0][0]           = dEdR*deltaR[0];
   forces[0][1]           = dEdR*deltaR[1];
   forces[0][2]           = dEdR*deltaR[2];

   dEdR                  *= -1.0;
   forces[1][0]           = dEdR*deltaR[0];
   forces[1][1]           = dEdR*deltaR[1];
   forces[1][2]           = dEdR*deltaR[2];

   double energy          = bondK*deltaIdeal2*(1.0 + bondCubic*deltaIdeal + bondQuartic*deltaIdeal2);
   return energy;
}

double AmoebaReferenceBondForce::calculateForceAndEnergy(int numBonds,
                                                         vector<Vec3>& particlePositions,
                                                         const std::vector<int>&   particle1,
                                                         const std::vector<int>&   particle2,
                                                         const std::vector<double>& length,
                                                         const std::vector<double>& kQuadratic,
                                                         double globalBondCubic,
                                                         double globalBondQuartic,
                                                         vector<Vec3>& forceData) const {
    double energy = 0.0; 
    for (int ii = 0; ii < numBonds; ii++) {
        int particle1Index = particle1[ii];
        int particle2Index = particle2[ii];
        double bondLength = length[ii];
        double bondK = kQuadratic[ii];
        Vec3 forces[2];

        energy += calculateBondIxn(particlePositions[particle1Index], particlePositions[particle2Index],
                                   bondLength, bondK, globalBondCubic, globalBondQuartic,
                                   forces);

        for (int jj = 0; jj < 3; jj++) {
            forceData[particle1Index][jj] += forces[0][jj];
            forceData[particle2Index][jj] += forces[1][jj];
        }

    }   
    return energy;
}
