
/* Portions copyright (c) 2010-2016 Stanford University and Simbios.
 * Contributors: Peter Eastman
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

#include "ReferenceCMAPTorsionIxn.h"
#include "ReferenceForce.h"

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   Constructor

   --------------------------------------------------------------------------------------- */

ReferenceCMAPTorsionIxn::ReferenceCMAPTorsionIxn(const vector<vector<vector<double> > >& coeff,
        const vector<int>& torsionMaps, const vector<vector<int> >& torsionIndices) :
        coeff(coeff), torsionMaps(torsionMaps), torsionIndices(torsionIndices), usePeriodic(false) {
}

void ReferenceCMAPTorsionIxn::setPeriodic(OpenMM::Vec3* vectors) {
    usePeriodic = true;
    boxVectors[0] = vectors[0];
    boxVectors[1] = vectors[1];
    boxVectors[2] = vectors[2];
}

/**---------------------------------------------------------------------------------------

   Calculate torsion interaction

   @param atomIndices      bond indices
   @param atomCoordinates  atom coordinates
   @param parameters       parameters
   @param forces           force array (forces added)
   @param energyByBond     bond energy
   @param energy           atom energy

   --------------------------------------------------------------------------------------- */

void ReferenceCMAPTorsionIxn::calculateIxn(vector<Vec3>& atomCoordinates, vector<Vec3>& forces, double* totalEnergy) const {
    for (unsigned int i = 0; i < torsionMaps.size(); i++)
        calculateOneIxn(i, atomCoordinates, forces, totalEnergy);
}

/**---------------------------------------------------------------------------------------

   Calculate the interaction due to a single torsion pair

   @param index            the index of the torsion
   @param atomCoordinates  atom coordinates
   @param forces           force array (forces added)
   @param totalEnergy      total energy

     --------------------------------------------------------------------------------------- */

void ReferenceCMAPTorsionIxn::calculateOneIxn(int index, vector<Vec3>& atomCoordinates, vector<Vec3>& forces,
                     double* totalEnergy) const {
    int map = torsionMaps[index];
    int a1 = torsionIndices[index][0];
    int a2 = torsionIndices[index][1];
    int a3 = torsionIndices[index][2];
    int a4 = torsionIndices[index][3];
    int b1 = torsionIndices[index][4];
    int b2 = torsionIndices[index][5];
    int b3 = torsionIndices[index][6];
    int b4 = torsionIndices[index][7];

    // Compute deltas between the various atoms involved.

    double deltaA[3][ReferenceForce::LastDeltaRIndex];
    double deltaB[3][ReferenceForce::LastDeltaRIndex];
    if (usePeriodic) {
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[a2], atomCoordinates[a1], boxVectors, deltaA[0]);
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[a2], atomCoordinates[a3], boxVectors, deltaA[1]);
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[a4], atomCoordinates[a3], boxVectors, deltaA[2]);
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[b2], atomCoordinates[b1], boxVectors, deltaB[0]);
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[b2], atomCoordinates[b3], boxVectors, deltaB[1]);
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[b4], atomCoordinates[b3], boxVectors, deltaB[2]);
    }
    else {
        ReferenceForce::getDeltaR(atomCoordinates[a2], atomCoordinates[a1], deltaA[0]);
        ReferenceForce::getDeltaR(atomCoordinates[a2], atomCoordinates[a3], deltaA[1]);
        ReferenceForce::getDeltaR(atomCoordinates[a4], atomCoordinates[a3], deltaA[2]);
        ReferenceForce::getDeltaR(atomCoordinates[b2], atomCoordinates[b1], deltaB[0]);
        ReferenceForce::getDeltaR(atomCoordinates[b2], atomCoordinates[b3], deltaB[1]);
        ReferenceForce::getDeltaR(atomCoordinates[b4], atomCoordinates[b3], deltaB[2]);
    }

    // Visual Studio complains if crossProduct declared as 'crossProduct[2][3]'

    double crossProductMemory[12];
    double* cpA[2];
    cpA[0] = crossProductMemory;
    cpA[1] = crossProductMemory + 3;
    double* cpB[2];
    cpB[0] = crossProductMemory + 6;
    cpB[1] = crossProductMemory + 9;

   // Compute the dihedral angles.

    double dotDihedral;
    double signOfAngle;
    double angleA =  getDihedralAngleBetweenThreeVectors(deltaA[0], deltaA[1], deltaA[2],
         cpA, &dotDihedral, deltaA[0], &signOfAngle, 1);
    double angleB =  getDihedralAngleBetweenThreeVectors(deltaB[0], deltaB[1], deltaB[2],
         cpB, &dotDihedral, deltaB[0], &signOfAngle, 1);
    angleA = fmod(angleA+2.0*M_PI, 2.0*M_PI);
    angleB = fmod(angleB+2.0*M_PI, 2.0*M_PI);

    // Identify which patch this is in.

    int size = (int) sqrt(coeff[map].size());
    double delta = 2*M_PI/size;
    int s = (int) fmin(angleA/delta, size-1);
    int t = (int) fmin(angleB/delta, size-1);
    const vector<double>& c = coeff[map][s+size*t];
    double da = angleA/delta-s;
    double db = angleB/delta-t;

    // Evaluate the spline to determine the energy and gradients.

    double energy = 0;
    double dEdA = 0;
    double dEdB = 0;
    for (int i = 3; i >= 0; i--) {
        energy = da*energy + ((c[i*4+3]*db + c[i*4+2])*db + c[i*4+1])*db + c[i*4+0];
        dEdA = db*dEdA + (3.0*c[i+3*4]*da + 2.0*c[i+2*4])*da + c[i+1*4];
        dEdB = da*dEdB + (3.0*c[i*4+3]*db + 2.0*c[i*4+2])*db + c[i*4+1];
    }
    dEdA /= delta;
    dEdB /= delta;
    if (totalEnergy != NULL)
        *totalEnergy += energy;

    // Apply the force to the first torsion.

    double forceFactors[4];
    double normCross1 = DOT3(cpA[0], cpA[0]);
    double normBC = deltaA[1][ReferenceForce::RIndex];
    forceFactors[0] = (-dEdA*normBC)/normCross1;
    double normCross2 = DOT3(cpA[1], cpA[1]);
    forceFactors[3] = (dEdA*normBC)/normCross2;
    forceFactors[1] = DOT3(deltaA[0], deltaA[1]);
    forceFactors[1] /= deltaA[1][ReferenceForce::R2Index];
    forceFactors[2] = DOT3(deltaA[2], deltaA[1]);
    forceFactors[2] /= deltaA[1][ReferenceForce::R2Index];
    for (int i = 0; i < 3; i++) {
        double f0 = forceFactors[0]*cpA[0][i];
        double f3 = forceFactors[3]*cpA[1][i];
        double s = forceFactors[1]*f0 - forceFactors[2]*f3;
        forces[a1][i] += f0;
        forces[a2][i] -= f0-s;
        forces[a3][i] -= f3+s;
        forces[a4][i] += f3;
    }

    // Apply the force to the second torsion.

    normCross1 = DOT3(cpB[0], cpB[0]);
    normBC = deltaB[1][ReferenceForce::RIndex];
    forceFactors[0] = (-dEdB*normBC)/normCross1;
    normCross2 = DOT3(cpB[1], cpB[1]);
    forceFactors[3] = (dEdB*normBC)/normCross2;
    forceFactors[1] = DOT3(deltaB[0], deltaB[1]);
    forceFactors[1] /= deltaB[1][ReferenceForce::R2Index];
    forceFactors[2] = DOT3(deltaB[2], deltaB[1]);
    forceFactors[2] /= deltaB[1][ReferenceForce::R2Index];
    for (int i = 0; i < 3; i++) {
        double f0 = forceFactors[0]*cpB[0][i];
        double f3 = forceFactors[3]*cpB[1][i];
        double s = forceFactors[1]*f0 - forceFactors[2]*f3;
        forces[b1][i] += f0;
        forces[b2][i] -= f0-s;
        forces[b3][i] -= f3+s;
        forces[b4][i] += f3;
    }
}

// ---------------------------------------------------------------------------------------

/**---------------------------------------------------------------------------------------

   This is present only because we must define it to subclass ReferenceBondIxn.  It is never called.

   --------------------------------------------------------------------------------------- */

void ReferenceCMAPTorsionIxn::calculateBondIxn(vector<int>& atomIndices, vector<Vec3>& atomCoordinates,
        vector<double>& parameters, vector<Vec3>& forces, double* totalEnergy, double* energyParamDerivs) {
}
