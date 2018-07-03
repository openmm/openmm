
/* Portions copyright (c) 2006-2016 Stanford University and Simbios.
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

#include "SimTKOpenMMUtilities.h"
#include "ReferenceLJCoulomb14.h"
#include "ReferenceForce.h"

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceLJCoulomb14 constructor

   --------------------------------------------------------------------------------------- */

ReferenceLJCoulomb14::ReferenceLJCoulomb14() {
}

/**---------------------------------------------------------------------------------------

   ReferenceLJCoulomb14 destructor

   --------------------------------------------------------------------------------------- */

ReferenceLJCoulomb14::~ReferenceLJCoulomb14() {
}

/**---------------------------------------------------------------------------------------

   Calculate LJ 1-4 ixn

   @param atomIndices      atom indices of 4 atoms in bond
   @param atomCoordinates  atom coordinates
   @param parameters       three parameters:
                                        parameters[0]= (c12/c6)**1/6  (sigma)
                                        parameters[1]= c6*c6/c12      (4*epsilon)
                                        parameters[2]= epsfac*q1*q2
   @param forces           force array (forces added to current values)
   @param totalEnergy      if not null, the energy will be added to this

   --------------------------------------------------------------------------------------- */

void ReferenceLJCoulomb14::calculateBondIxn(vector<int>& atomIndices, vector<Vec3>& atomCoordinates,
                                     vector<double>& parameters, vector<Vec3>& forces,
                                     double* totalEnergy, double* energyParamDerivs) {
   double deltaR[2][ReferenceForce::LastDeltaRIndex];

   // get deltaR, R2, and R between 2 atoms

   int atomAIndex = atomIndices[0];
   int atomBIndex = atomIndices[1];
   ReferenceForce::getDeltaR(atomCoordinates[atomBIndex], atomCoordinates[atomAIndex], deltaR[0]);  

   double inverseR  = 1.0/(deltaR[0][ReferenceForce::RIndex]);
   double sig2      = inverseR*parameters[0];
          sig2     *= sig2;
   double sig6      = sig2*sig2*sig2;

   double dEdR      = parameters[1]*(12.0*sig6 - 6.0)*sig6;
          dEdR     += ONE_4PI_EPS0*parameters[2]*inverseR;
          dEdR     *= inverseR*inverseR;

   // accumulate forces

   for (int ii = 0; ii < 3; ii++) {
      double force        = dEdR*deltaR[0][ii];
      forces[atomAIndex][ii] += force;
      forces[atomBIndex][ii] -= force;
   }

   // accumulate energies

   if (totalEnergy != NULL)
       *totalEnergy += parameters[1]*(sig6 - 1.0)*sig6 + (ONE_4PI_EPS0*parameters[2]*inverseR);
}
