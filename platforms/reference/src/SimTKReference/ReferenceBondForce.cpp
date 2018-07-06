
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
#include "ReferenceBondForce.h"

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceBondForce constructor

   --------------------------------------------------------------------------------------- */

ReferenceBondForce::ReferenceBondForce() {
}

/**---------------------------------------------------------------------------------------

   ReferenceBondForce destructor

   --------------------------------------------------------------------------------------- */

ReferenceBondForce::~ReferenceBondForce() {
}

/**---------------------------------------------------------------------------------------

   Calculate forces/energy for bonds

   @param numberOfBonds    number of bonds
   @param atomIndices      indices of atoms participating in bond ixn: atomIndices[bondIndex][indices]
   @param atomCoordinates  atom coordinates: atomCoordinates[atomIndex][3]
   @param parameters       parameters: parameters[bondIndex][*]; contents of array 
                           depend on ixn
   @param forces           force array (forces added to current values): forces[atomIndex][3]
   @param totalEnergy      if not null, the energy will be added to this
   @param ReferenceBondIxn ixn to be calculated

   --------------------------------------------------------------------------------------- */

void ReferenceBondForce::calculateForce(int numberOfBonds, vector<vector<int> >& atomIndices,
                                        vector<Vec3>& atomCoordinates,
                                        vector<vector<double> >& parameters,
                                        vector<Vec3>& forces, 
                                        double *totalEnergy, 
                                        ReferenceBondIxn& referenceBondIxn) {

   for (int ii = 0; ii < numberOfBonds; ii++) {

      // calculate bond ixn

      referenceBondIxn.calculateBondIxn(atomIndices[ii], atomCoordinates, parameters[ii], 
                                        forces, totalEnergy, NULL);
   }
}

