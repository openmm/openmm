
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

#include <cstring>
#include <sstream>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceStochasticDynamics.h"
#include "ReferenceVirtualSites.h"
#include "openmm/OpenMMException.h"

#include <cstdio>

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceStochasticDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param friction       friction coefficient
   @param temperature    temperature

   --------------------------------------------------------------------------------------- */

ReferenceStochasticDynamics::ReferenceStochasticDynamics(int numberOfAtoms,
                                                         double deltaT, double friction,
                                                         double temperature) : 
           ReferenceDynamics(numberOfAtoms, deltaT, temperature), friction(friction) {
   xPrime.resize(numberOfAtoms);
   inverseMasses.resize(numberOfAtoms);
}

/**---------------------------------------------------------------------------------------

   ReferenceStochasticDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceStochasticDynamics::~ReferenceStochasticDynamics() {
}

/**---------------------------------------------------------------------------------------

   Get friction coefficient

   --------------------------------------------------------------------------------------- */

double ReferenceStochasticDynamics::getFriction() const {
   return friction;
}

/**---------------------------------------------------------------------------------------

   First SD update; based on code in update.c do_update_sd() Gromacs 3.1.4

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param inverseMasses       inverse atom masses
   @param xPrime              xPrime

   --------------------------------------------------------------------------------------- */

void ReferenceStochasticDynamics::updatePart1(int numberOfAtoms, vector<Vec3>& atomCoordinates,
                                              vector<Vec3>& velocities,
                                              vector<Vec3>& forces, vector<double>& inverseMasses,
                                              vector<Vec3>& xPrime) {
   // perform first update

   double dt = getDeltaT();
   double friction = getFriction();
   const double vscale = exp(-dt*friction);
   const double fscale = (friction == 0 ? dt : (1-vscale)/friction);
   const double kT = BOLTZ*getTemperature();
   const double noisescale = sqrt(kT*(1-vscale*vscale));

   for (int ii = 0; ii < numberOfAtoms; ii++) {
       if (inverseMasses[ii] != 0.0) {
           double sqrtInvMass = sqrt(inverseMasses[ii]);
           for (int jj = 0; jj < 3; jj++) {
               velocities[ii][jj]  = vscale*velocities[ii][jj] + fscale*inverseMasses[ii]*forces[ii][jj] + noisescale*sqrtInvMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
           }
       }
   }
}

/**---------------------------------------------------------------------------------------

   Second update; based on code in update.c do_update_sd() w/ bFirstHalf = false in Gromacs 3.1.4

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   --------------------------------------------------------------------------------------- */

void ReferenceStochasticDynamics::updatePart2(int numberOfAtoms, vector<Vec3>& atomCoordinates,
                                              vector<Vec3>& velocities,
                                              vector<Vec3>& forces, vector<double>& inverseMasses,
                                              vector<Vec3>& xPrime) {
   // perform second update

   for (int ii = 0; ii < numberOfAtoms; ii++) {
       if (inverseMasses[ii] != 0.0)
            xPrime[ii] = atomCoordinates[ii]+velocities[ii]*getDeltaT();
   }
}

void ReferenceStochasticDynamics::updatePart3(int numberOfAtoms, vector<Vec3>& atomCoordinates,
                                              vector<Vec3>& velocities, vector<double>& inverseMasses,
                                              vector<Vec3>& xPrime) {
   double invStepSize = 1.0/getDeltaT();
   for (int i = 0; i < numberOfAtoms; ++i)
       if (inverseMasses[i] != 0) {
            velocities[i] = (xPrime[i]-atomCoordinates[i])*invStepSize;
            atomCoordinates[i] = xPrime[i];
       }
}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing stochastic dynamics update of coordinates
   and velocities

   @param system              the System to be integrated
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   --------------------------------------------------------------------------------------- */

void ReferenceStochasticDynamics::update(const OpenMM::System& system, vector<Vec3>& atomCoordinates,
                                          vector<Vec3>& velocities, vector<Vec3>& forces, vector<double>& masses, double tolerance) {
   // first-time-through initialization

   int numberOfAtoms = system.getNumParticles();
   if (getTimeStep() == 0) {
      // invert masses

      for (int ii = 0; ii < numberOfAtoms; ii++) {
         if (masses[ii] == 0.0)
             inverseMasses[ii] = 0.0;
         else
             inverseMasses[ii] = 1.0/masses[ii];
      }
   }

   // 1st update

   updatePart1(numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses, xPrime);

   // 2nd update

   updatePart2(numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses, xPrime);

   ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
   if (referenceConstraintAlgorithm)
      referenceConstraintAlgorithm->apply(atomCoordinates, xPrime, inverseMasses, tolerance);

   // copy xPrime -> atomCoordinates

   updatePart3(numberOfAtoms, atomCoordinates, velocities, inverseMasses, xPrime);

   ReferenceVirtualSites::computePositions(system, atomCoordinates);
   incrementTimeStep();
}
