
/* Portions copyright (c) 2006-2013 Stanford University and Simbios.
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
#include <algorithm>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceVariableStochasticDynamics.h"
#include "ReferenceVirtualSites.h"
#include "openmm/OpenMMException.h"

#include <cstdio>

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceVariableStochasticDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param tau            viscosity(?)
   @param temperature    temperature
   @param accuracy       required accuracy

   --------------------------------------------------------------------------------------- */

ReferenceVariableStochasticDynamics::ReferenceVariableStochasticDynamics(int numberOfAtoms,
                                                          RealOpenMM tau, RealOpenMM temperature,
                                                          RealOpenMM accuracy) :
           ReferenceDynamics(numberOfAtoms, 0.0f, temperature), _tau(tau), _accuracy(accuracy) {
   if (tau <= 0) {
      std::stringstream message;
      message << "illegal tau value: " << tau;
      throw OpenMMException(message.str());
   }
   xPrime.resize(numberOfAtoms);
   inverseMasses.resize(numberOfAtoms);
}

/**---------------------------------------------------------------------------------------

   ReferenceVariableStochasticDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceVariableStochasticDynamics::~ReferenceVariableStochasticDynamics() {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceVariableStochasticDynamics::~ReferenceVariableStochasticDynamics";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Get the required accuracy

   @return accuracy

 --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceVariableStochasticDynamics::getAccuracy() const {
    return _accuracy;
}

/**---------------------------------------------------------------------------------------

   Set the required accuracy

 --------------------------------------------------------------------------------------- */

void ReferenceVariableStochasticDynamics::setAccuracy(RealOpenMM accuracy) {
    _accuracy = accuracy;
}

/**---------------------------------------------------------------------------------------

   Get tau

   @return tau

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceVariableStochasticDynamics::getTau() const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceVariableStochasticDynamics::getTau";

   // ---------------------------------------------------------------------------------------

   return _tau;
}

/**---------------------------------------------------------------------------------------

   First SD update; based on code in update.c do_update_sd() Gromacs 3.1.4

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses
   @param inverseMasses       inverse atom masses
   @param xPrime              xPrime
   @param maxStepSize         maximum time step

   --------------------------------------------------------------------------------------- */

void ReferenceVariableStochasticDynamics::updatePart1(int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                              vector<RealVec>& velocities,
                                              vector<RealVec>& forces, vector<RealOpenMM>& masses, vector<RealOpenMM>& inverseMasses,
                                              vector<RealVec>& xPrime, RealOpenMM maxStepSize) {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceVariableStochasticDynamics::updatePart1";

   // ---------------------------------------------------------------------------------------


   // first-time-through initialization

   if (getTimeStep() == 0) {
      // invert masses

      for (int ii = 0; ii < numberOfAtoms; ii++) {
         if (masses[ii] == 0)
             inverseMasses[ii] = 0;
         else
             inverseMasses[ii] = 1/masses[ii];
      }
   }

   // Select the step size to use
    RealOpenMM error = 0;
    for (int i = 0; i < numberOfAtoms; ++i) {
        for (int j = 0; j < 3; ++j) {
            RealOpenMM xerror = inverseMasses[i]*forces[i][j];
            error += xerror*xerror;
        }
    }
    error = SQRT(error/(numberOfAtoms*3));
    RealOpenMM newStepSize = SQRT(getAccuracy()/error);
    if (getDeltaT() > 0.0f)
        newStepSize = std::min(newStepSize, getDeltaT()*2.0f); // For safety, limit how quickly dt can increase.
    if (newStepSize > getDeltaT() && newStepSize < 1.2f*getDeltaT())
        newStepSize = getDeltaT(); // Keeping dt constant between steps improves the behavior of the integrator.
    if (newStepSize > maxStepSize)
        newStepSize = maxStepSize;
    setDeltaT(newStepSize);
 
    // perform first update

   RealOpenMM tau = getTau();
   const RealOpenMM vscale = EXP(-getDeltaT()/tau);
   const RealOpenMM fscale = (1-vscale)*tau;
   const RealOpenMM kT = BOLTZ*getTemperature();
   const RealOpenMM noisescale = SQRT(2*kT/tau)*SQRT(0.5*(1-vscale*vscale)*tau);

   for (int ii = 0; ii < numberOfAtoms; ii++) {
       if (masses[ii] != 0) {
           RealOpenMM sqrtInvMass = SQRT(inverseMasses[ii]);
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

void ReferenceVariableStochasticDynamics::updatePart2(int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                              vector<RealVec>& velocities,
                                              vector<RealVec>& forces, vector<RealOpenMM>& inverseMasses,
                                              vector<RealVec>& xPrime) {

   // ---------------------------------------------------------------------------------------

   //static const char* methodName  = "\nReferenceVariableStochasticDynamics::updatePart2";

   // ---------------------------------------------------------------------------------------

   // perform second update

   for (int ii = 0; ii < numberOfAtoms; ii++) {
       if (inverseMasses[ii] != 0.0)
           for (int jj = 0; jj < 3; jj++)
               xPrime[ii][jj] = atomCoordinates[ii][jj]+getDeltaT()*velocities[ii][jj];
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

void ReferenceVariableStochasticDynamics::update(const OpenMM::System& system, vector<RealVec>& atomCoordinates,
                                          vector<RealVec>& velocities,
                                          vector<RealVec>& forces, vector<RealOpenMM>& masses, RealOpenMM maxStepSize, RealOpenMM tolerance) {

   // ---------------------------------------------------------------------------------------

   //static const char* methodName      = "\nReferenceVariableStochasticDynamics::update";

   // ---------------------------------------------------------------------------------------

   // 1st update

   int numberOfAtoms = system.getNumParticles();
   updatePart1(numberOfAtoms, atomCoordinates, velocities, forces, masses, inverseMasses, xPrime, maxStepSize);

   // 2nd update

   updatePart2(numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses, xPrime);

   ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
   if (referenceConstraintAlgorithm)
      referenceConstraintAlgorithm->apply(atomCoordinates, xPrime, inverseMasses, tolerance);

   // copy xPrime -> atomCoordinates

   RealOpenMM invStepSize = 1.0/getDeltaT();
   for (int ii = 0; ii < numberOfAtoms; ii++) {
       if (masses[ii] != 0.0) {
           velocities[ii] = (xPrime[ii]-atomCoordinates[ii])*invStepSize;
           atomCoordinates[ii] = xPrime[ii];
       }
   }

   ReferenceVirtualSites::computePositions(system, atomCoordinates);
   incrementTimeStep();
}
