
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

#include <cstring>
#include <sstream>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceDynamics.h"

#include <cstdio>

using std::vector;
using namespace OpenMM;


/**---------------------------------------------------------------------------------------

   ReferenceDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param temperature    temperature

   --------------------------------------------------------------------------------------- */

ReferenceDynamics::ReferenceDynamics(int numberOfAtoms,  double deltaT, double temperature) : 
                  _numberOfAtoms(numberOfAtoms), _deltaT(deltaT), _temperature(temperature) {

   _timeStep = 0;
   _ownReferenceConstraint = false;
   _referenceConstraint    = NULL;
}

/**---------------------------------------------------------------------------------------

   ReferenceDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceDynamics::~ReferenceDynamics() {
   if (_ownReferenceConstraint) {
      delete _referenceConstraint;
   }
}

/**---------------------------------------------------------------------------------------

   Get number of atoms

   @return number of atoms

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::getNumberOfAtoms() const {
   return _numberOfAtoms;
}

/**---------------------------------------------------------------------------------------

   Get time step

   @return time step

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::getTimeStep() const {
   return _timeStep;
}

/**---------------------------------------------------------------------------------------

   Increment time step

   @return incremented time step

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::incrementTimeStep() {
   return (++_timeStep);
}

/**---------------------------------------------------------------------------------------

   Get delta t

   @return deltaT

   --------------------------------------------------------------------------------------- */

double ReferenceDynamics::getDeltaT() const {
   return _deltaT;
}

/**---------------------------------------------------------------------------------------

   Set delta t

   --------------------------------------------------------------------------------------- */

void ReferenceDynamics::setDeltaT(double deltaT) {
   _deltaT = deltaT;
}

/**---------------------------------------------------------------------------------------

   Get temperature

   @return temperature

   --------------------------------------------------------------------------------------- */

double ReferenceDynamics::getTemperature() const {
   return _temperature;
}

/**---------------------------------------------------------------------------------------

   Get ReferenceConstraint

   @return ReferenceConstraint  object

   --------------------------------------------------------------------------------------- */

ReferenceConstraintAlgorithm* ReferenceDynamics::getReferenceConstraintAlgorithm() const {
   return _referenceConstraint;
}

/**---------------------------------------------------------------------------------------

   Set ReferenceConstraint

   @param referenceConstraint  ReferenceConstraint object

   --------------------------------------------------------------------------------------- */

void ReferenceDynamics::setReferenceConstraintAlgorithm(ReferenceConstraintAlgorithm* referenceConstraint) {
   // delete if own

   if (_referenceConstraint && _ownReferenceConstraint) {
      delete _referenceConstraint;
   }

   _referenceConstraint = referenceConstraint;
   _ownReferenceConstraint = 0;
}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing dynamics update of coordinates
   and velocities

   @param system              the System to be integrated
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses
   @param tolerance           the constraint tolerance

   --------------------------------------------------------------------------------------- */

void ReferenceDynamics::update(const OpenMM::System& system, vector<Vec3>& atomCoordinates,
                               vector<Vec3>& velocities, vector<Vec3>& forces, vector<double>& masses, double tolerance) {
}
