
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

ReferenceDynamics::ReferenceDynamics(int numberOfAtoms,  RealOpenMM deltaT, RealOpenMM temperature) : 
                  _numberOfAtoms(numberOfAtoms), _deltaT(deltaT), _temperature(temperature) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceDynamics::ReferenceDynamics";

   static const RealOpenMM one        =  1.0;

   // ---------------------------------------------------------------------------------------

   _timeStep             = 0;

   _ownReferenceConstraint = false;
   _referenceConstraint    = NULL;
}

/**---------------------------------------------------------------------------------------

   ReferenceDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceDynamics::~ReferenceDynamics() {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceDynamics::~ReferenceDynamics";

   // ---------------------------------------------------------------------------------------

   if (_ownReferenceConstraint) {
      delete _referenceConstraint;
   }
}

/**---------------------------------------------------------------------------------------

   Get number of atoms

   @return number of atoms

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::getNumberOfAtoms() const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getNumberOfAtoms";

   // ---------------------------------------------------------------------------------------

   return _numberOfAtoms;
}

/**---------------------------------------------------------------------------------------

   Get time step

   @return time step

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::getTimeStep() const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getTimeStep";

   // ---------------------------------------------------------------------------------------

   return _timeStep;
}

/**---------------------------------------------------------------------------------------

   Increment time step

   @return incremented time step

   --------------------------------------------------------------------------------------- */

int ReferenceDynamics::incrementTimeStep() {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getTimeStep";

   // ---------------------------------------------------------------------------------------

   return (++_timeStep);
}

/**---------------------------------------------------------------------------------------

   Get delta t

   @return deltaT

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceDynamics::getDeltaT() const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getDeltaT";

   // ---------------------------------------------------------------------------------------

   return _deltaT;
}

/**---------------------------------------------------------------------------------------

   Set delta t

   --------------------------------------------------------------------------------------- */

void ReferenceDynamics::setDeltaT(RealOpenMM deltaT) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::setDeltaT";

   // ---------------------------------------------------------------------------------------

   _deltaT = deltaT;
}

/**---------------------------------------------------------------------------------------

   Get temperature

   @return temperature

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceDynamics::getTemperature() const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getTemperature";

   // ---------------------------------------------------------------------------------------

   return _temperature;
}

/**---------------------------------------------------------------------------------------

   Get ReferenceConstraint

   @return ReferenceConstraint  object

   --------------------------------------------------------------------------------------- */

ReferenceConstraintAlgorithm* ReferenceDynamics::getReferenceConstraintAlgorithm() const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::getReferenceConstraint";

   // ---------------------------------------------------------------------------------------

   return _referenceConstraint;
}

/**---------------------------------------------------------------------------------------

   Set ReferenceConstraint

   @param referenceConstraint  ReferenceConstraint object

   --------------------------------------------------------------------------------------- */

void ReferenceDynamics::setReferenceConstraintAlgorithm(ReferenceConstraintAlgorithm* referenceConstraint) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceDynamics::setReferenceConstraint";

   // ---------------------------------------------------------------------------------------

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

void ReferenceDynamics::update(const OpenMM::System& system, vector<RealVec>& atomCoordinates,
                               vector<RealVec>& velocities, vector<RealVec>& forces, vector<RealOpenMM>& masses, RealOpenMM tolerance) {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceDynamics::update";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   // ---------------------------------------------------------------------------------------
}
