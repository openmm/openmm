
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

#include <string.h>
#include <sstream>

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceConstraint.h"
#include "ReferenceDynamics.h"

/**---------------------------------------------------------------------------------------

   ReferenceConstraint constructor

   --------------------------------------------------------------------------------------- */

ReferenceConstraint::ReferenceConstraint( ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName = "\nReferenceConstraint::ReferenceConstraint";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceConstraint destructor

   --------------------------------------------------------------------------------------- */

ReferenceConstraint::~ReferenceConstraint( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceConstraint::~ReferenceConstraint";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   ReferenceShakeConstraint constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param tau            viscosity
   @param temperature    temperature

   --------------------------------------------------------------------------------------- */

ReferenceShakeConstraint::ReferenceShakeConstraint( int atomIndex1, int atomIndex2, 
                                                    RealOpenMM constraintDistance, 
                                                    RealOpenMM inverseMass1,
                                                    RealOpenMM inverseMass2 ) : ReferenceConstraint( ) {

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nReferenceShakeConstraint::ReferenceShakeConstraint";

   // ---------------------------------------------------------------------------------------

   // insure heavy atom is in 0 slot

   if( inverseMass1 > inverseMass2 ){
      int temp         = atomIndex1;
      atomIndex1       = atomIndex2;
      atomIndex2       = temp;

      RealOpenMM tempR = inverseMass1;
      inverseMass1     = inverseMass2;
      inverseMass2     = tempR;
   }

   _atomIndices[0]     = atomIndex1;
   _atomIndices[1]     = atomIndex2;

   _inverseMasses[0]   = inverseMass1;
   _inverseMasses[1]   = inverseMass2;

   _constraintDistance = constraintDistance;

}

/**---------------------------------------------------------------------------------------

   ReferenceShakeConstraint destructor

   --------------------------------------------------------------------------------------- */

ReferenceShakeConstraint::~ReferenceShakeConstraint( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceShakeConstraint::~ReferenceShakeConstraint";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Get constraint distance

   @return constraintDistance

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceShakeConstraint::getConstraintDistance( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceShakeConstraint::getConstraintDistance";

   // ---------------------------------------------------------------------------------------

   return _constraintDistance;
}

/**---------------------------------------------------------------------------------------

   Get inverse mass of heavy atom

   @return inverse mass

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceShakeConstraint::getHeavyAtomInverseMass( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceShakeConstraint::getHeavyAtomInverseMass";

   // ---------------------------------------------------------------------------------------

   return _inverseMasses[0];
}

/**---------------------------------------------------------------------------------------

   Get inverse mass of light atom

   @return inverse mass

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceShakeConstraint::getLightAtomInverseMass( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceShakeConstraint::getLightAtomInverseMass";

   // ---------------------------------------------------------------------------------------

   return _inverseMasses[1];
}

/**---------------------------------------------------------------------------------------

   Get index of heavy atom

   @return index

   --------------------------------------------------------------------------------------- */

int ReferenceShakeConstraint::getHeavyAtomIndex( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceShakeConstraint::getHeavyAtomIndex";

   // ---------------------------------------------------------------------------------------

   return _atomIndices[0];
}

/**---------------------------------------------------------------------------------------

   Get index of light atom

   @return index

   --------------------------------------------------------------------------------------- */

int ReferenceShakeConstraint::getLightAtomIndex( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceShakeConstraint::getLightAtomIndex";

   // ---------------------------------------------------------------------------------------

   return _atomIndices[1];
}
