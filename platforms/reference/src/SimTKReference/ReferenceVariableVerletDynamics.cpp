
/* Portions copyright (c) 2006-2009 Stanford University and Simbios.
 * Contributors: Peter Eastman, Pande Group
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
#include "ReferenceVariableVerletDynamics.h"

/**---------------------------------------------------------------------------------------

   ReferenceVariableVerletDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         initial delta t for dynamics
   @param accuracy       required accuracy

   --------------------------------------------------------------------------------------- */

ReferenceVariableVerletDynamics::ReferenceVariableVerletDynamics( int numberOfAtoms,
                                                          RealOpenMM deltaT, RealOpenMM accuracy ) :
           ReferenceDynamics( numberOfAtoms, deltaT, 0.0 ), _accuracy(accuracy) {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceVariableVerletDynamics::ReferenceVariableVerletDynamics";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   // ---------------------------------------------------------------------------------------

   allocate2DArrays( numberOfAtoms, 3, Max2DArrays );
   allocate1DArrays( numberOfAtoms, Max1DArrays );

}

/**---------------------------------------------------------------------------------------

   ReferenceVariableVerletDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceVariableVerletDynamics::~ReferenceVariableVerletDynamics( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceVariableVerletDynamics::~ReferenceVariableVerletDynamics";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Get the required accuracy

   @return accuracy

 --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceVariableVerletDynamics::getAccuracy( void ) const {
    return _accuracy;
}

/**---------------------------------------------------------------------------------------

   Set the required accuracy

 --------------------------------------------------------------------------------------- */

void ReferenceVariableVerletDynamics::setAccuracy( RealOpenMM accuracy ) {
    _accuracy = accuracy;
}

/**---------------------------------------------------------------------------------------

   Get the actual size of the last step that was taken

   @return step size

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceVariableVerletDynamics::getLastStepSize( void ) const {
    return _lastStepSize;
}

/**---------------------------------------------------------------------------------------

   Print parameters

   @param message             message

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceVariableVerletDynamics::printParameters( std::stringstream& message ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceVariableVerletDynamics::printParameters";

   // ---------------------------------------------------------------------------------------

   // print parameters

   ReferenceDynamics::printParameters( message );

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing Verlet dynamics update of coordinates
   and velocities

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceVariableVerletDynamics::update( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                          RealOpenMM** velocities,
                                          RealOpenMM** forces, RealOpenMM* masses ){

    // ---------------------------------------------------------------------------------------

    static const char* methodName      = "\nReferenceVariableVerletDynamics::update";

    static const RealOpenMM zero       =  0.0;
    static const RealOpenMM one        =  1.0;

    static int debug                   =  0;

    // ---------------------------------------------------------------------------------------

    // get work arrays

    RealOpenMM** xPrime = get2DArrayAtIndex( xPrime2D );
    RealOpenMM* inverseMasses = get1DArrayAtIndex( InverseMasses );

    // first-time-through initialization

    if( getTimeStep() == 0 ){

       std::stringstream message;
       message << methodName;
       int errors = 0;

       // invert masses

       for( int ii = 0; ii < numberOfAtoms; ii++ ){
          if( masses[ii] <= zero ){
             message << "mass at atom index=" << ii << " (" << masses[ii] << ") is <= 0" << std::endl;
             errors++;
          } else {
             inverseMasses[ii] = one/masses[ii];
          }
       }

       // exit if errors

       if( errors ){
          SimTKOpenMMLog::printError( message );
       }
    }

    // Try different step sizes until the accuracy is acceptable.

    bool success = false;
    RealOpenMM maxStepSize = 5.0f*getDeltaT();
    while (!success) {
        // Perform the integration and estimate the error.

        _lastStepSize = getDeltaT();
        RealOpenMM error = zero;
        for (int i = 0; i < numberOfAtoms; ++i) {
            for (int j = 0; j < 3; ++j) {
                RealOpenMM xref = atomCoordinates[i][j] + velocities[i][j]*getDeltaT();
                RealOpenMM vPrime = velocities[i][j] + inverseMasses[i]*forces[i][j]*getDeltaT();
                xPrime[i][j] = atomCoordinates[i][j] + vPrime*getDeltaT();
                RealOpenMM xerror = xPrime[i][j]-xref;
                error += xerror*xerror;
            }
        }
        error = SQRT(error/(numberOfAtoms*3));

        // Select a new step size.

        const RealOpenMM Safety = 0.9f, MinShrink = 0.1f;
        const RealOpenMM HysteresisLow = 0.9f, HysteresisHigh = 1.0f, ErrorOrder = 2.0f;

        RealOpenMM newStepSize = Safety*getDeltaT()*POW(getAccuracy()/error, 1.0f/ErrorOrder);
        if (newStepSize > getDeltaT()) {
            if (newStepSize < HysteresisHigh*getDeltaT())
                newStepSize = getDeltaT();
        }
        if (newStepSize < getDeltaT() && error <= getAccuracy())
            newStepSize = getDeltaT();
        newStepSize = std::min(newStepSize, maxStepSize);
        newStepSize = std::max(newStepSize, MinShrink*getDeltaT());
        if (newStepSize < getDeltaT())
            newStepSize = std::min(newStepSize, HysteresisLow*getDeltaT());
        success = (newStepSize >= getDeltaT());
        if (success) {
            ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
            if (referenceConstraintAlgorithm)
                success = (referenceConstraintAlgorithm->apply(numberOfAtoms, atomCoordinates, xPrime, inverseMasses) != ErrorReturn);
            if (!success) {
                newStepSize *= 0.5f;
                maxStepSize = newStepSize;
            }
        }
        setDeltaT(newStepSize);
   }

   // Update the positions and velocities.

   RealOpenMM velocityScale = static_cast<RealOpenMM>( 1.0/_lastStepSize );
   for (int i = 0; i < numberOfAtoms; ++i) {
       for (int j = 0; j < 3; ++j) {
           velocities[i][j] = velocityScale*(xPrime[i][j] - atomCoordinates[i][j]);
           atomCoordinates[i][j] = xPrime[i][j];
       }
   }

   incrementTimeStep();

   return ReferenceDynamics::DefaultReturn;

}


