
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

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceStochasticDynamics.h"

#include <cstdio>

/**---------------------------------------------------------------------------------------

   ReferenceStochasticDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param tau            viscosity(?)
   @param temperature    temperature

   --------------------------------------------------------------------------------------- */

ReferenceStochasticDynamics::ReferenceStochasticDynamics( int numberOfAtoms,
                                                          RealOpenMM deltaT, RealOpenMM tau,
                                                          RealOpenMM temperature ) : 
           ReferenceDynamics( numberOfAtoms, deltaT, temperature ), _tau( tau ) {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceStochasticDynamics::ReferenceStochasticDynamics";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   // ---------------------------------------------------------------------------------------

   // ensure tau is not zero -- if it is print warning message

   if( _tau == zero ){

      std::stringstream message;
      message << methodName;
      message << " input tau value=" << tau << " is invalid -- setting to 1.";
      SimTKOpenMMLog::printError( message );

      _tau = one;
     
   }

   allocate2DArrays( numberOfAtoms, 3, Max2DArrays );
   allocate1DArrays( numberOfAtoms, Max1DArrays );
   
}

/**---------------------------------------------------------------------------------------

   ReferenceStochasticDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceStochasticDynamics::~ReferenceStochasticDynamics( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceStochasticDynamics::~ReferenceStochasticDynamics";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Get tau

   @return tau

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceStochasticDynamics::getTau( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceStochasticDynamics::getTau";

   // ---------------------------------------------------------------------------------------

   return _tau;
}

/**---------------------------------------------------------------------------------------

   First SD update; based on code in update.c do_update_sd() Gromacs 3.1.4

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param inverseMasses       inverse atom masses
   @param xPrime              xPrime

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceStochasticDynamics::updatePart1( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                              RealOpenMM** velocities,
                                              RealOpenMM** forces, RealOpenMM* inverseMasses,
                                              RealOpenMM** xPrime ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName  = "\nReferenceStochasticDynamics::updatePart1";

   // ---------------------------------------------------------------------------------------

   // perform first update

   RealOpenMM tau = getTau();
   const RealOpenMM vscale = EXP(-getDeltaT()/tau);
   const RealOpenMM fscale = (1-vscale)*tau;
   const RealOpenMM kT = BOLTZ*getTemperature();
   const RealOpenMM noisescale = SQRT(2*kT/tau)*SQRT(0.5*(1-vscale*vscale)*tau);

   for (int ii = 0; ii < numberOfAtoms; ii++) {
      RealOpenMM sqrtInvMass = SQRT(inverseMasses[ii]);
      for (int jj = 0; jj < 3; jj++) {
           velocities[ii][jj]  = vscale*velocities[ii][jj] + fscale*inverseMasses[ii]*forces[ii][jj] + noisescale*sqrtInvMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
      }
   }

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Second update; based on code in update.c do_update_sd() w/ bFirstHalf = false in Gromacs 3.1.4

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceStochasticDynamics::updatePart2( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                              RealOpenMM** velocities,
                                              RealOpenMM** forces, RealOpenMM* inverseMasses,
                                              RealOpenMM** xPrime ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName  = "\nReferenceStochasticDynamics::updatePart2";

   // ---------------------------------------------------------------------------------------

   // perform second update

   for (int ii = 0; ii < numberOfAtoms; ii++) {
      for (int jj = 0; jj < 3; jj++) {
         xPrime[ii][jj] = atomCoordinates[ii][jj]+getDeltaT()*velocities[ii][jj];
      }
   }

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing stochastic dynamics update of coordinates
   and velocities

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceStochasticDynamics::update( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                          RealOpenMM** velocities, RealOpenMM** forces, RealOpenMM* masses ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceStochasticDynamics::update";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   // ---------------------------------------------------------------------------------------

   // get work arrays

   RealOpenMM** xPrime          = get2DArrayAtIndex( xPrime2D );
   RealOpenMM* inverseMasses    = get1DArrayAtIndex( InverseMasses );

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

   // 1st update

   updatePart1( numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses, xPrime );

   // 2nd update

   updatePart2( numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses, xPrime );

   ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
   if( referenceConstraintAlgorithm ){
      referenceConstraintAlgorithm->apply( numberOfAtoms, atomCoordinates, xPrime, inverseMasses );
   }

   // copy xPrime -> atomCoordinates

   RealOpenMM invStepSize = 1.0/getDeltaT();
   for (int i = 0; i < numberOfAtoms; ++i)
       for (int j = 0; j < 3; ++j) {
           velocities[i][j] = invStepSize*(xPrime[i][j]-atomCoordinates[i][j]);
           atomCoordinates[i][j] = xPrime[i][j];
       }

   incrementTimeStep();

   return ReferenceDynamics::DefaultReturn;

}
