
/* Portions copyright (c) 2006-2008 Stanford University and Simbios.
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
#include "ReferenceBrownianDynamics.h"

/**---------------------------------------------------------------------------------------

   ReferenceBrownianDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param friction       friction coefficient
   @param temperature    temperature

   --------------------------------------------------------------------------------------- */

ReferenceBrownianDynamics::ReferenceBrownianDynamics( int numberOfAtoms,
                                                          RealOpenMM deltaT, RealOpenMM friction,
                                                          RealOpenMM temperature ) : 
           ReferenceDynamics( numberOfAtoms, deltaT, temperature ), friction( friction ) {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceBrownianDynamics::ReferenceBrownianDynamics";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   // ---------------------------------------------------------------------------------------

   if( friction <= zero ){

      std::stringstream message;
      message << methodName;
      message << " input frction value=" << friction << " is invalid -- setting to 1.";
      SimTKOpenMMLog::printError( message );

      this->friction = one;
     
   }
   allocate2DArrays( numberOfAtoms, 3, Max2DArrays );
   allocate1DArrays( numberOfAtoms, Max1DArrays );
   
}

/**---------------------------------------------------------------------------------------

   ReferenceBrownianDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceBrownianDynamics::~ReferenceBrownianDynamics( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceBrownianDynamics::~ReferenceBrownianDynamics";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Get the friction coefficient

   @return friction

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceBrownianDynamics::getFriction( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceBrownianDynamics::getFriction";

   // ---------------------------------------------------------------------------------------

   return friction;
}

/**---------------------------------------------------------------------------------------

   Print parameters

   @param message             message

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceBrownianDynamics::printParameters( std::stringstream& message ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceBrownianDynamics::printParameters";

   // ---------------------------------------------------------------------------------------

   // print parameters

   ReferenceDynamics::printParameters( message );
   message << " friction=" << getFriction();
   message << " T=" << getTemperature();

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing Brownian dynamics update of coordinates
   and velocities

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceBrownianDynamics::update( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                          RealOpenMM** velocities,
                                          RealOpenMM** forces, RealOpenMM* masses ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceBrownianDynamics::update";

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
   
   // Perform the integration.
   
   const RealOpenMM noiseAmplitude = sqrt(2.0*BOLTZ*getTemperature()*getDeltaT()/getFriction());
   const RealOpenMM forceScale = getDeltaT()/getFriction();
   for (int i = 0; i < numberOfAtoms; ++i) {
       for (int j = 0; j < 3; ++j) {
           xPrime[i][j] = atomCoordinates[i][j] + forceScale*forces[i][j] + noiseAmplitude*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
       }
   }
   ReferenceShakeAlgorithm* referenceShakeAlgorithm = getReferenceShakeAlgorithm();
   if( referenceShakeAlgorithm )
      referenceShakeAlgorithm->applyShake( numberOfAtoms, atomCoordinates, xPrime, inverseMasses );
   
   // Update the positions and velocities.
   
   RealOpenMM velocityScale = 1.0/getDeltaT();
   for (int i = 0; i < numberOfAtoms; ++i) {
       for (int j = 0; j < 3; ++j) {
           velocities[i][j] = velocityScale*(xPrime[i][j] - atomCoordinates[i][j]);
           atomCoordinates[i][j] = xPrime[i][j];
       }
   }

   incrementTimeStep();

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Write state

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses
   @param state               0 if initial state; otherwise nonzero
   @param baseFileName        base file name

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceBrownianDynamics::writeState( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                             RealOpenMM** velocities,
                                             RealOpenMM** forces, RealOpenMM* masses,
                                             int state, const std::string& baseFileName ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceBrownianDynamics::writeState";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;
   static const int threeI            =  3;

   // ---------------------------------------------------------------------------------------

   std::stringstream stateFileName;

   stateFileName << baseFileName;
   stateFileName << "_Step" << getTimeStep();
   // stateFileName << "_State" << state;
   stateFileName << ".txt";

   // ---------------------------------------------------------------------------------------

   // open file -- return if unsuccessful

   FILE* stateFile = NULL;
#ifdef WIN32
   fopen_s( &stateFile, stateFileName.str().c_str(), "w" );
#else
   stateFile = fopen( stateFileName.str().c_str(), "w" );
#endif

   // ---------------------------------------------------------------------------------------

   // diagnostics

   if( stateFile != NULL ){
      std::stringstream message;
      message << methodName;
      message << " Opened file=<" << stateFileName.str() << ">.\n";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName;
      message << " could not open file=<" << stateFileName.str() << "> -- abort output.\n";
      SimTKOpenMMLog::printMessage( message );
      return ReferenceDynamics::ErrorReturn;
   }   

   // ---------------------------------------------------------------------------------------

   StringVector scalarNameI;
   IntVector scalarI;

   StringVector scalarNameR;
   RealOpenMMVector scalarR;

   StringVector scalarNameR1;
   RealOpenMMPtrVector scalarR1;

   StringVector scalarNameR2;
   RealOpenMMPtrPtrVector scalarR2;

   scalarI.push_back( getNumberOfAtoms() );
   scalarNameI.push_back( "Atoms" );

   scalarI.push_back( getTimeStep() );
   scalarNameI.push_back( "Timestep" );

   if( state == 0 || state == -1 ){

      scalarR.push_back( getDeltaT() );
      scalarNameR.push_back( "delta_t" );

      scalarR.push_back( getTemperature() );
      scalarNameR.push_back( "T" );

      scalarR.push_back( getFriction() );
      scalarNameR.push_back( "Friction" );

      scalarR1.push_back( masses );
      scalarNameR1.push_back( "mass" );

      scalarR2.push_back( atomCoordinates );
      scalarNameR2.push_back( "coord" );

      scalarR2.push_back( velocities );
      scalarNameR2.push_back( "velocities" );

      scalarR2.push_back( forces );
      scalarNameR2.push_back( "forces" );

      if( state == -1 ){

         RealOpenMM** xPrime = get2DArrayAtIndex( xPrime2D );

         scalarR2.push_back( xPrime );
         scalarNameR2.push_back( "xPrime" );
      }
      
   } else {

      scalarR2.push_back( atomCoordinates );
      scalarNameR2.push_back( "coord" );

      scalarR2.push_back( velocities );
      scalarNameR2.push_back( "velocities" );

   }

   writeStateToFile( stateFile, scalarNameI, scalarI, scalarNameR, scalarR, getNumberOfAtoms(), scalarNameR1, scalarR1, threeI, scalarNameR2, scalarR2 ); 

   (void) fclose( stateFile );

   return ReferenceDynamics::DefaultReturn;

}
