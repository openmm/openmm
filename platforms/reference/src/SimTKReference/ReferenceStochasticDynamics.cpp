
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
#include "ReferenceStochasticDynamics.h"

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

   // insure tau is not zero -- if it is print warning message

   if( _tau == zero ){

      std::stringstream message;
      message << methodName;
      message << " input tau value=" << tau << " is invalid -- setting to 1.";
      SimTKOpenMMLog::printError( message );

      _tau = one;
     
   }
   _setFixedParameters( );

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

   Set fixed parameters

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceStochasticDynamics::_setFixedParameters( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceStochasticDynamics::_setFixedParameters";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;
   static const RealOpenMM two        =  2.0;
   static const RealOpenMM three      =  3.0;
   static const RealOpenMM four       =  4.0;
   static const RealOpenMM half       =  0.5;

   // ---------------------------------------------------------------------------------------

   _fixedParameters[GDT]      = getDeltaT()/getTau();
   _fixedParameters[EPH]      = EXP(  half*_fixedParameters[GDT] );
   _fixedParameters[EMH]      = EXP( -half*_fixedParameters[GDT] );
   _fixedParameters[EM]       = EXP(      -_fixedParameters[GDT] );
   _fixedParameters[EP]       = EXP(       _fixedParameters[GDT] );

   if( _fixedParameters[GDT] >= (RealOpenMM) 0.1 ){

      RealOpenMM term1        = _fixedParameters[EPH] - one;
                 term1       *= term1;
      _fixedParameters[B]     = _fixedParameters[GDT]*(_fixedParameters[EP] - one) - four*term1;

      _fixedParameters[C]     = _fixedParameters[GDT] - three + four*_fixedParameters[EMH] - _fixedParameters[EM];
      _fixedParameters[D]     = two - _fixedParameters[EPH] - _fixedParameters[EMH];

    } else {

      // this has not been debugged

      RealOpenMM term1        = half*_fixedParameters[GDT];
      RealOpenMM term2        = term1*term1;
      RealOpenMM term4        = term2*term2;

      RealOpenMM third        = (RealOpenMM) ( 1.0/3.0 );
      RealOpenMM o7_9         = (RealOpenMM) ( 7.0/9.0 );
      RealOpenMM o1_12        = (RealOpenMM) ( 1.0/12.0 );
      RealOpenMM o17_90       = (RealOpenMM) ( 17.0/90.0 );
      RealOpenMM o7_30        = (RealOpenMM) ( 7.0/30.0 );
      RealOpenMM o31_1260     = (RealOpenMM) ( 31.0/1260.0 );
      RealOpenMM o_360        = (RealOpenMM) ( 1.0/360.0 );

      _fixedParameters[B]     = term4*( third  + term1*( third + term1*( o17_90 + term1*o7_9 )));
      _fixedParameters[C]     = term2*term1*( two*third + term1*( -half + term1*( o7_30 + term1*(-o1_12 + term1*o31_1260 ))));
      _fixedParameters[D]     = term2*( -one + term2*(-o1_12 - term2*o_360));
   }    

   RealOpenMM kT        = ((RealOpenMM) BOLTZ)*getTemperature();

   _fixedParameters[V]  = SQRT( kT*( one - _fixedParameters[EM]) );
   _fixedParameters[X]  = getTau()*SQRT( kT*_fixedParameters[C] );
   _fixedParameters[Yv] = SQRT( kT*_fixedParameters[B]/_fixedParameters[C] );
   _fixedParameters[Yx] = getTau()*SQRT( kT*_fixedParameters[B]/(one - _fixedParameters[EM]) );

   return ReferenceDynamics::DefaultReturn;

};

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

   Get array of fixed parameters indexed by 'FixedParameters' enums

   @return array

   --------------------------------------------------------------------------------------- */
   
const RealOpenMM* ReferenceStochasticDynamics::getFixedParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceStochasticDynamics::getFixedParameters";

   // ---------------------------------------------------------------------------------------

   return _fixedParameters;
}

/**---------------------------------------------------------------------------------------

   Print parameters

   @param message             message

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceStochasticDynamics::printParameters( std::stringstream& message ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceStochasticDynamics::printParameters";
   static const char* parameterNames[MaxFixedParameters] = { "gdt", "ep", "eph", "emh", "em", "B", "C", "D",
                                                             "V", "X", "Yv", "Yx" };

   // ---------------------------------------------------------------------------------------

   // print parameters

   ReferenceDynamics::printParameters( message );
   message << " tau=" << getTau();
   message << " T=" << getTemperature();
   int cut = 3;
   for( int ii = 0; ii < MaxFixedParameters; ii++ ){
      message << " " << parameterNames[ii] << "=" << _fixedParameters[ii];
      if( cut++ > 5 ){
         cut = 0;
         message << std::endl;
      }
   }

   return ReferenceDynamics::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   First SD update; based on code in update.c do_update_sd() Gromacs 3.1.4

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param inverseMasses       inverse atom masses
   @param xPrime              xPrime
   @param oldVelocities       previous velocities
   @param xVector             xVector
   @param vVector             vVector

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceStochasticDynamics::updatePart1( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                              RealOpenMM** velocities,
                                              RealOpenMM** forces, RealOpenMM* inverseMasses,
                                              RealOpenMM** xPrime, RealOpenMM** oldVelocities,
                                              RealOpenMM** xVector, RealOpenMM** vVector ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceStochasticDynamics::updatePart1";

   static const RealOpenMM one    =  1.0;
   static int debug               = 0;

   // ---------------------------------------------------------------------------------------

   // perform first update

   const RealOpenMM*  fixedParameters = getFixedParameters();
   RealOpenMM   tau                   = getTau();
   RealOpenMM   fix1                  = tau*(fixedParameters[EPH] - fixedParameters[EMH]);

   for( int ii = 0; ii < numberOfAtoms; ii++ ){

      RealOpenMM sqrtInvMass = SQRT( inverseMasses[ii] );

      for( int jj = 0; jj < 3; jj++ ){

         oldVelocities[ii][jj] = velocities[ii][jj];

         RealOpenMM Vmh        = xVector[ii][jj]*fixedParameters[D]/(tau*fixedParameters[C]) +
                                 sqrtInvMass*fixedParameters[Yv]*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();

           vVector[ii][jj]     = sqrtInvMass*fixedParameters[V]*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
          
           velocities[ii][jj]  = oldVelocities[ii][jj]*fixedParameters[EM] +
                                 inverseMasses[ii]*forces[ii][jj]*tau*( one - fixedParameters[EM]) +
                                 vVector[ii][jj] - fixedParameters[EM]*Vmh;
            
           xPrime[ii][jj]      = atomCoordinates[ii][jj] + velocities[ii][jj]*fix1;
      }
   }

   // diagnostics

   if( debug ){
      int maxPrint = 5;
      std::stringstream message;
      message << methodName << " Post SD1 atoms=" << numberOfAtoms << "\n";
      RealOpenMM** oldVelocities   = get2DArrayAtIndex( OldV );
      for( int ii = 0; ii < maxPrint; ii++ ){
         message << " mI=" << inverseMasses[ii];
        	RealOpenMM sqrtInvMass = SQRT( inverseMasses[ii] );
         message << " sdpc[" <<	sqrtInvMass*fixedParameters[Yv] << " " << sqrtInvMass*fixedParameters[V];
         message << " x[";
         SimTKOpenMMUtilities::formatRealStringStream( message, atomCoordinates[ii], 3, one );
         message << "] xp[";
         SimTKOpenMMUtilities::formatRealStringStream( message, xPrime[ii], 3, one );
         message << "] v[";
         SimTKOpenMMUtilities::formatRealStringStream( message, velocities[ii], 3, one );
         message << "] vV[";
         SimTKOpenMMUtilities::formatRealStringStream( message, vVector[ii], 3, one );
         message << "] ov[";
         SimTKOpenMMUtilities::formatRealStringStream( message, oldVelocities[ii], 3, one );
         message << "] xV[";
         message << "] f[";
         SimTKOpenMMUtilities::formatRealStringStream( message, forces[ii], 3, one );
         message << "] xV[";
         SimTKOpenMMUtilities::formatRealStringStream( message, xVector[ii], 3, one );
         message << "]\n";
      }
      SimTKOpenMMLog::printMessage( message );
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
                                              RealOpenMM** xPrime, RealOpenMM** oldVelocities,
                                              RealOpenMM** xVector, RealOpenMM** vVector ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceStochasticDynamics::updatePart2";

   static const RealOpenMM one    =  1.0;
   static int debug               = 0;

   // ---------------------------------------------------------------------------------------

   // perform second update

   const RealOpenMM*  fixedParameters = getFixedParameters();
   RealOpenMM   tau                   = getTau();
   RealOpenMM   fix1                  = tau*(fixedParameters[EPH] - fixedParameters[EMH]);
                fix1                  = one/fix1;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){

      RealOpenMM sqrtInvMass = SQRT( inverseMasses[ii] );

      for( int jj = 0; jj < 3; jj++ ){

           velocities[ii][jj]  = (xPrime[ii][jj] - atomCoordinates[ii][jj])*fix1;
         RealOpenMM Xmh        = vVector[ii][jj]*tau*fixedParameters[D]/(fixedParameters[EM] - one) +
                                 sqrtInvMass*fixedParameters[Yx]*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();

           xVector[ii][jj]     = sqrtInvMass*fixedParameters[X]*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
            
           xPrime[ii][jj]     += xVector[ii][jj] - Xmh;
      }
   }

   // diagnostics

   if( debug ){
      int maxPrint = 5;
      std::stringstream message;
      message << methodName << " Post SD2 atoms=" << numberOfAtoms << "\n";
      RealOpenMM** oldVelocities   = get2DArrayAtIndex( OldV );
      for( int ii = 0; ii < maxPrint; ii++ ){
         message << " mI=" << inverseMasses[ii];
        	RealOpenMM sqrtInvMass = SQRT( inverseMasses[ii] );
         message << " sdpc[" <<	sqrtInvMass*fixedParameters[Yx] << " " << sqrtInvMass*fixedParameters[X];
         message << " x[";
         message << " x[";
         SimTKOpenMMUtilities::formatRealStringStream( message, atomCoordinates[ii], 3, one );
         message << "] xp[";
         SimTKOpenMMUtilities::formatRealStringStream( message, xPrime[ii], 3, one );
         message << "] v[";
         SimTKOpenMMUtilities::formatRealStringStream( message, velocities[ii], 3, one );
         message << "] vV[";
         SimTKOpenMMUtilities::formatRealStringStream( message, vVector[ii], 3, one );
         message << "] oV[";
         SimTKOpenMMUtilities::formatRealStringStream( message, oldVelocities[ii], 3, one );
         message << "] xV[";
         SimTKOpenMMUtilities::formatRealStringStream( message, xVector[ii], 3, one );
         message << "]\n";
      }
      SimTKOpenMMLog::printMessage( message );
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
                                          RealOpenMM** velocities,
                                          RealOpenMM** forces, RealOpenMM* masses ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceStochasticDynamics::update";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   static int debug                   =  0;

   // ---------------------------------------------------------------------------------------

   // get work arrays

   RealOpenMM** xPrime          = get2DArrayAtIndex( xPrime2D );
   RealOpenMM** oldVelocities   = get2DArrayAtIndex( OldV );
   RealOpenMM** xVector         = get2DArrayAtIndex( X2D );
   RealOpenMM** vVector         = get2DArrayAtIndex( V2D );

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

      // set xVector 

      const RealOpenMM*  fixedParameters = getFixedParameters();
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         RealOpenMM sqrtInverseMass = SQRT( inverseMasses[ii] )*fixedParameters[X];
         for( int jj = 0; jj < 3; jj++ ){
            xVector[ii][jj] = sqrtInverseMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
         }
      }

      // exit if errors

      if( errors ){
         SimTKOpenMMLog::printError( message );
      }
   }

   // 1st update

   updatePart1( numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses,
                xPrime, oldVelocities, xVector, vVector );

   //writeState( numberOfAtoms, atomCoordinates, velocities, forces, masses, -1 , "Sd1" );

   ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
   if( referenceConstraintAlgorithm ){

/*
      std::stringstream message;
      message << methodName;
      message << " calling constrain1\n";
      SimTKOpenMMLog::printMessage( message );
*/

      referenceConstraintAlgorithm->apply( numberOfAtoms, atomCoordinates, xPrime,
                                           inverseMasses );

   }

   // 2nd update

   if( debug ){
      int maxPrint = 5;
      std::stringstream message;
      message << methodName << " Pre SD2 atoms=" << numberOfAtoms << std::endl;
      RealOpenMM** oldVelocities   = get2DArrayAtIndex( OldV );
      for( int ii = 0; ii < maxPrint; ii++ ){
         message << " x[";
         SimTKOpenMMUtilities::formatRealStringStream( message, atomCoordinates[ii], 3, one );
         message << "] xp[";
         SimTKOpenMMUtilities::formatRealStringStream( message, xPrime[ii], 3, one );
         message << "] v[";
         SimTKOpenMMUtilities::formatRealStringStream( message, velocities[ii], 3, one );
         message << "] ov[";
         SimTKOpenMMUtilities::formatRealStringStream( message, oldVelocities[ii], 3, one );
         message << "]\n";
      }
      SimTKOpenMMLog::printMessage( message );
   }

   updatePart2( numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses,
                xPrime, oldVelocities, xVector, vVector );

   if( debug ){
      int maxPrint = 5;
      std::stringstream message;
      message << methodName << " Post SD2 atoms=" << numberOfAtoms << "\n";
      RealOpenMM** oldVelocities   = get2DArrayAtIndex( OldV );
      for( int ii = 0; ii < maxPrint; ii++ ){
         message << " x[";
         SimTKOpenMMUtilities::formatRealStringStream( message, atomCoordinates[ii], 3, one );
         message << "] xp[";
         SimTKOpenMMUtilities::formatRealStringStream( message, xPrime[ii], 3, one );
         message << "] v[";
         SimTKOpenMMUtilities::formatRealStringStream( message, velocities[ii], 3, one );
         message << "] ov[";
         SimTKOpenMMUtilities::formatRealStringStream( message, oldVelocities[ii], 3, one );
         message << "]\n";
      }
      SimTKOpenMMLog::printMessage( message );
   }

   if( referenceConstraintAlgorithm ){

/*
      std::stringstream message;
      message << methodName;
      message << " calling constrain2\n";
      SimTKOpenMMLog::printMessage( message );
*/

      referenceConstraintAlgorithm->apply( numberOfAtoms, atomCoordinates, xPrime,
                                           inverseMasses );
   }

   // copy xPrime -> atomCoordinates

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      atomCoordinates[ii][0] = xPrime[ii][0];
      atomCoordinates[ii][1] = xPrime[ii][1];
      atomCoordinates[ii][2] = xPrime[ii][2];
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

int ReferenceStochasticDynamics::writeState( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                             RealOpenMM** velocities,
                                             RealOpenMM** forces, RealOpenMM* masses,
                                             int state, const std::string& baseFileName ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceStochasticDynamics::writeState";

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

      scalarR.push_back( getTau() );
      scalarNameR.push_back( "Tau" );

      scalarR1.push_back( masses );
      scalarNameR1.push_back( "mass" );

      scalarR2.push_back( atomCoordinates );
      scalarNameR2.push_back( "coord" );

      scalarR2.push_back( velocities );
      scalarNameR2.push_back( "velocities" );

      scalarR2.push_back( forces );
      scalarNameR2.push_back( "forces" );

      if( state == -1 ){

         RealOpenMM** xPrime          = get2DArrayAtIndex( xPrime2D );
         RealOpenMM** oldVelocities   = get2DArrayAtIndex( OldV );
         RealOpenMM** xVector         = get2DArrayAtIndex( X2D );
         RealOpenMM** vVector         = get2DArrayAtIndex( V2D );

         scalarR2.push_back( xPrime );
         scalarNameR2.push_back( "xPrime" );

         scalarR2.push_back( oldVelocities);
         scalarNameR2.push_back( "vold" );

         scalarR2.push_back( xVector );
         scalarNameR2.push_back( "xVector" );

         scalarR2.push_back( vVector );
         scalarNameR2.push_back( "vVector" );
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
