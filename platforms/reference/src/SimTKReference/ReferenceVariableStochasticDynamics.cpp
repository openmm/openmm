
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
#include "ReferenceVariableStochasticDynamics.h"

#include <cstdio>

/**---------------------------------------------------------------------------------------

   ReferenceVariableStochasticDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param tau            viscosity(?)
   @param temperature    temperature
   @param accuracy       required accuracy

   --------------------------------------------------------------------------------------- */

ReferenceVariableStochasticDynamics::ReferenceVariableStochasticDynamics( int numberOfAtoms,
                                                          RealOpenMM tau, RealOpenMM temperature,
                                                          RealOpenMM accuracy ) :
           ReferenceDynamics(numberOfAtoms, 0.0f, temperature), _tau(tau), _accuracy(accuracy) {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceVariableStochasticDynamics::ReferenceVariableStochasticDynamics";

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

   allocate2DArrays( numberOfAtoms, 3, Max2DArrays );
   allocate1DArrays( numberOfAtoms, Max1DArrays );

}

/**---------------------------------------------------------------------------------------

   ReferenceVariableStochasticDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceVariableStochasticDynamics::~ReferenceVariableStochasticDynamics( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceVariableStochasticDynamics::~ReferenceVariableStochasticDynamics";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Get the required accuracy

   @return accuracy

 --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceVariableStochasticDynamics::getAccuracy( void ) const {
    return _accuracy;
}

/**---------------------------------------------------------------------------------------

   Set the required accuracy

 --------------------------------------------------------------------------------------- */

void ReferenceVariableStochasticDynamics::setAccuracy( RealOpenMM accuracy ) {
    _accuracy = accuracy;
}

/**---------------------------------------------------------------------------------------

   Set fixed parameters

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceVariableStochasticDynamics::_setFixedParameters( RealOpenMM timeStep, RealOpenMM prevTimeStep ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceVariableStochasticDynamics::_setFixedParameters";

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;
   static const RealOpenMM two        =  2.0;
   static const RealOpenMM three      =  3.0;
   static const RealOpenMM four       =  4.0;
   static const RealOpenMM half       =  0.5;

   // ---------------------------------------------------------------------------------------

   _fixedParameters[GDT]      = timeStep/getTau();
   _fixedParameters[EPH]      = EXP(  half*_fixedParameters[GDT] );
   _fixedParameters[EMH]      = EXP( -half*_fixedParameters[GDT] );
   _fixedParameters[EM]       = EXP(      -_fixedParameters[GDT] );
   _fixedParameters[EM_V]     = EXP(      -half*(timeStep+prevTimeStep)/getTau() );
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

RealOpenMM ReferenceVariableStochasticDynamics::getTau( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceVariableStochasticDynamics::getTau";

   // ---------------------------------------------------------------------------------------

   return _tau;
}

/**---------------------------------------------------------------------------------------

   Get array of fixed parameters indexed by 'FixedParameters' enums

   @return array

   --------------------------------------------------------------------------------------- */

const RealOpenMM* ReferenceVariableStochasticDynamics::getFixedParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceVariableStochasticDynamics::getFixedParameters";

   // ---------------------------------------------------------------------------------------

   return _fixedParameters;
}

/**---------------------------------------------------------------------------------------

   Print parameters

   @param message             message

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceVariableStochasticDynamics::printParameters( std::stringstream& message ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceVariableStochasticDynamics::printParameters";
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
   @param masses              atom masses
   @param inverseMasses       inverse atom masses
   @param xPrime              xPrime
   @param oldVelocities       previous velocities
   @param xVector             xVector
   @param vVector             vVector
   @param maxStepSize         maximum time step

   @return ReferenceDynamics::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceVariableStochasticDynamics::updatePart1( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                              RealOpenMM** velocities,
                                              RealOpenMM** forces, RealOpenMM* masses, RealOpenMM* inverseMasses,
                                              RealOpenMM** xPrime, RealOpenMM** oldVelocities,
                                              RealOpenMM** xVector, RealOpenMM** vVector,
                                              RealOpenMM maxStepSize ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceVariableStochasticDynamics::updatePart1";

   static const RealOpenMM zero   =  0.0;
   static const RealOpenMM one    =  1.0;
   static int debug               = 0;

   // ---------------------------------------------------------------------------------------


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

   // Select the step size to use
    RealOpenMM error = zero;
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
   _setFixedParameters(newStepSize, getDeltaT());
    setDeltaT(newStepSize);
 
   if( getTimeStep() == 0 ){

      // Initialize xVector

      const RealOpenMM*  fixedParameters = getFixedParameters();
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         RealOpenMM sqrtInverseMass = SQRT( inverseMasses[ii] )*fixedParameters[X];
         for( int jj = 0; jj < 3; jj++ ){
            xVector[ii][jj] = sqrtInverseMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
         }
      }
   }

    // perform first update

   const RealOpenMM*  fixedParameters = getFixedParameters();
   RealOpenMM   tau                   = getTau();
   RealOpenMM   fix1                  = tau*(fixedParameters[EPH] - fixedParameters[EMH]);
   if (fix1 == zero)
       fix1 = getDeltaT();

   for( int ii = 0; ii < numberOfAtoms; ii++ ){

      RealOpenMM sqrtInvMass = SQRT( inverseMasses[ii] );

      for( int jj = 0; jj < 3; jj++ ){

         oldVelocities[ii][jj] = velocities[ii][jj];

         RealOpenMM Vmh        = xVector[ii][jj]*fixedParameters[D]/(tau*fixedParameters[C]) +
                                 sqrtInvMass*fixedParameters[Yv]*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();

           vVector[ii][jj]     = sqrtInvMass*fixedParameters[V]*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();

           RealOpenMM vPrime   = oldVelocities[ii][jj]*fixedParameters[EM_V] +
                                 inverseMasses[ii]*forces[ii][jj]*tau*(one - fixedParameters[EM_V]) +
                                 vVector[ii][jj] - fixedParameters[EM]*Vmh;

           xPrime[ii][jj]      = atomCoordinates[ii][jj] + vPrime*fix1;
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

int ReferenceVariableStochasticDynamics::updatePart2( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                              RealOpenMM** velocities,
                                              RealOpenMM** forces, RealOpenMM* inverseMasses,
                                              RealOpenMM** xPrime, RealOpenMM** oldVelocities,
                                              RealOpenMM** xVector, RealOpenMM** vVector ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nReferenceVariableStochasticDynamics::updatePart2";

   static const RealOpenMM zero   =  0.0;
   static const RealOpenMM one    =  1.0;
   static int debug               = 0;

   // ---------------------------------------------------------------------------------------

   // perform second update

   const RealOpenMM*  fixedParameters = getFixedParameters();
   RealOpenMM   tau                   = getTau();
   RealOpenMM   fix1                  = tau*(fixedParameters[EPH] - fixedParameters[EMH]);
   if (fix1 == zero)
       fix1 = getDeltaT();
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

int ReferenceVariableStochasticDynamics::update( int numberOfAtoms, RealOpenMM** atomCoordinates,
                                          RealOpenMM** velocities,
                                          RealOpenMM** forces, RealOpenMM* masses, RealOpenMM maxStepSize ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceVariableStochasticDynamics::update";

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

   // 1st update

   updatePart1( numberOfAtoms, atomCoordinates, velocities, forces, masses, inverseMasses,
                xPrime, oldVelocities, xVector, vVector, maxStepSize );

   ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
   if( referenceConstraintAlgorithm ){
      referenceConstraintAlgorithm->apply( numberOfAtoms, atomCoordinates, xPrime,
                                           inverseMasses );
   }

   // 2nd update

   updatePart2( numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses,
                xPrime, oldVelocities, xVector, vVector );

   if( referenceConstraintAlgorithm ){
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
