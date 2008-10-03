/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include <sstream>
#include "BrookLangevinDynamics.h"
#include "BrookPlatform.h"
#include "OpenMMException.h"
#include "BrookStreamImpl.h"
#include "gpu/kshakeh.h"
#include "gpu/kupdatesd.h"
#include "gpu/kcommon.h"

// use random number generator

#include "../../reference/src/SimTKUtilities/SimTKOpenMMUtilities.h"

using namespace OpenMM;
using namespace std;

/** 
 *
 * Constructor
 * 
 */

BrookLangevinDynamics::BrookLangevinDynamics( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookLangevinDynamics::BrookLangevinDynamics";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat one                     = (BrookOpenMMFloat)  1.0;
   BrookOpenMMFloat oneMinus                = (BrookOpenMMFloat) -1.0;

// ---------------------------------------------------------------------------------------

   _numberOfAtoms             = -1;

   // mark stream dimension variables as unset

   _sdAtomStreamWidth         = -1;
   _sdAtomStreamHeight        = -1;
   _sdAtomStreamSize          = -1;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _sdStreams[ii]   = NULL;
   }

   for( int ii = 0; ii < MaxDerivedParameters; ii++ ){
      _derivedParameters[ii]   = oneMinus;
   }

   _temperature = oneMinus;
   _stepSize    = oneMinus;
   _tau         = oneMinus;

   // setup inverse sqrt masses

   _inverseSqrtMasses = NULL;

   // set randomNumber seed 

   _randomNumberSeed  = 1393;

   //_randomNumberSeed = randomNumberSeed ? randomNumberSeed : 1393;
   //SimTKOpenMMUtilities::setRandomNumberSeed( randomNumberSeed );
}   
 
/** 
 * Destructor
 * 
 */

BrookLangevinDynamics::~BrookLangevinDynamics( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookLangevinDynamics::~BrookLangevinDynamics";

// ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      delete _sdStreams[ii];
   }

   delete[] _inverseSqrtMasses;

}

/** 
 * Get tau
 * 
 * @return  tau
 *
 */

BrookOpenMMFloat BrookLangevinDynamics::getTau( void ) const {
   return _tau;
}

/** 
 * Get friction
 * 
 * @return  friction
 *
 */

BrookOpenMMFloat BrookLangevinDynamics::getFriction( void ) const {
   static const BrookOpenMMFloat zero = (BrookOpenMMFloat) 0.0; 
   static const BrookOpenMMFloat one  = (BrookOpenMMFloat) 1.0; 
   return ( (_tau == zero) ? zero : (one/_tau) );
}

/** 
 * Get temperature
 * 
 * @return  temperature
 *
 */

BrookOpenMMFloat BrookLangevinDynamics::getTemperature( void ) const {
   return _temperature;
}

/** 
 * Get stepSize
 * 
 * @return  stepSize
 *
 */

BrookOpenMMFloat BrookLangevinDynamics::getStepSize( void ) const {
   return _stepSize;
}

/** 
 * Set tau
 * 
 * @param tau   new tau value
 *
 * @return      DefaultReturnValue
 *
 */

int BrookLangevinDynamics::_setTau( BrookOpenMMFloat tau ){
   _tau = tau;
   return DefaultReturnValue;
}

/** 
 * Set friction = 1/tau
 * 
 * @param friction   new friction value
 *
 * @return      DefaultReturnValue
 *
 */

int BrookLangevinDynamics::_setFriction( BrookOpenMMFloat friction ){
   _tau   = (BrookOpenMMFloat) ( (friction != 0.0) ? 1.0/friction : 0.0);
   return DefaultReturnValue;
}

/** 
 * Set temperature
 * 
 * @parameter   temperature
 *
 * @return      DefaultReturnValue
 *
 */

int BrookLangevinDynamics::_setTemperature( BrookOpenMMFloat temperature ){
   _temperature = temperature;
   return DefaultReturnValue;
}

/** 
 * Set stepSize
 * 
 * @param   stepSize
 *
 * @return      DefaultReturnValue
 *
 */

int BrookLangevinDynamics::_setStepSize( BrookOpenMMFloat stepSize ){
   _stepSize = stepSize;
   return DefaultReturnValue;
}

/** 
 * Update derived parameters
 * 
 * @return  DefaultReturnValue
 *
 * @throw   if tau too small
 *
 */

int BrookLangevinDynamics::_updateDerivedParameters( void ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nBrookLangevinDynamics::_updateDerivedParameters";

   static const BrookOpenMMFloat zero       =  0.0;
   static const BrookOpenMMFloat one        =  1.0;
   static const BrookOpenMMFloat two        =  2.0;
   static const BrookOpenMMFloat three      =  3.0;
   static const BrookOpenMMFloat four       =  4.0;
   static const BrookOpenMMFloat half       =  0.5;

   float epsilon                            = 1.0e-08f;

   // ---------------------------------------------------------------------------------------

   BrookOpenMMFloat tau         = getTau();
   BrookOpenMMFloat temperature = getTemperature();
   BrookOpenMMFloat stepSize    = getStepSize();

   if(  fabsf( float( tau ) ) < epsilon ){
      std::stringstream message;
      message << methodName << " tau=" << tau << " too small.";
      throw OpenMMException( message.str() );
   }   


   _derivedParameters[GDT]      = stepSize/tau;

   _derivedParameters[EPH]      = EXP(  half*_derivedParameters[GDT] );
   _derivedParameters[EMH]      = EXP( -half*_derivedParameters[GDT] );
   _derivedParameters[EM]       = EXP(      -_derivedParameters[GDT] );
   _derivedParameters[EP]       = EXP(       _derivedParameters[GDT] );

   if( _derivedParameters[GDT] >= (BrookOpenMMFloat) 0.1 ){

      BrookOpenMMFloat term1    = _derivedParameters[EPH] - one;
                 term1         *= term1;
      _derivedParameters[B]     = _derivedParameters[GDT]*(_derivedParameters[EP] - one) - four*term1;

      _derivedParameters[C]     = _derivedParameters[GDT] - three + four*_derivedParameters[EMH] - _derivedParameters[EM];
      _derivedParameters[D]     = two - _derivedParameters[EPH] - _derivedParameters[EMH];

    } else {

      BrookOpenMMFloat term1        = half*_derivedParameters[GDT];
      BrookOpenMMFloat term2        = term1*term1;
      BrookOpenMMFloat term4        = term2*term2;

      BrookOpenMMFloat third        = (BrookOpenMMFloat) ( 1.0/3.0 );
      BrookOpenMMFloat o7_9         = (BrookOpenMMFloat) ( 7.0/9.0 );
      BrookOpenMMFloat o1_12        = (BrookOpenMMFloat) ( 1.0/12.0 );
      BrookOpenMMFloat o17_90       = (BrookOpenMMFloat) ( 17.0/90.0 );
      BrookOpenMMFloat o7_30        = (BrookOpenMMFloat) ( 7.0/30.0 );
      BrookOpenMMFloat o31_1260     = (BrookOpenMMFloat) ( 31.0/1260.0 );
      BrookOpenMMFloat o_360        = (BrookOpenMMFloat) ( 1.0/360.0 );

      _derivedParameters[B]         = term4*( third  + term1*( third + term1*( o17_90 + term1*o7_9 )));
      _derivedParameters[C]         = term2*term1*( two*third + term1*( -half + term1*( o7_30 + term1*(-o1_12 + term1*o31_1260 ))));
      _derivedParameters[D]         = term2*( -one + term2*(-o1_12 - term2*o_360));
   }    

   BrookOpenMMFloat kT        = ((BrookOpenMMFloat) BOLTZ)*temperature;

   _derivedParameters[V]      = SQRT( kT*( one - _derivedParameters[EM]) );
   _derivedParameters[X]      = tau*SQRT( kT*_derivedParameters[C] );
   _derivedParameters[Yv]     = SQRT( kT*_derivedParameters[B]/_derivedParameters[C] );
   _derivedParameters[Yx]     = tau*SQRT( kT*_derivedParameters[B]/(one - _derivedParameters[EM]) );

   _derivedParameters[Sd1pc1] = tau*( one - _derivedParameters[EM] );
   _derivedParameters[Sd1pc2] = tau*( _derivedParameters[EPH] - _derivedParameters[EMH] );
   _derivedParameters[Sd1pc3] = _derivedParameters[D]/(tau*_derivedParameters[C] );
   _derivedParameters[Sd2pc1] = one/_derivedParameters[Sd1pc2];
   _derivedParameters[Sd2pc2] = tau*_derivedParameters[D]/( _derivedParameters[EM] - one );
       
   return DefaultReturnValue;

}

/** 
 * Update parameters -- only way parameters can be set
 * 
 * @param  temperature     temperature
 * @param  friction        friction
 * @param  step size       step size
 *
 * @return   solute dielectric
 *
 */

int BrookLangevinDynamics::updateParameters( double temperature, double friction, double stepSize ){

   // ---------------------------------------------------------------------------------------

   static int showUpdate         = 1;
   static int maxShowUpdate      = 3;
   static const char* methodName = "\nBrookLangevinDynamics::updateParameters";

   // ---------------------------------------------------------------------------------------

   _setStepSize(    (BrookOpenMMFloat)  stepSize );
   _setFriction(    (BrookOpenMMFloat)  friction );
   //_setTau(    (BrookOpenMMFloat)  friction );
   _setTemperature( (BrookOpenMMFloat)  temperature );

   _updateDerivedParameters( );
   _updateSdStreams( );

   // show update

   if( showUpdate && getLog() && (showUpdate++ < maxShowUpdate) ){
      std::string contents = getContentsString( );
      (void) fprintf( getLog(), "%s contents\n%s", methodName, contents.c_str() );
      (void) fflush( getLog() );

   }

   return DefaultReturnValue;

};

/** 
 * Update
 * 
 * @param  positions                   atom positions
 * @param  velocities                  atom velocities
 * @param  forces                      atom forces
 * @param  brookShakeAlgorithm         BrookShakeAlgorithm reference
 * @param  brookRandomNumberGenerator  BrookRandomNumberGenerator reference
 *
 * @return  DefaultReturnValue
 *
 */

int BrookLangevinDynamics::update( Stream& positions, Stream& velocities,
                                     const Stream& forces,
                                     BrookShakeAlgorithm& brookShakeAlgorithm,
                                     BrookRandomNumberGenerator& brookRandomNumberGenerator ){

// ---------------------------------------------------------------------------------------

   // unused Shake parameter

   float omega                   = 1.0f;

   static const char* methodName = "\nBrookLangevinDynamics::update";

   static const int PrintOn      = 0;

// ---------------------------------------------------------------------------------------

   const BrookOpenMMFloat* derivedParameters           = getDerivedParameters();

   BrookStreamImpl& positionStream                     = dynamic_cast<BrookStreamImpl&>       (positions.getImpl());
   BrookStreamImpl& velocityStream                     = dynamic_cast<BrookStreamImpl&>       (velocities.getImpl());
   const BrookStreamImpl& forceStreamC                 = dynamic_cast<const BrookStreamImpl&> (forces.getImpl());
   BrookStreamImpl& forceStream                        = const_cast<BrookStreamImpl&> (forceStreamC);

   if( (1 || PrintOn) && getLog() ){

      static int showAux = 1;

      if( PrintOn ){
         (void) fprintf( getLog(), "%s shake=%d\n", methodName, brookShakeAlgorithm.getNumberOfConstraints() );
         (void) fflush( getLog() );
      }

      // show update
   
      if( showAux ){
         showAux = 0;

         std::string contents             = brookRandomNumberGenerator.getContentsString( );
         (void) fprintf( getLog(), "%s RNG contents\n%s", methodName, contents.c_str() );
         // brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( getLog() ); 

         contents             = brookShakeAlgorithm.getContentsString( );
         (void) fprintf( getLog(), "%s Shake contents\n%s", methodName, contents.c_str() );

         (void) fflush( getLog() );
      }
   }

   // first integration step

   kupdate_sd1_fix1(
         (float) getLangevinDynamicsAtomStreamWidth(),
         (float) brookRandomNumberGenerator.getRandomNumberStreamWidth(),
         (float) brookRandomNumberGenerator.getRvStreamOffset(),
         derivedParameters[EM],
         derivedParameters[Sd1pc1],
         derivedParameters[Sd1pc2],
         derivedParameters[Sd1pc3],
         getSDPC1Stream()->getBrookStream(),
         brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->getBrookStream(),
         getSD2XStream()->getBrookStream(),
         positionStream.getBrookStream(),
         forceStream.getBrookStream(),
         velocityStream.getBrookStream(),
         getInverseMassStream()->getBrookStream(),
         getSD1VStream()->getBrookStream(),
         getVPrimeStream()->getBrookStream(),
         getXPrimeStream()->getBrookStream() 
         );
 
   // diagnostics

   if( PrintOn ){
      (void) fprintf( getLog(), "\nPost kupdate_sd1_fix1: atomStrW=%3d rngStrW=%3d rngOff=%5d "
                                "EM=%12.5e Sd1pc[]=[%12.5e %12.5e %12.5e]",
                                getLangevinDynamicsAtomStreamWidth(),
                                brookRandomNumberGenerator.getRandomNumberStreamWidth(),
                                brookRandomNumberGenerator.getRvStreamOffset(),
                                derivedParameters[EM], derivedParameters[Sd1pc1], derivedParameters[Sd1pc2], derivedParameters[Sd1pc3] );

      (void) fprintf( getLog(), "\nSDPC1Stream\n" );
      getSDPC1Stream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nSD2XStream\n" );
      getSD2XStream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nInverseMassStream\n" );
      getInverseMassStream()->printToFile( getLog() );

      //StreamImpl& positionStreamImpl               = positionStream.getImpl();
      //const BrookStreamImpl brookPositions         = dynamic_cast<BrookStreamImpl&> (positionStreamImpl);
      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalPos->printToFile( getLog() );

      (void) fprintf( getLog(), "\nForceStream\n" );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamImpl();
      brookStreamInternalF->printToFile( getLog() );

      BrookStreamInternal* brookStreamInternalV   = velocityStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nVelocityStream\n" );
      brookStreamInternalV->printToFile( getLog() );

      (void) fprintf( getLog(), "\nSD1VStream\n" );
      getSD1VStream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nVPrimeStream\n" );
      getVPrimeStream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nXPrimeStream\n" );
      getXPrimeStream()->printToFile( getLog() ); 

      (void) fprintf( getLog(), "\nRvStreamIndex=%d\n", brookRandomNumberGenerator.getRvStreamIndex() );
      // brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( getLog() ); 
   }   

   // advance random number cursor

   brookRandomNumberGenerator.advanceGVCursor( 2*getNumberOfAtoms() );

   // first Shake step

   if( brookShakeAlgorithm.getNumberOfConstraints() > 0 ){
      kshakeh_fix1( 
                    10.0f,
                    (float) getLangevinDynamicsAtomStreamWidth(),
                    brookShakeAlgorithm.getInverseHydrogenMass(),
                    omega, 
                    brookShakeAlgorithm.getShakeAtomIndicesStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeAtomParameterStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream() );
   
      // first Shake gather
   
      kshakeh_update1_fix1(
                    (float) getLangevinDynamicsAtomStreamWidth(),
                    derivedParameters[Sd2pc1],
                    brookShakeAlgorithm.getShakeInverseMapStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    getVPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream(),
                    getXPrimeStream()->getBrookStream() );
   }

   // second integration step

   kupdate_sd2_fix1(
         (float) getLangevinDynamicsAtomStreamWidth(),
         (float) brookRandomNumberGenerator.getRandomNumberStreamWidth(),
         (float) brookRandomNumberGenerator.getRvStreamOffset(),
         derivedParameters[Sd2pc1],
         derivedParameters[Sd2pc2],
         getSDPC2Stream()->getBrookStream(),
         brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->getBrookStream(),
         getSD1VStream()->getBrookStream(),
         positionStream.getBrookStream(),
         getXPrimeStream()->getBrookStream(), 
         getVPrimeStream()->getBrookStream(),
         getSD2XStream()->getBrookStream(),
         velocityStream.getBrookStream(),
         getXPrimeStream()->getBrookStream() 
         );

   // diagnostics

   if( PrintOn ){
      (void) fprintf( getLog(), "\nPost kupdate_sd2_fix1: atomStrW=%3d rngStrW=%3d rngOff=%5d "
                                "Sd2pc[]=[%12.5e %12.5e]",
                                getLangevinDynamicsAtomStreamWidth(),
                                brookRandomNumberGenerator.getRandomNumberStreamWidth(),
                                brookRandomNumberGenerator.getRvStreamOffset(),
                                derivedParameters[Sd2pc1], derivedParameters[Sd2pc2] );

      (void) fprintf( getLog(), "\nSDPC2Stream\n" );
      getSDPC2Stream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nSD2XStream\n" );
      getSD1VStream()->printToFile( getLog() );

      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalPos->printToFile( getLog() );

      (void) fprintf( getLog(), "\nVPrimeStream\n" );
      getVPrimeStream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nXPrimeStream\n" );
      getXPrimeStream()->printToFile( getLog() ); 

      (void) fprintf( getLog(), "\ngetSD2XStream\n" );
      getSD2XStream()->printToFile( getLog() );

      BrookStreamInternal* brookStreamInternalVel  = velocityStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nVelocityStream\n" );
      brookStreamInternalVel->printToFile( getLog() );

      (void) fprintf( getLog(), "\nRvStreamIndex=%d\n", brookRandomNumberGenerator.getRvStreamIndex() );
//      brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( getLog() ); 
   }   

   // advance random number cursor

   brookRandomNumberGenerator.advanceGVCursor( 2*getNumberOfAtoms() );

   // second Shake step

   if( brookShakeAlgorithm.getNumberOfConstraints() > 0 ){
      kshakeh_fix1( 
                    10.0f,
                    (float) getLangevinDynamicsAtomStreamWidth(),
                    brookShakeAlgorithm.getInverseHydrogenMass(),
                    omega, 
                    brookShakeAlgorithm.getShakeAtomIndicesStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeAtomParameterStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream() );
   
      // second Shake gather
   
      kshakeh_update2_fix1( 
                    (float) getLangevinDynamicsAtomStreamWidth(),
                    brookShakeAlgorithm.getShakeInverseMapStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream(),
                    positionStream.getBrookStream() );
   
   } else {
      //kadd3( getXPrimeStream()->getBrookStream(), positionStream.getBrookStream() );
      ksetStr3( getXPrimeStream()->getBrookStream(), positionStream.getBrookStream() );
   }

   return DefaultReturnValue;

};

/**
 *
 * Get array of derived parameters indexed by 'DerivedParameters' enums
 *
 * @return array
 *
 */
   
const BrookOpenMMFloat* BrookLangevinDynamics::getDerivedParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nBrookLangevinDynamics::getDerivedParameters";

   // ---------------------------------------------------------------------------------------

   return _derivedParameters;
}

/** 
 * Get Atom stream size
 *
 * @return  Atom stream size
 *
 */

int BrookLangevinDynamics::getLangevinDynamicsAtomStreamSize( void ) const {
   return _sdAtomStreamSize;
}

/** 
 * Get atom stream width
 *
 * @return  atom stream width
 *
 */

int BrookLangevinDynamics::getLangevinDynamicsAtomStreamWidth( void ) const {
   return _sdAtomStreamWidth;
}

/** 
 * Get atom stream height
 *
 * @return atom stream height
 */

int BrookLangevinDynamics::getLangevinDynamicsAtomStreamHeight( void ) const {
   return _sdAtomStreamHeight;
}

/** 
 * Get SDPC1 stream 
 *
 * @return  SDPC1 stream
 *
 */

BrookFloatStreamInternal* BrookLangevinDynamics::getSDPC1Stream( void ) const {
   return _sdStreams[SDPC1Stream];
}

/** 
 * Get SDPC2 stream 
 *
 * @return  SDPC2 stream
 *
 */

BrookFloatStreamInternal* BrookLangevinDynamics::getSDPC2Stream( void ) const {
   return _sdStreams[SDPC2Stream];
}

/** 
 * Get SD2X stream 
 *
 * @return  SD2X stream
 *
 */

BrookFloatStreamInternal* BrookLangevinDynamics::getSD2XStream( void ) const {
   return _sdStreams[SD2XStream];
}

/** 
 * Get SD1V stream 
 *
 * @return  SD1V stream
 *
 */

BrookFloatStreamInternal* BrookLangevinDynamics::getSD1VStream( void ) const {
   return _sdStreams[SD1VStream];
}

/** 
 * Get VPrime stream 
 *
 * @return  Vprime stream
 *
 */

BrookFloatStreamInternal* BrookLangevinDynamics::getVPrimeStream( void ) const {
   return _sdStreams[VPrimeStream];
}

/** 
 * Get XPrime stream 
 *
 * @return  Xprime stream
 *
 */

BrookFloatStreamInternal* BrookLangevinDynamics::getXPrimeStream( void ) const {
   return _sdStreams[XPrimeStream];
}

/** 
 * Get InverseMass stream 
 *
 * @return  inverse mass stream
 *
 */

BrookFloatStreamInternal* BrookLangevinDynamics::getInverseMassStream( void ) const {
   return _sdStreams[InverseMassStream];
}

/** 
 * Initialize stream dimensions
 * 
 * @param numberOfAtoms             number of atoms
 * @param platform                  platform
 *      
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookLangevinDynamics::_initializeStreamSizes( int numberOfAtoms, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookLangevinDynamics::_initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _sdAtomStreamSize     = getAtomStreamSize( platform );
   _sdAtomStreamWidth    = getAtomStreamWidth( platform );
   _sdAtomStreamHeight   = getAtomStreamHeight( platform );

   return DefaultReturnValue;
}

/** 
 * Initialize stream dimensions
 * 
 * @param numberOfAtoms             number of atoms
 * @param platform                  platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

//std::string BrookLangevinDynamics::_getDerivedParametersString( BrookLangevinDynamics::DerivedParameters  derivedParametersIndex ) const {
std::string BrookLangevinDynamics::_getDerivedParametersString( int derivedParametersIndex ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookLangevinDynamics::_getDerivedParametersString";

// ---------------------------------------------------------------------------------------

   std::string returnString;
 
   switch( derivedParametersIndex ){
 
      case GDT:
         returnString = "GDT";
         break;
 
      case EPH:
         returnString = "EPH";
         break;
 
      case EMH:
         returnString = "EMH";
         break;
 
      case EP:
         returnString = "EP";
         break;
 
      case EM:
         returnString = "EM";
         break;
 
      case B:
         returnString = "B";
         break;
 
      case C:
         returnString = "C";
         break;
 
      case D:
         returnString = "D";
         break;
 
      case V:
         returnString = "V";
         break;
 
      case X:
         returnString = "X";
         break;
 
      case Yv:
         returnString = "Yv";
         break;
 
      case Yx:
         returnString = "Yx";
         break;
 
      case Sd1pc1:
         returnString = "Sd1pc1";
         break;
 
      case Sd1pc2:
         returnString = "Sd1pc2";
         break;
 
      case Sd1pc3:
         returnString = "Sd1pc3";
         break;
 
      case Sd2pc1:
         returnString = "Sd2pc1";
         break;
 
      case Sd2pc2:
         returnString = "Sd2pc2";
         break;
 
      default:
         returnString = "Unknown";
         break;
    }
 
    return returnString;
}

/** 
 * Initialize streams
 * 
 * @param platform                  platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookLangevinDynamics::_initializeStreams( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookLangevinDynamics::_initializeStreams";

   BrookOpenMMFloat dangleValue            = (BrookOpenMMFloat) 0.0;

// ---------------------------------------------------------------------------------------

   int sdAtomStreamSize             = getLangevinDynamicsAtomStreamSize();
   int sdAtomStreamWidth            = getLangevinDynamicsAtomStreamWidth();

    _sdStreams[SDPC1Stream]         = new BrookFloatStreamInternal( BrookCommon::SDPC1Stream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float2, dangleValue );

    _sdStreams[SDPC2Stream]         = new BrookFloatStreamInternal( BrookCommon::SDPC2Stream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float2, dangleValue );

    _sdStreams[SD2XStream]          = new BrookFloatStreamInternal( BrookCommon::SD2XStream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[SD1VStream]          = new BrookFloatStreamInternal( BrookCommon::SD1VStream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[VPrimeStream]        = new BrookFloatStreamInternal( BrookCommon::VPrimeStream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[XPrimeStream]        = new BrookFloatStreamInternal( BrookCommon::XPrimeStream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[InverseMassStream]   = new BrookFloatStreamInternal( BrookCommon::InverseMassStream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float, dangleValue );

   return DefaultReturnValue;
}

/** 
 * Update sd streams -- called after parameters change
 * 
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookLangevinDynamics::_updateSdStreams( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName       = "BrookLangevinDynamics::_updateSdStreams";

// ---------------------------------------------------------------------------------------

   int sdAtomStreamSize                      = getLangevinDynamicsAtomStreamSize();

   BrookOpenMMFloat* sdpc[2];
   for( int ii = 0; ii < 2; ii++ ){
      sdpc[ii] = new BrookOpenMMFloat[2*sdAtomStreamSize];
      memset( sdpc[ii], 0, 2*sdAtomStreamSize*sizeof( BrookOpenMMFloat ) ); 
   }
   BrookOpenMMFloat* inverseMass = new BrookOpenMMFloat[sdAtomStreamSize];
   memset( inverseMass, 0, sdAtomStreamSize*sizeof( BrookOpenMMFloat ) ); 

   const BrookOpenMMFloat* derivedParameters = getDerivedParameters( );
   int numberOfAtoms                         = getNumberOfAtoms();
   int index                                 = 0;
   for( int ii = 0; ii < numberOfAtoms; ii++, index += 2 ){

      sdpc[0][index]      = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[Yv]) );
      sdpc[0][index+1]    = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[V])  );

      sdpc[1][index]      = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[Yx]) );
      sdpc[1][index+1]    = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[X])  );

      inverseMass[ii]     = _inverseSqrtMasses[ii]*_inverseSqrtMasses[ii];

   }

   _sdStreams[SDPC1Stream]->loadFromArray( sdpc[0] );
   _sdStreams[SDPC2Stream]->loadFromArray( sdpc[1] );
   _sdStreams[InverseMassStream]->loadFromArray( inverseMass );

   for( int ii = 0; ii < 2; ii++ ){
      delete[] sdpc[ii];
   }
   delete[] inverseMass;

   // initialize SD2X

   BrookOpenMMFloat* sd2x = new BrookOpenMMFloat[3*sdAtomStreamSize];
   //SimTKOpenMMUtilities::setRandomNumberSeed( (uint32_t) getRandomNumberSeed() );

   memset( sd2x, 0, 3*sdAtomStreamSize*sizeof( BrookOpenMMFloat ) ); 

   index = 0;
   for( int ii = 0; ii < numberOfAtoms; ii++, index += 3 ){
      sd2x[index]        = _inverseSqrtMasses[ii]*derivedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
      sd2x[index+1]      = _inverseSqrtMasses[ii]*derivedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
      sd2x[index+2]      = _inverseSqrtMasses[ii]*derivedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
   }
   
   _sdStreams[SD2XStream]->loadFromArray( sd2x );

   delete[] sd2x;

   return DefaultReturnValue;

}

/** 
 * Set masses 
 * 
 * @param masses             atomic masses
 *
 */

int BrookLangevinDynamics::_setInverseSqrtMasses( const std::vector<double>& masses ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookLangevinDynamics::_setInverseSqrtMasses";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat one                     = (BrookOpenMMFloat)  1.0;

// ---------------------------------------------------------------------------------------

   // setup inverse sqrt masses

   _inverseSqrtMasses = new BrookOpenMMFloat[masses.size()];
   int index          = 0;
   for( std::vector<double>::const_iterator ii = masses.begin(); ii != masses.end(); ii++, index++ ){
      if( *ii != 0.0 ){
         BrookOpenMMFloat value    = static_cast<BrookOpenMMFloat>(*ii);
         _inverseSqrtMasses[index] = ( SQRT( one/value ) );
      } else {
         _inverseSqrtMasses[index] = zero;
      }
   }

   return DefaultReturnValue;
}   
 
/*  
 * Setup of LangevinDynamics parameters
 *
 * @param masses                masses
 * @param platform              Brook platform
 *
 * @return nonzero value if error
 *
 * */
    
int BrookLangevinDynamics::setup( const std::vector<double>& masses, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookLangevinDynamics::setup";

// ---------------------------------------------------------------------------------------

   const BrookPlatform brookPlatform            = dynamic_cast<const BrookPlatform&> (platform);
   setLog( brookPlatform.getLog() );

   int numberOfAtoms  = (int) masses.size();
   setNumberOfAtoms( numberOfAtoms );

   // set stream sizes and then create streams

   _initializeStreamSizes( numberOfAtoms, platform );
   _initializeStreams( platform );

   _setInverseSqrtMasses( masses );

   return DefaultReturnValue;
}

/* 
 * Get contents of object
 *
 * @param level   level of dump
 *
 * @return string containing contents
 *
 * */

std::string BrookLangevinDynamics::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookLangevinDynamics::getContentsString";

   static const unsigned int MAX_LINE_CHARS = 256;
   char value[MAX_LINE_CHARS];
   static const char* Set                   = "Set";
   static const char* NotSet                = "Not set";

// ---------------------------------------------------------------------------------------

   std::stringstream message;
   std::string tab   = "   ";

#ifdef WIN32
#define LOCAL_SPRINTF(a,b,c) sprintf_s( (a), MAX_LINE_CHARS, (b), (c) );   
#else
#define LOCAL_SPRINTF(a,b,c) sprintf( (a), (b), (c) );   
#endif

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfAtoms() );
   message << _getLine( tab, "Number of atoms:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamWidth() );
   message << _getLine( tab, "Atom stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamHeight() );
   message << _getLine( tab, "Atom stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamSize() );
   message << _getLine( tab, "Atom stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getTau() );
   message << _getLine( tab, "Tau:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getTemperature() );
   message << _getLine( tab, "Temperature:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getStepSize() );
   message << _getLine( tab, "Step size:", value ); 

   const BrookOpenMMFloat* derivedParameters = getDerivedParameters();
   for( int ii = 0; ii < MaxDerivedParameters; ii++ ){
      (void) LOCAL_SPRINTF( value, "%.5e", derivedParameters[ii] );
      message << _getLine( tab, _getDerivedParametersString( ii ), value ); 
   }

   (void) LOCAL_SPRINTF( value, "%.5f", getTemperature() );
   message << _getLine( tab, "Temperature:", value ); 

   message << _getLine( tab, "Log:",                  (getLog()                ? Set : NotSet) ); 

   message << _getLine( tab, "SDPC1:",                (getSDPC1Stream()        ? Set : NotSet) ); 
   message << _getLine( tab, "SDPC2:",                (getSDPC2Stream()        ? Set : NotSet) ); 
   message << _getLine( tab, "SD2X:",                 (getSD2XStream()         ? Set : NotSet) ); 
   message << _getLine( tab, "SD1V:",                 (getSD1VStream()         ? Set : NotSet) ); 
 
   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      message << std::endl;
      if( _sdStreams[ii] ){
         message << _sdStreams[ii]->getContentsString( );
      }
   }

#undef LOCAL_SPRINTF

   return message.str();
}
