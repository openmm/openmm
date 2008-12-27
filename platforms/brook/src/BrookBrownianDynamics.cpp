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
#include "BrookBrownianDynamics.h"
#include "BrookPlatform.h"
#include "OpenMMException.h"
#include "BrookStreamImpl.h"
#include "gpu/kshakeh.h"
#include "gpu/kupdatebd.h"
#include "gpu/kcommon.h"

// use random number generator

#include "../../reference/src/SimTKUtilities/SimTKOpenMMUtilities.h"
#include "../../reference/src/SimTKUtilities/SimTKOpenMMRealType.h"

using namespace OpenMM;
using namespace std;

/** 
 *
 * Constructor
 * 
 */

BrookBrownianDynamics::BrookBrownianDynamics( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookBrownianDynamics::BrookBrownianDynamics";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat one                     = (BrookOpenMMFloat)  1.0;
   BrookOpenMMFloat oneMinus                = (BrookOpenMMFloat) -1.0;

// ---------------------------------------------------------------------------------------

   _numberOfParticles             = -1;

   // mark stream dimension variables as unset

   _bdParticleStreamWidth         = -1;
   _bdParticleStreamHeight        = -1;
   _bdParticleStreamSize          = -1;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _streams[ii]   = NULL;
   }

   _temperature      = oneMinus;
   _stepSize         = oneMinus;
   _tau              = oneMinus;
   _noiseAmplitude   = zero;
   _forceScale       = zero;

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

BrookBrownianDynamics::~BrookBrownianDynamics( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookBrownianDynamics::~BrookBrownianDynamics";

// ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      delete _streams[ii];
   }

   delete[] _inverseSqrtMasses;

}

/** 
 * Get tau
 * 
 * @return  tau
 *
 */

BrookOpenMMFloat BrookBrownianDynamics::getTau( void ) const {
   return _tau;
}

/** 
 * Get friction
 * 
 * @return  friction
 *
 */

BrookOpenMMFloat BrookBrownianDynamics::getFriction( void ) const {
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

BrookOpenMMFloat BrookBrownianDynamics::getTemperature( void ) const {
   return _temperature;
}

/** 
 * Get stepSize
 * 
 * @return  stepSize
 *
 */

BrookOpenMMFloat BrookBrownianDynamics::getStepSize( void ) const {
   return _stepSize;
}

/** 
 * Get force scale
 * 
 * @return  force scale
 *
 */

BrookOpenMMFloat BrookBrownianDynamics::getForceScale( void ) const {
   return _forceScale;
}

/** 
 * Get noise amplitude
 * 
 * @return  noise amplitude
 *
 */

BrookOpenMMFloat BrookBrownianDynamics::getNoiseAmplitude( void ) const {
   return _noiseAmplitude;
}

/** 
 * Set tau
 * 
 * @param tau   new tau value
 *
 * @return      DefaultReturnValue
 *
 */

int BrookBrownianDynamics::_setTau( BrookOpenMMFloat tau ){
   _tau = tau;
   _setNoiseAmplitude();
   _setForceScale();
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

int BrookBrownianDynamics::_setFriction( BrookOpenMMFloat friction ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookBrownianDynamics::_setFriction";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat one                     = (BrookOpenMMFloat)  1.0;

// ---------------------------------------------------------------------------------------

   _tau   = (BrookOpenMMFloat) ( (friction != zero) ? one/friction : zero);

   _setNoiseAmplitude();
   _setForceScale();

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

int BrookBrownianDynamics::_setTemperature( BrookOpenMMFloat temperature ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookBrownianDynamics::_setTemperature";

// ---------------------------------------------------------------------------------------

   _temperature = temperature;
   _setNoiseAmplitude();
   return DefaultReturnValue;
}

/** 
 * Set stepSize
 * 
 * @param   stepSize
 *
 * @return      DefaultReturnValue
 *
 * @throw OpenMMException, if stepSize <= 0
 *
 */

int BrookBrownianDynamics::_setStepSize( BrookOpenMMFloat stepSize ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookBrownianDynamics::_setStepSize";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat two                     = (BrookOpenMMFloat)  2.0;

// ---------------------------------------------------------------------------------------

   if( stepSize <= zero ){
      std::stringstream message;
      message << methodName << " step size=" << stepSize << " is invalid.";
      throw OpenMMException( message.str() );
   }    

   _stepSize = stepSize;
   _setNoiseAmplitude();
   _setForceScale();

   return DefaultReturnValue;
}

/** 
 * Set noise amplitude
 * 
 * @return  DefaultReturnValue
 *
 */

int BrookBrownianDynamics::_setNoiseAmplitude( void ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookBrownianDynamics::_setNoiseAmplitude";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat two                     = (BrookOpenMMFloat)  2.0;

// ---------------------------------------------------------------------------------------

   _noiseAmplitude =  _temperature*_stepSize*_tau;

   // skip sqrt if value is negative -- presumably not all values updated

   if( _noiseAmplitude > zero ){
      RealOpenMM intermediate = (RealOpenMM) BOLTZ*( (RealOpenMM) (two*_noiseAmplitude) );
      _noiseAmplitude         = (BrookOpenMMFloat) SQRT( intermediate );
   }
   return DefaultReturnValue;
}

/** 
 * Set force scale
 * 
 * @return  DefaultReturnValue
 *
 */

int BrookBrownianDynamics::_setForceScale( void ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookBrownianDynamics::_setForceScale";

// ---------------------------------------------------------------------------------------

   _forceScale = _stepSize*_tau;
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

int BrookBrownianDynamics::updateParameters( double temperature, double friction, double stepSize ){

   // ---------------------------------------------------------------------------------------

   static int showUpdate         = 1;
   static int maxShowUpdate      = 3;
   static const char* methodName = "\nBrookBrownianDynamics::updateParameters";

   // ---------------------------------------------------------------------------------------

   _setStepSize(    (BrookOpenMMFloat)  stepSize );
   _setFriction(    (BrookOpenMMFloat)  friction );
   //_setTau(    (BrookOpenMMFloat)  friction );
   _setTemperature( (BrookOpenMMFloat)  temperature );

   _updateStreams( );

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
 * @param  positions                   particle positions
 * @param  velocities                  particle velocities
 * @param  forces                      particle forces
 * @param  brookShakeAlgorithm         BrookShakeAlgorithm reference
 * @param  brookRandomNumberGenerator  BrookRandomNumberGenerator reference
 *
 * @return  DefaultReturnValue
 *
 */

int BrookBrownianDynamics::update( Stream& positions, Stream& velocities,
                                   const Stream& forces,
                                   BrookShakeAlgorithm& brookShakeAlgorithm,
                                   BrookRandomNumberGenerator& brookRandomNumberGenerator ){

// ---------------------------------------------------------------------------------------

   // unused Shake parameter

   float omega                   = 1.0f;
   float one                     = 1.0f;

   static const char* methodName = "\nBrookBrownianDynamics::update";

   static const int PrintOn      = 0;

// ---------------------------------------------------------------------------------------

   BrookStreamImpl& positionStream                     = dynamic_cast<BrookStreamImpl&>       (positions.getImpl());
   BrookStreamImpl& velocityStream                     = dynamic_cast<BrookStreamImpl&>       (velocities.getImpl());
   const BrookStreamImpl& forceStreamC                 = dynamic_cast<const BrookStreamImpl&> (forces.getImpl());
   BrookStreamImpl& forceStream                        = const_cast<BrookStreamImpl&> (forceStreamC);

   if( (1 || PrintOn) && getLog() ){

      static int showAux = 1;

      if( showAux ){
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

   // diagnostics

   if( (1 || PrintOn) ){
      (void) fprintf( getLog(), "\nPre kintegrate_bd: %d rngStrW=%3d rngOff=%5d "
                                "ForceScale=%12.5e NoiseAmplitude=%12.5e\n",
                                getBrownianDynamicsParticleStreamWidth(),
                                brookRandomNumberGenerator.getRandomNumberStreamWidth(),
                                brookRandomNumberGenerator.getRvStreamOffset(),
                                getForceScale(),  getNoiseAmplitude() );

      // (void) fprintf( getLog(), "\nInverseMassStream\n" );
      //getInverseMassStream()->printToFile( getLog() );

      //StreamImpl& positionStreamImpl               = positionStream.getImpl();
      //const BrookStreamImpl brookPositions         = dynamic_cast<BrookStreamImpl&> (positionStreamImpl);
/*
      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalPos->printToFile( getLog() );
*/

      double forceSum[3];
      BrookStreamInternal* brookStreamInternalFF       = forceStream.getBrookStreamImpl();
      BrookFloatStreamInternal* brookStreamInternalF   = dynamic_cast<BrookFloatStreamInternal*> (brookStreamInternalFF);
      brookStreamInternalF->sumByDimension( getNumberOfParticles(), forceSum );
      (void) fprintf( getLog(), "\nForceStream [%18.10e %18.10e %18.10e]\n", forceSum[0], forceSum[1], forceSum[2] );
      brookStreamInternalF->printToFile( getLog() );
/*

      (void) fprintf( getLog(), "\nXPrimeStream\n" );
      getXPrimeStream()->printToFile( getLog() ); 

      (void) fprintf( getLog(), "\nRvStreamIndex=%d\n", brookRandomNumberGenerator.getRvStreamIndex() );
*/
      // brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( getLog() ); 

         BrookStreamInternal* brookStreamInternalVel  = velocityStream.getBrookStreamImpl();
         (void) fprintf( getLog(), "\nVelocityStream\n" );
         brookStreamInternalVel->printToFile( getLog() ); 
   }   

   // integration step -- deltas returned in XPrime

   kintegrate_bd(
         (float) getBrownianDynamicsParticleStreamWidth(),
         (float) brookRandomNumberGenerator.getRandomNumberStreamWidth(),
         (float) brookRandomNumberGenerator.getRvStreamOffset(),
         getForceScale(), getNoiseAmplitude(),
         brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->getBrookStream(),
         positionStream.getBrookStream(), forceStream.getBrookStream(),
         getXPrimeStream()->getBrookStream() );
 
   // diagnostics

   if( PrintOn ){
      (void) fprintf( getLog(), "\nPost kintegrate_bd: %d rngStrW=%3d rngOff=%5d "
                                "ForceScale=%12.5e NoiseAmplitude=%12.5e\n",
                                getBrownianDynamicsParticleStreamWidth(),
                                brookRandomNumberGenerator.getRandomNumberStreamWidth(),
                                brookRandomNumberGenerator.getRvStreamOffset(),
                                getForceScale(),  getNoiseAmplitude() );

      // (void) fprintf( getLog(), "\nInverseMassStream\n" );
      //getInverseMassStream()->printToFile( getLog() );

      //StreamImpl& positionStreamImpl               = positionStream.getImpl();
      //const BrookStreamImpl brookPositions         = dynamic_cast<BrookStreamImpl&> (positionStreamImpl);
      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalPos->printToFile( getLog() );

      (void) fprintf( getLog(), "\nForceStream\n" );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamImpl();
      brookStreamInternalF->printToFile( getLog() );

      (void) fprintf( getLog(), "\nXPrimeStream\n" );
      getXPrimeStream()->printToFile( getLog() ); 

      (void) fprintf( getLog(), "\nRvStreamIndex=%d\n", brookRandomNumberGenerator.getRvStreamIndex() );
      // brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( getLog() ); 
   }   

   // advance random number cursor

   brookRandomNumberGenerator.advanceGVCursor( getNumberOfParticles() );

   // Shake

   if( brookShakeAlgorithm.getNumberOfConstraints() > 0 ){
      kshakeh_fix1( 
                    25.0f,
                    (float) getBrownianDynamicsParticleStreamWidth(),
                    brookShakeAlgorithm.getInverseHydrogenMass(),
                    brookShakeAlgorithm.getShakeParticleIndicesStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeParticleParameterStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream() );
   
      // second Shake gather
   
      kshakeh_update2_fix1( 
                    (float) getBrownianDynamicsParticleStreamWidth(),
                    brookShakeAlgorithm.getShakeInverseMapStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream(),
                    positionStream.getBrookStream() );
   

      // diagnostics
   
      if( PrintOn ){
   
         (void) fprintf( getLog(), "\nPost kupdate_bd:\n" );
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
         (void) fprintf( getLog(), "\nPositionStream\n" );
         brookStreamInternalPos->printToFile( getLog() );
   
         (void) fprintf( getLog(), "\nXPrimeStream\n" );
         getXPrimeStream()->printToFile( getLog() ); 
   
         BrookStreamInternal* brookStreamInternalVel  = velocityStream.getBrookStreamImpl();
         (void) fprintf( getLog(), "\nVelocityStream\n" );
         brookStreamInternalVel->printToFile( getLog() ); 
   
      }   
   
   } else {

      // update step
   
      float velocityScale = (float) (one/getStepSize());
      kupdate_bd2( velocityScale, getXPrimeStream()->getBrookStream(),
                   positionStream.getBrookStream(), velocityStream.getBrookStream(), positionStream.getBrookStream() );

      // diagnostics
   
      if( (1 || PrintOn) ){
   
         (void) fprintf( getLog(), "\nPost kupdate_bd2: velocityScale=%12.5e\n", velocityScale );

         (void) fprintf( getLog(), "\nXPrimeStream\n" );
         getXPrimeStream()->printToFile( getLog() ); 

         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
         (void) fprintf( getLog(), "\nPositionStream\n" );
         brookStreamInternalPos->printToFile( getLog() );
   
         BrookStreamInternal* brookStreamInternalVel  = velocityStream.getBrookStreamImpl();
         (void) fprintf( getLog(), "\nVelocityStream\n" );
         brookStreamInternalVel->printToFile( getLog() ); 
   
      }   
   
   }
 
   return DefaultReturnValue;

};

/** 
 * Get Particle stream size
 *
 * @return  Particle stream size
 *
 */

int BrookBrownianDynamics::getBrownianDynamicsParticleStreamSize( void ) const {
   return _bdParticleStreamSize;
}

/** 
 * Get particle stream width
 *
 * @return  particle stream width
 *
 */

int BrookBrownianDynamics::getBrownianDynamicsParticleStreamWidth( void ) const {
   return _bdParticleStreamWidth;
}

/** 
 * Get particle stream height
 *
 * @return particle stream height
 */

int BrookBrownianDynamics::getBrownianDynamicsParticleStreamHeight( void ) const {
   return _bdParticleStreamHeight;
}

/** 
 * Get XPrime stream 
 *
 * @return  Xprime stream
 *
 */

BrookFloatStreamInternal* BrookBrownianDynamics::getXPrimeStream( void ) const {
   return _streams[XPrimeStream];
}

/** 
 * Get InverseMass stream 
 *
 * @return  inverse mass stream
 *
 */

BrookFloatStreamInternal* BrookBrownianDynamics::getInverseMassStream( void ) const {
   return _streams[InverseMassStream];
}

/** 
 * Initialize stream dimensions
 * 
 * @param numberOfParticles             number of particles
 * @param platform                  platform
 *      
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookBrownianDynamics::_initializeStreamSizes( int numberOfParticles, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookBrownianDynamics::_initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _bdParticleStreamSize     = getParticleStreamSize( platform );
   _bdParticleStreamWidth    = getParticleStreamWidth( platform );
   _bdParticleStreamHeight   = getParticleStreamHeight( platform );

   return DefaultReturnValue;
}

/** 
 * Initialize streams
 * 
 * @param platform                  platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookBrownianDynamics::_initializeStreams( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookBrownianDynamics::_initializeStreams";

   BrookOpenMMFloat dangleValue            = (BrookOpenMMFloat) 0.0;

// ---------------------------------------------------------------------------------------

   int sdParticleStreamSize             = getBrownianDynamicsParticleStreamSize();
   int sdParticleStreamWidth            = getBrownianDynamicsParticleStreamWidth();

    _streams[VPrimeStream]          = new BrookFloatStreamInternal( BrookCommon::VPrimeStream,
                                                                    sdParticleStreamSize, sdParticleStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _streams[XPrimeStream]          = new BrookFloatStreamInternal( BrookCommon::XPrimeStream,
                                                                    sdParticleStreamSize, sdParticleStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _streams[InverseMassStream]     = new BrookFloatStreamInternal( BrookCommon::InverseMassStream,
                                                                    sdParticleStreamSize, sdParticleStreamWidth,
                                                                    BrookStreamInternal::Float, dangleValue );

   return DefaultReturnValue;
}

/** 
 * Update streams -- called after parameters change
 * 
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookBrownianDynamics::_updateStreams( void ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName       = "BrookBrownianDynamics::_updateStreams";

// ---------------------------------------------------------------------------------------

   //int particleStreamSize                        = getBrownianDynamicsParticleStreamSize();

   return DefaultReturnValue;

}

/** 
 * Set masses 
 * 
 * @param masses             particle masses
 *
 */

int BrookBrownianDynamics::_setInverseSqrtMasses( const std::vector<double>& masses ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookBrownianDynamics::_setInverseSqrtMasses";

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
 * Setup of BrownianDynamics parameters
 *
 * @param masses                masses
 * @param platform              Brook platform
 *
 * @return nonzero value if error
 *
 * */
    
int BrookBrownianDynamics::setup( const std::vector<double>& masses, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookBrownianDynamics::setup";

// ---------------------------------------------------------------------------------------

   const BrookPlatform brookPlatform            = dynamic_cast<const BrookPlatform&> (platform);
   setLog( brookPlatform.getLog() );

   int numberOfParticles  = (int) masses.size();
   setNumberOfParticles( numberOfParticles );

   // set stream sizes and then create streams

   _initializeStreamSizes( numberOfParticles, platform );
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

std::string BrookBrownianDynamics::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookBrownianDynamics::getContentsString";

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

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfParticles() );
   message << _getLine( tab, "Number of particles:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamWidth() );
   message << _getLine( tab, "Particle stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamHeight() );
   message << _getLine( tab, "Particle stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamSize() );
   message << _getLine( tab, "Particle stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getTau() );
   message << _getLine( tab, "Tau:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getTemperature() );
   message << _getLine( tab, "Temperature:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getStepSize() );
   message << _getLine( tab, "Step size:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5e", getForceScale() );
   message << _getLine( tab, "Force scale:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5e", getNoiseAmplitude() );
   message << _getLine( tab, "Noise amplitude:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getTemperature() );
   message << _getLine( tab, "Temperature:", value ); 

   message << _getLine( tab, "Log:",                  (getLog()                ? Set : NotSet) ); 

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      message << std::endl;
      if( _streams[ii] ){
         message << _streams[ii]->getContentsString( );
      }
   }

#undef LOCAL_SPRINTF

   return message.str();
}
