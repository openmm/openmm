/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -i------------------------------------------------------------------------- *
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
#include "kernels/kshakeh.h"
#include "kernels/kupdatesd.h"
#include "kernels/kcommon.h"

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

   _numberOfParticles             = -1;
   _internalStepCount             = 0;

   // mark stream dimension variables as unset

   _sdParticleStreamWidth         = -1;
   _sdParticleStreamHeight        = -1;
   _sdParticleStreamSize          = -1;

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
   static const BrookOpenMMFloat zero = static_cast<BrookOpenMMFloat>( 0.0 ); 
   static const BrookOpenMMFloat one  = static_cast<BrookOpenMMFloat>( 1.0 ); 
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
   _tau   = static_cast<BrookOpenMMFloat>( (friction != 0.0) ? 1.0/friction : 0.0);
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
 * @throw OpenMMException if tau too small
 *
 */

int BrookLangevinDynamics::_updateDerivedParameters( void ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName      = "\nBrookLangevinDynamics::_updateDerivedParameters";

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

   if( fabsf( (float) tau ) < epsilon ){
      std::stringstream message;
      message << methodName << " tau=" << tau << " too small.";
      throw OpenMMException( message.str() );
   }   


   _derivedParameters[GDT]      = stepSize/tau;

   _derivedParameters[EPH]      = EXP(  half*_derivedParameters[GDT] );
   _derivedParameters[EMH]      = EXP( -half*_derivedParameters[GDT] );
   _derivedParameters[EM]       = EXP(      -_derivedParameters[GDT] );
   _derivedParameters[EP]       = EXP(       _derivedParameters[GDT] );

   if( _derivedParameters[GDT] >= static_cast<BrookOpenMMFloat>( 0.1 ) ){

      BrookOpenMMFloat term1    = _derivedParameters[EPH] - one;
                 term1         *= term1;
      _derivedParameters[B]     = _derivedParameters[GDT]*(_derivedParameters[EP] - one) - four*term1;

      _derivedParameters[C]     = _derivedParameters[GDT] - three + four*_derivedParameters[EMH] - _derivedParameters[EM];
      _derivedParameters[D]     = two - _derivedParameters[EPH] - _derivedParameters[EMH];

    } else {

      BrookOpenMMFloat term1        = half*_derivedParameters[GDT];
      BrookOpenMMFloat term2        = term1*term1;
      BrookOpenMMFloat term4        = term2*term2;

      BrookOpenMMFloat third        = static_cast<BrookOpenMMFloat>( ( 1.0/3.0 ) );
      BrookOpenMMFloat o7_9         = static_cast<BrookOpenMMFloat>( ( 7.0/9.0 ) );
      BrookOpenMMFloat o1_12        = static_cast<BrookOpenMMFloat>( ( 1.0/12.0 ) );
      BrookOpenMMFloat o17_90       = static_cast<BrookOpenMMFloat>( ( 17.0/90.0 ) );
      BrookOpenMMFloat o7_30        = static_cast<BrookOpenMMFloat>( ( 7.0/30.0 ) );
      BrookOpenMMFloat o31_1260     = static_cast<BrookOpenMMFloat>( ( 31.0/1260.0 ) );
      BrookOpenMMFloat o_360        = static_cast<BrookOpenMMFloat>( ( 1.0/360.0 ) );

      _derivedParameters[B]         = term4*( third  + term1*( third + term1*( o17_90 + term1*o7_9 )));
      _derivedParameters[C]         = term2*term1*( two*third + term1*( -half + term1*( o7_30 + term1*(-o1_12 + term1*o31_1260 ))));
      _derivedParameters[D]         = term2*( -one + term2*(-o1_12 - term2*o_360));
   }    

   BrookOpenMMFloat kT        = static_cast<BrookOpenMMFloat>( BOLTZ )*temperature;

   _derivedParameters[V]      = SQRT( kT*( one - _derivedParameters[EM]) );
   _derivedParameters[X]      = tau*SQRT( kT*_derivedParameters[C] );
   _derivedParameters[Yv]     = SQRT( kT*_derivedParameters[B]/_derivedParameters[C] );
   _derivedParameters[Yx]     = tau*SQRT( kT*_derivedParameters[B]/(one - _derivedParameters[EM]) );

   _derivedParameters[Sd1pc1] = tau*( one - _derivedParameters[EM] );
   _derivedParameters[Sd1pc2] = tau*( _derivedParameters[EPH] - _derivedParameters[EMH] );
   if(  fabsf( _derivedParameters[Sd1pc2] ) < 1.0e-06 ){
      _derivedParameters[Sd1pc2] = tau*_derivedParameters[GDT];
   }
   _derivedParameters[Sd1pc3] = _derivedParameters[D]/(tau*_derivedParameters[C] );

   // if tau was greater than 20000, then _derivedParameters[Sd1pc2] was zero
   // using tau*_derivedParameters[GDT] as approximation to tau*( _derivedParameters[EPH] - _derivedParameters[EMH] ) for values in this
   // range

   if( fabsf( _derivedParameters[Sd1pc2] ) > 1.0e-10 ){
      _derivedParameters[Sd2pc1] = one/_derivedParameters[Sd1pc2];
   } else {
      std::stringstream message;
      message << methodName << " Sd1pc2=" << _derivedParameters[Sd1pc2] << " is too small (Sd2pc1 is 1/Sd1pc2)";
      throw OpenMMException( message.str() );
   }
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

   static int showUpdate               = 1;
   static int maxShowUpdate            = 3;
   static const std::string methodName = "\nBrookLangevinDynamics::updateParameters";

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
      (void) fprintf( getLog(), "%s contents\n%s", methodName.c_str(), contents.c_str() );
      (void) fflush( getLog() );

   }

   return DefaultReturnValue;

}

/**
 *
 * Get array of derived parameters indexed by 'DerivedParameters' enums
 *
 * @return array
 *
 */
   
const BrookOpenMMFloat* BrookLangevinDynamics::getDerivedParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName  = "\nBrookLangevinDynamics::getDerivedParameters";

   // ---------------------------------------------------------------------------------------

   return _derivedParameters;
}

/** 
 * Get Particle stream size
 *
 * @return  Particle stream size
 *
 */

int BrookLangevinDynamics::getLangevinDynamicsParticleStreamSize( void ) const {
   return _sdParticleStreamSize;
}

/** 
 * Get particle stream width
 *
 * @return  particle stream width
 *
 */

int BrookLangevinDynamics::getLangevinDynamicsParticleStreamWidth( void ) const {
   return _sdParticleStreamWidth;
}

/** 
 * Get particle stream height
 *
 * @return particle stream height
 */

int BrookLangevinDynamics::getLangevinDynamicsParticleStreamHeight( void ) const {
   return _sdParticleStreamHeight;
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
 * @param numberOfParticles         number of particles
 * @param platform                  platform
 *      
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookLangevinDynamics::_initializeStreamSizes( int numberOfParticles, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookLangevinDynamics::_initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _sdParticleStreamSize     = getParticleStreamSize( platform );
   _sdParticleStreamWidth    = getParticleStreamWidth( platform );
   _sdParticleStreamHeight   = getParticleStreamHeight( platform );

   return DefaultReturnValue;
}

/** 
 * Initialize stream dimensions
 * 
 * @param numberOfParticles         number of particles
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

   int sdParticleStreamSize         = getLangevinDynamicsParticleStreamSize();
   int sdParticleStreamWidth        = getLangevinDynamicsParticleStreamWidth();

    _sdStreams[SDPC1Stream]         = new BrookFloatStreamInternal( BrookCommon::SDPC1Stream,
                                                                    sdParticleStreamSize, sdParticleStreamWidth,
                                                                    BrookStreamInternal::Float2, dangleValue );

    _sdStreams[SDPC2Stream]         = new BrookFloatStreamInternal( BrookCommon::SDPC2Stream,
                                                                    sdParticleStreamSize, sdParticleStreamWidth,
                                                                    BrookStreamInternal::Float2, dangleValue );

    _sdStreams[SD2XStream]          = new BrookFloatStreamInternal( BrookCommon::SD2XStream,
                                                                    sdParticleStreamSize, sdParticleStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[SD1VStream]          = new BrookFloatStreamInternal( BrookCommon::SD1VStream,
                                                                    sdParticleStreamSize, sdParticleStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[VPrimeStream]        = new BrookFloatStreamInternal( BrookCommon::VPrimeStream,
                                                                    sdParticleStreamSize, sdParticleStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[XPrimeStream]        = new BrookFloatStreamInternal( BrookCommon::XPrimeStream,
                                                                    sdParticleStreamSize, sdParticleStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[InverseMassStream]   = new BrookFloatStreamInternal( BrookCommon::InverseMassStream,
                                                                    sdParticleStreamSize, sdParticleStreamWidth,
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

   int sdParticleStreamSize                  = getLangevinDynamicsParticleStreamSize();

   // create and initialize sdpc streams

   BrookOpenMMFloat* sdpc[2];
   for( int ii = 0; ii < 2; ii++ ){
      sdpc[ii] = new BrookOpenMMFloat[2*sdParticleStreamSize];
      memset( sdpc[ii], 0, 2*sdParticleStreamSize*sizeof( BrookOpenMMFloat ) ); 
   }
   BrookOpenMMFloat* inverseMass = new BrookOpenMMFloat[sdParticleStreamSize];
   memset( inverseMass, 0, sdParticleStreamSize*sizeof( BrookOpenMMFloat ) ); 

   const BrookOpenMMFloat* derivedParameters = getDerivedParameters( );
   int numberOfParticles                     = getNumberOfParticles();
   int index                                 = 0;
   for( int ii = 0; ii < numberOfParticles; ii++, index += 2 ){

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

   BrookOpenMMFloat* sd2x = new BrookOpenMMFloat[3*sdParticleStreamSize];
   //SimTKOpenMMUtilities::setRandomNumberSeed( (uint32_t) getRandomNumberSeed() );

   memset( sd2x, 0, 3*sdParticleStreamSize*sizeof( BrookOpenMMFloat ) ); 

   index = 0;
   int useFixedRandomValue = 0;
   if( useFixedRandomValue ){

      // diagnostics only!

      BrookOpenMMFloat fixedRandomValue = static_cast<BrookOpenMMFloat>( 0.1 );
      for( int ii = 0; ii < numberOfParticles; ii++, index += 3 ){
         BrookOpenMMFloat value = _inverseSqrtMasses[ii]*derivedParameters[X]*fixedRandomValue;
         sd2x[index]            = value;
         sd2x[index+1]          = value;
         sd2x[index+2]          = value;
      }

      // print message letting user know non-random value being used

      FILE* log = getLog() ? getLog() : stderr;
      (void) fprintf( log, "%s using fixed 'random value'=%.3f to initialize sd2x\n", methodName.c_str(), fixedRandomValue );

   } else {

      for( int ii = 0; ii < numberOfParticles; ii++, index += 3 ){
         sd2x[index]        = _inverseSqrtMasses[ii]*derivedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
         sd2x[index+1]      = _inverseSqrtMasses[ii]*derivedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
         sd2x[index+2]      = _inverseSqrtMasses[ii]*derivedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
      }
   }
   
   _sdStreams[SD2XStream]->loadFromArray( sd2x );

   delete[] sd2x;

   return DefaultReturnValue;

}

/** 
 * Set masses 
 * 
 * @param masses             particle masses
 *
 */

int BrookLangevinDynamics::_setInverseSqrtMasses( const std::vector<double>& masses ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookLangevinDynamics::_setInverseSqrtMasses";

   BrookOpenMMFloat zero                    = static_cast<BrookOpenMMFloat>( 0.0 );
   BrookOpenMMFloat one                     = static_cast<BrookOpenMMFloat>( 1.0 );

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

   const BrookPlatform& brookPlatform       = dynamic_cast<const BrookPlatform&> (platform);
   setLog( brookPlatform.getLog() );

   int numberOfParticles  = (int) masses.size();
   setNumberOfParticles( numberOfParticles );

   // set stream sizes and then create streams

   _initializeStreamSizes( numberOfParticles, platform );
   _initializeStreams( platform );

   _setInverseSqrtMasses( masses );

   return DefaultReturnValue;
}

/** 
 * Get T
 * 
 * @param velocities             velocities
 * @param inverseMassStream      inverse masses
 * @param numberOfConstraints    number of constraints
 *
 * @return   temperature
 */

float BrookLangevinDynamics::getTemperature( BrookStreamInternal* velocities, BrookFloatStreamInternal* inverseMassStream,
                                             int numberOfConstraints ) const {
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookLangevinDynamics::getTemperature";

// ---------------------------------------------------------------------------------------

   void* dataArrayV           = velocities->getData( 1 );
   float* velocitiesI         = (float*) dataArrayV;

   void* inverseMassStreamV   = inverseMassStream->getData( 1 );
   float* inverseMassStreamI  = (float*) inverseMassStreamV;

   float ke                   = 0.0f;
   int index                  = 0;

   int numberOfParticles      = getNumberOfParticles();
   for( int ii = 0; ii < numberOfParticles; ii++, index += 3 ){
      ke    += (velocitiesI[index]*velocitiesI[index] + velocitiesI[index+1]*velocitiesI[index+1] + velocitiesI[index+2]*velocitiesI[index+2] )/inverseMassStreamI[ii];
   }

   int degreesOfFreedom = 3*getNumberOfParticles() - numberOfConstraints;
   float denominator    = 1.0f/( ((float) BOLTZ)*((float) ( degreesOfFreedom )) );
//(void) fprintf( stderr, "%s ke=%.5e T=%.3f dof=%d\n", methodName.c_str(), ke, (ke*denominator),  degreesOfFreedom ); 
   ke *= denominator;

   return ke;
}

/** 
 * Remove velocity com (diagnostics)
 * 
 * @param velocities             velocities
 * @param inverseMassStream      inverse masses
 *
 * @return DefaultReturnValue  
 */

int BrookLangevinDynamics::removeCom( BrookStreamInternal* velocities, BrookFloatStreamInternal* inverseMassStream ) const {
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookLangevinDynamics::removeCom";

// ---------------------------------------------------------------------------------------

   void* dataArrayV           = velocities->getData( 1 );
   float* velocitiesI         = (float*) dataArrayV;

   void* inverseMassStreamV   = inverseMassStream->getData( 1 );
   float* inverseMassStreamI  = (float*) inverseMassStreamV;

   float totalMass            = 0.0f;
   float com[3]               = { 0.0f, 0.0f, 0.0f };
   int index                  = 0;

   for( int ii = 0; ii < getNumberOfParticles(); ii++ ){
      float mass   = 1.0f/inverseMassStreamI[ii];
      totalMass   += mass;
      com[0]      += mass*velocitiesI[index];
      com[1]      += mass*velocitiesI[index+1];
      com[2]      += mass*velocitiesI[index+2];
      index       += 3;
   }
   totalMass   = 1.0f/totalMass;
   com[0]     *= totalMass;
   com[1]     *= totalMass;
   com[2]     *= totalMass;

   index                = 0;
   double* newVelocities = new double[velocities->getStreamSize()*velocities->getWidth()];
   memset( newVelocities, 0, sizeof( double )*velocities->getStreamSize()*velocities->getWidth() );
   for( int ii = 0; ii < getNumberOfParticles(); ii++ ){
      newVelocities[index]    = (double) velocitiesI[index]   - com[0];
      newVelocities[index+1]  = (double) velocitiesI[index+1] - com[1];
      newVelocities[index+2]  = (double) velocitiesI[index+2] - com[2];
      index                  += 3;
   }

   velocities->loadFromArray( newVelocities );

   dataArrayV          = velocities->getData( 1 );
   velocitiesI         = (float*) dataArrayV;

/*
   (void) fprintf( stderr, "%s readback\n", methodName.c_str() );
   for( int ii = 0; ii < velocities->getStreamSize()*3; ii += 3 ){
      (void) fprintf( stderr, "%s %d velocitiesI[%14.5e %14.5e %14.5e]\n", methodName.c_str(), ii/3, velocitiesI[ii], velocitiesI[ii+1], velocitiesI[ii+2] ); 
   }
*/

   velocities->loadFromArray( newVelocities );

   delete[] newVelocities;

   return DefaultReturnValue;
}

/** 
 * Reset velocities (diagnostics)
 * 
 * @param velocities             velocities
 *
 * @return DefaultReturnValue  
 */

int BrookLangevinDynamics::resetVelocities( BrookStreamInternal* velocities ) const {
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookLangevinDynamics::resetVelocities";

// ---------------------------------------------------------------------------------------

   // reset velocities to determinisitic values
   // note use of double instead of float for the load array

   double* newVelocities = new double[velocities->getStreamSize()*velocities->getWidth()];
   memset( newVelocities, 0, sizeof( double )*velocities->getStreamSize()*velocities->getWidth() );

   for( int ii = 1; ii <= 3*getNumberOfParticles(); ii++ ){
      int jj                = ii % 10;
      double sign           = jj % 2 ? 0.1 : -0.1;
      newVelocities[ii-1]   = sign*( (double) jj); 
   }

   velocities->loadFromArray( newVelocities );

   delete[] newVelocities;

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

int BrookLangevinDynamics::update( BrookStreamImpl& positionStream, BrookStreamImpl& velocityStream,
                                   BrookStreamImpl& forceStream,
                                   BrookShakeAlgorithm& brookShakeAlgorithm,
                                   BrookRandomNumberGenerator& brookRandomNumberGenerator ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nBrookLangevinDynamics::update";
   int printOn                         = 0;
   FILE* log;

// ---------------------------------------------------------------------------------------

   _internalStepCount++;

//setLog( stderr );

   if( printOn && getLog() ){
      log = getLog();
   } else { 
      printOn = 0;
   }

   const BrookOpenMMFloat* derivedParameters           = getDerivedParameters();

   if( printOn ){

      static int showAux = 1;

      if( printOn && showAux ){
         (void) fprintf( log, "%s step=%d shake=%d\n", methodName.c_str(), _internalStepCount, brookShakeAlgorithm.getNumberOfConstraints() );
         (void) fflush( log );
      }

      // show update
   
      if( showAux ){
         showAux = 0;

         std::string contents             = getContentsString( );
         (void) fprintf( log, "%s step=%d contents\n%s", methodName.c_str(), _internalStepCount, contents.c_str() );

         contents             = brookRandomNumberGenerator.getContentsString( );
         (void) fprintf( log, "%s step=%d RNG contents\n%s", methodName.c_str(), _internalStepCount, contents.c_str() );
         // brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( log ); 

         contents             = brookShakeAlgorithm.getContentsString( );
         (void) fprintf( log, "%s step=%d Shake contents\n%s", methodName.c_str(), _internalStepCount, contents.c_str() );
         (void) fflush( log );
      }

   }

   // diagnostics

   if( 0 && _internalStepCount == 1 ){ 
      //resetVelocities( velocityStream.getBrookStreamInternal() );
      std::string velocityFileName = "kupdate_sd1_strV.in";
      std::string rvName0          = "kupdate_sd1_0.fgauss.in";
      std::string rvName1          = "kupdate_sd1_1.fgauss.in";
      std::string sd2xName         = "kupdate_sd1_strSD2X.in";
      velocityStream.getBrookStreamInternal()->loadStreamGivenFileName( velocityFileName );
      brookRandomNumberGenerator.getRandomNumberStream( 0 )->loadStreamGivenFileName( rvName0 );
      brookRandomNumberGenerator.getRandomNumberStream( 1 )->loadStreamGivenFileName( rvName1 );
      getSD2XStream()->loadStreamGivenFileName( sd2xName );
   }

   // more diagnostics

   if( 0 && (_internalStepCount % 10) == 0 ){
      FILE*	log1 = stderr;
      (void) fprintf( log1, "\nVelocityStream %d XX\n", _internalStepCount ); fflush( log1 );
      BrookStreamInternal* brookStreamInternalPos  = velocityStream.getBrookStreamInternal();
      float temperature = getTemperature( brookStreamInternalPos, getInverseMassStream(), brookShakeAlgorithm.getNumberOfConstraints() );
//      removeCom( brookStreamInternalPos, getInverseMassStream() );
      float temperaturePost = getTemperature( brookStreamInternalPos, getInverseMassStream(), brookShakeAlgorithm.getNumberOfConstraints() );
      (void) fprintf( log1, "\nVelocityStream %d Tp=%.3f %.3f\n", _internalStepCount, temperature, temperaturePost );
      (void) fprintf( log1, "\n%s step=%d Post kupdate_sd1_fix1: particleStrW=%3d rngStrW=%3d rngOff=%5d "
                                "EM=%12.5e Sd1pc[]=[%12.5e %12.5e %12.5e]", methodName.c_str(), _internalStepCount,
                                getLangevinDynamicsParticleStreamWidth(),
                                brookRandomNumberGenerator.getRandomNumberStreamWidth(),
                                brookRandomNumberGenerator.getRvStreamOffset(),
                                derivedParameters[EM], derivedParameters[Sd1pc1], derivedParameters[Sd1pc2], derivedParameters[Sd1pc3] );
      (void) fprintf( log1, "\nRvStreamIndex=%d step=%d\n", brookRandomNumberGenerator.getRvStreamIndex(), _internalStepCount );
      fflush( log1 );
   }

   // first integration step

   kupdate_sd1_fix1(
         (float) getLangevinDynamicsParticleStreamWidth(),
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
         getXPrimeStream()->getBrookStream() );
 
   // diagnostics

   if( 0 && printOn ){
      (void) fprintf( log, "\n%s step=%d Post kupdate_sd1_fix1: particleStrW=%3d rngStrW=%3d rngOff=%5d "
                                "EM=%12.5e Sd1pc[]=[%12.5e %12.5e %12.5e]", methodName.c_str(), _internalStepCount,
                                getLangevinDynamicsParticleStreamWidth(),
                                brookRandomNumberGenerator.getRandomNumberStreamWidth(),
                                brookRandomNumberGenerator.getRvStreamOffset(),
                                derivedParameters[EM], derivedParameters[Sd1pc1], derivedParameters[Sd1pc2], derivedParameters[Sd1pc3] );

      if( _internalStepCount == 1 ){
         (void) fprintf( log, "\nSDPC1Stream fixed input sd1 step=%d\n", _internalStepCount );
         getSDPC1Stream()->printToFile( log );
   
         (void) fprintf( log, "\nSD2XStream fixed input sd1 step=%d\n", _internalStepCount );
         getSD2XStream()->printToFile( log );
   
         (void) fprintf( log, "\nInverseMassStream fixed input sd1 step=%d\n", _internalStepCount );
         getInverseMassStream()->printToFile( log );
   
      }

      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
      (void) fprintf( log, "\nPositionStream input sd1 step=%d\n", _internalStepCount );
      brookStreamInternalPos->printToFile( log );

      (void) fprintf( log, "\nForceStream input sd1 step=%d\n", _internalStepCount );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamInternal();
      brookStreamInternalF->printToFile( log );

      BrookStreamInternal* brookStreamInternalV   = velocityStream.getBrookStreamInternal();
      (void) fprintf( log, "\nVelocityStream input sd1 step=%d\n", _internalStepCount );
      brookStreamInternalV->printToFile( log );

      (void) fprintf( log, "\nSD1VStream output sd1 step=%d\n", _internalStepCount );
      getSD1VStream()->printToFile( log );

      (void) fprintf( log, "\nVPrimeStream output sd1 step=%d\n", _internalStepCount );
      getVPrimeStream()->printToFile( log );

      (void) fprintf( log, "\nXPrimeStream output sd1 step=%d\n", _internalStepCount );
      getXPrimeStream()->printToFile( log ); 

      if( _internalStepCount == 1 ){
            std::vector<BrookStreamInternal*> streams;
            streams.push_back( brookStreamInternalPos );
            streams.push_back( getXPrimeStream() );
            std::stringstream fileNameBaseS;
            fileNameBaseS << "Brook_Sd1PreShk_" << _internalStepCount << ".txt"; 
            BrookStreamInternal::printStreamsToFile( fileNameBaseS.str(), streams );
      }
      (void) fprintf( log, "\nRvStreamIndex=%d step=%d\n", brookRandomNumberGenerator.getRvStreamIndex(), _internalStepCount );
      // brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( log ); 
   }   

   // advance random number cursor

   brookRandomNumberGenerator.advanceGVCursor( 2*getNumberOfParticles() );

   // first Shake step

   if( brookShakeAlgorithm.getNumberOfConstraints() > 0 ){

      kshakeh_fix1( 
                    (float) brookShakeAlgorithm.getMaxIterations(),
                    (float) getLangevinDynamicsParticleStreamWidth(),
                    brookShakeAlgorithm.getShakeTolerance(),
                    brookShakeAlgorithm.getShakeParticleIndicesStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeParticleParameterStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream() );
   
      // first Shake gather
   
      kshakeh_update1_fix1(
                    (float) brookShakeAlgorithm.getShakeConstraintStreamWidth(),
                    brookShakeAlgorithm.getShakeInverseMapStream()->getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream(),
                    getXPrimeStream()->getBrookStream() );

      if( 0 && printOn ){

         (void) fprintf( log, "\n%s Post kshakeh_update1: sw=%d ShkCnstStrW=%3d tol=%.3f maxIt=%d",
                         methodName.c_str(), getLangevinDynamicsParticleStreamWidth(),
                         brookShakeAlgorithm.getShakeConstraintStreamWidth(),
                         brookShakeAlgorithm.getShakeTolerance(), brookShakeAlgorithm.getMaxIterations() );
   
         if( _internalStepCount == 1 ){
            (void) fprintf( log, "\nShakeInverseMapStream fixed input sd shake1 at step=%d\n" );
            brookShakeAlgorithm.getShakeInverseMapStream()->printToFile( log );
         }
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
         (void) fprintf( log, "\nPositionStream input sd shake1 at step=%d\n", _internalStepCount );
         brookStreamInternalPos->printToFile( log );
   
         (void) fprintf( log, "\nXPrimeStream output sd shake1 at step=%d\n", _internalStepCount );
         getXPrimeStream()->printToFile( log );  

         (void) fprintf( log, "\nShakeXCons0 output sd shake1 at step=%d\n", _internalStepCount );
         brookShakeAlgorithm.getShakeXCons0Stream()->printToFile( log );  

         (void) fprintf( log, "\nShakeXCons1 output sd shake1 at step=%d\n", _internalStepCount );
         brookShakeAlgorithm.getShakeXCons1Stream()->printToFile( log );  

         (void) fprintf( log, "\nShakeXCons2 output sd shake1 at step=%d\n", _internalStepCount );
         brookShakeAlgorithm.getShakeXCons2Stream()->printToFile( log );  

         (void) fprintf( log, "\nShakeXCons3 output sd shake1 at step=%d\n", _internalStepCount );
         brookShakeAlgorithm.getShakeXCons3Stream()->printToFile( log );  

         if( _internalStepCount < 2 ){
            std::vector<BrookStreamInternal*> streams;
            streams.push_back( brookStreamInternalPos );
            streams.push_back( getXPrimeStream() );
            std::stringstream fileNameBaseS;
            fileNameBaseS << "Brook_Sd1PostShk_" << _internalStepCount << ".txt"; 
            BrookStreamInternal::printStreamsToFile( fileNameBaseS.str(), streams );
         }
      }   

   }

   // second integration step

   kupdate_sd2_fix1(
         (float) getLangevinDynamicsParticleStreamWidth(),
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

   if( 0 && printOn ){
      (void) fprintf( log, "\n%s step=%d Post kupdate_sd2_fix1: particleStrW=%3d rngStrW=%3d rngOff=%5d "
                                "Sd2pc[]=[%12.5e %12.5e]", methodName.c_str(), _internalStepCount,
                                getLangevinDynamicsParticleStreamWidth(),
                                brookRandomNumberGenerator.getRandomNumberStreamWidth(),
                                brookRandomNumberGenerator.getRvStreamOffset(),
                                derivedParameters[Sd2pc1], derivedParameters[Sd2pc2] );

      (void) fprintf( log, "\nSDPC2Stream input sd2 step=%d\n", _internalStepCount );
      getSDPC2Stream()->printToFile( log );

      (void) fprintf( log, "\ngetSD1VStream input sd2 step=%d\n", _internalStepCount );
      getSD1VStream()->printToFile( log );

      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
      (void) fprintf( log, "\nPositionStream input sd2 step=%d \n", _internalStepCount );
      brookStreamInternalPos->printToFile( log );

      (void) fprintf( log, "\nVPrimeStream input sd2 step=%d\n", _internalStepCount );
      getVPrimeStream()->printToFile( log );

      (void) fprintf( log, "\nXPrimeStream output sd2 step=%d\n", _internalStepCount );
      getXPrimeStream()->printToFile( log ); 

      (void) fprintf( log, "\nSD2XStream output sd2 step=%d\n", _internalStepCount );
      getSD2XStream()->printToFile( log );

      BrookStreamInternal* brookStreamInternalVel = velocityStream.getBrookStreamInternal();
      (void) fprintf( log, "\nVelocityStream output sd2 step=%d\n", _internalStepCount );
      brookStreamInternalVel->printToFile( log );

      (void) fprintf( log, "\nRvStreamIndex=%d step=%d\n", brookRandomNumberGenerator.getRvStreamIndex(), _internalStepCount );
//      brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( log ); 
 
      if( _internalStepCount < 3 ){
         std::vector<BrookStreamInternal*> streams;
         streams.push_back( brookStreamInternalVel );
         streams.push_back( getXPrimeStream() );
         std::stringstream fileNameBaseS;
         fileNameBaseS << "Brook_Sd2Out_" << _internalStepCount << ".txt"; 
         BrookStreamInternal::printStreamsToFile( fileNameBaseS.str(), streams );
      }
   }   

   // advance random number cursor

   brookRandomNumberGenerator.advanceGVCursor( 2*getNumberOfParticles() );

   // second Shake step

   if( brookShakeAlgorithm.getNumberOfConstraints() > 0 ){

      kshakeh_fix1( 
                    (float) brookShakeAlgorithm.getMaxIterations(),
                    (float) getLangevinDynamicsParticleStreamWidth(),
                    brookShakeAlgorithm.getShakeTolerance(),
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
                    (float) brookShakeAlgorithm.getShakeConstraintStreamWidth(),
                    brookShakeAlgorithm.getShakeInverseMapStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream(),
                    positionStream.getBrookStream() );
   
      // diagnostics
   
      if( 0 && printOn ){
         (void) fprintf( log, "\n%s step=%d Post kshakeh_update2_fix1: ShakeConstraintStreamWidth=%3d rngStrW=%3d rngOff=%5d "
                              "Sd2pc[]=[%12.5e %12.5e]\n", methodName.c_str(), _internalStepCount,
                              brookShakeAlgorithm.getShakeConstraintStreamWidth(),
                              brookRandomNumberGenerator.getRandomNumberStreamWidth(),
                              brookRandomNumberGenerator.getRvStreamOffset(),
                              derivedParameters[Sd2pc1], derivedParameters[Sd2pc2] );
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
         std::string violationString;
         brookShakeAlgorithm.checkConstraints( brookStreamInternalPos, violationString, 0.0001f );
         (void) fprintf( log, "\nPositionStream output sd shake2 step=%d %s\n", _internalStepCount, violationString.c_str() );
         brookStreamInternalPos->printToFile( log );
   
         BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamInternal();
         std::vector<std::vector<double> > forceStatistics;
         brookStreamInternalF->getStatistics( forceStatistics, getNumberOfParticles() );

         std::stringstream tag;
         tag << _internalStepCount << " Fxx "; 
         std::string stats = brookStreamInternalF->printStatistics( tag.str(), forceStatistics );
         (void) fprintf( log, "\nStep %d Force stats:\n%s", _internalStepCount, stats.c_str() );
         brookStreamInternalF->printToFile( log );

         brookStreamInternalPos  = velocityStream.getBrookStreamInternal();
         std::vector<std::vector<double> > velocityStatistics;
         brookStreamInternalPos->getStatistics( velocityStatistics, getNumberOfParticles() );
         std::stringstream tagV;
         tagV << _internalStepCount << " Vxx "; 
         stats = brookStreamInternalPos->printStatistics( tagV.str(), velocityStatistics );
         (void) fprintf( log, "\nStep %d Velocity stats:\n%s", _internalStepCount, stats.c_str() );
         float temperature = getTemperature( brookStreamInternalPos, getInverseMassStream(), brookShakeAlgorithm.getNumberOfConstraints() );
         (void) fprintf( log, "\nVelocityStream %d T=%.3f\n", _internalStepCount, temperature );
         brookStreamInternalPos->printToFile( log );
   
         (void) fprintf( log, "\nXPrimeStream input sd shake2 step=%d\n", _internalStepCount );
         getXPrimeStream()->printToFile( log ); 
   
//      brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( log ); 
      }   

   } else {

      // no constraints

      if( 0 && printOn ){
         (void) fprintf( log, "\n%s Pre ksetStr3 (no constraints)", methodName.c_str() );
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
         (void) fprintf( log, "\nPositionStream\n" );
         brookStreamInternalPos->printToFile( log );
      }

      kadd3( getXPrimeStream()->getBrookStream(), positionStream.getBrookStream(), positionStream.getBrookStream() );

      // diagnostics
   
      if( printOn ){
         (void) fprintf( log, "\n%s Post ksetStr3 (no constraints)", methodName.c_str() );

         (void) fprintf( log, "\nXPrimeStream\n" );
         getXPrimeStream()->printToFile( log );
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
         (void) fprintf( log, "\nPositionStream\n" );
         brookStreamInternalPos->printToFile( log );
   
//      brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( log ); 
      }   
   }

   // diagnostics

   if( 0 && (_internalStepCount % 1000) == 0 ){
      FILE*	log1     = stderr;
      float  epsilon = 1.0e-01f;

      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
      BrookStreamInternal* brookStreamInternalVel  = velocityStream.getBrookStreamInternal();
      BrookStreamInternal* brookStreamInternalFrc  = forceStream.getBrookStreamInternal();

      // check for nan and infinities

      int coordinateNans                           = brookStreamInternalPos->checkForNans( );
      int velocityNans                             = brookStreamInternalVel->checkForNans( );
      int forceNans                                = brookStreamInternalFrc->checkForNans( );
      int abort                                    = abs( coordinateNans ) + abs( velocityNans ) + abs( forceNans );

      // Shake violations

      std::string violationString;
      int constraintViolations                     = brookShakeAlgorithm.checkConstraints( brookStreamInternalPos, violationString, 0.0001f );

      abort                                       += abs( constraintViolations );

      // check T consistent w/ specified value

      float temperature                            = getTemperature( brookStreamInternalVel, getInverseMassStream(), brookShakeAlgorithm.getNumberOfConstraints() );
      if( fabsf( temperature - getTemperature() ) > 2.0f*temperature ){
         abort++;
      }

      // force sums ~ 0?
      std::vector<float> sums;
      brookStreamInternalFrc->sumColumns( sums );

      // check if should abort

      (void) fprintf( log1, "%d T=%.3f Nans: x=%d v=%d f=%d ", _internalStepCount, temperature, coordinateNans, velocityNans, forceNans );
      (void) fprintf( log1, " Fsum[" );
      for( int ii = 0; ii < 3; ii++ ){
         if( fabsf( sums[ii] ) > epsilon ){
            abort++;
         }
         (void) fprintf( log1, "%12.4e ", sums[ii] );
      }
      (void) fprintf( log1, "] %s abort=%d\n", violationString.c_str(), abort );
      (void) fflush( log1 );
      if( abort ){

         int nans[2];
         nans[0] = brookRandomNumberGenerator.getRandomNumberStream( 0 )->checkForNans();
         nans[1] = brookRandomNumberGenerator.getRandomNumberStream( 1 )->checkForNans();
         (void) fprintf( log1, "Aborting: Nans rng: active index=%d %d %d\n", brookRandomNumberGenerator.getRvStreamIndex(), nans[0], nans[1] );

         brookStreamInternalPos->printToFile( log1 );
         brookStreamInternalVel->printToFile( log1 );
         brookStreamInternalFrc->printToFile( log1 );

         brookRandomNumberGenerator.getRandomNumberStream( brookRandomNumberGenerator.getRvStreamIndex() )->printToFile( log1 );

         exit(1);
      }
  
/*
      std::vector<std::vector<double> > velocityStatistics;
      brookStreamInternalPos->getStatistics( velocityStatistics, getNumberOfParticles() );
      std::stringstream tagV;
      tagV << _internalStepCount << " Vxx "; 
      std::string stats = brookStreamInternalPos->printStatistics( tagV.str(), velocityStatistics );
      (void) fprintf( log1, "\nStep %d Velocity stats:\n%s", _internalStepCount, stats.c_str() );
*/
      //removeCom( brookStreamInternalPos, getInverseMassStream() );
   }

   return DefaultReturnValue;

}
