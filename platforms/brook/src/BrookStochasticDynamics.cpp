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
#include "BrookStochasticDynamics.h"
#include "BrookPlatform.h"
#include "OpenMMException.h"
#include "BrookStreamImpl.h"

// use random number generator

#include "SimTKOpenMMUtilities.h"

using namespace OpenMM;
using namespace std;

/** 
 *
 * Constructor
 * 
 */

BrookStochasticDynamics::BrookStochasticDynamics( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStochasticDynamics::BrookStochasticDynamics";

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

BrookStochasticDynamics::~BrookStochasticDynamics( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStochasticDynamics::~BrookStochasticDynamics";

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

BrookOpenMMFloat BrookStochasticDynamics::getTau( void ) const {
   return _tau;
}

/** 
 * Get friction
 * 
 * @return  friction
 *
 */

BrookOpenMMFloat BrookStochasticDynamics::getFriction( void ) const {
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

BrookOpenMMFloat BrookStochasticDynamics::getTemperature( void ) const {
   return _temperature;
}

/** 
 * Get stepSize
 * 
 * @return  stepSize
 *
 */

BrookOpenMMFloat BrookStochasticDynamics::getStepSize( void ) const {
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

int BrookStochasticDynamics::_setTau( BrookOpenMMFloat tau ){
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

int BrookStochasticDynamics::_setFriction( BrookOpenMMFloat friction ){
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

int BrookStochasticDynamics::_setTemperature( BrookOpenMMFloat temperature ){
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

int BrookStochasticDynamics::_setStepSize( BrookOpenMMFloat stepSize ){
   _stepSize = stepSize;
   return DefaultReturnValue;
}

/** 
 * Update derived parameters
 * 
 * @return  DefaultReturnValue
 *
 */

int BrookStochasticDynamics::_updateDerivedParameters( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nBrookStochasticDynamics::_updateDerivedParameters";

   static const BrookOpenMMFloat zero       =  0.0;
   static const BrookOpenMMFloat one        =  1.0;
   static const BrookOpenMMFloat two        =  2.0;
   static const BrookOpenMMFloat three      =  3.0;
   static const BrookOpenMMFloat four       =  4.0;
   static const BrookOpenMMFloat half       =  0.5;

   // ---------------------------------------------------------------------------------------

   BrookOpenMMFloat tau         = getTau();
   BrookOpenMMFloat temperature = getTemperature();
   BrookOpenMMFloat stepSize    = getStepSize();

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

};

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

int BrookStochasticDynamics::updateParameters( double temperature, double friction, double stepSize ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nBrookStochasticDynamics::updateParameters";

   // ---------------------------------------------------------------------------------------

   _setStepSize(    (BrookOpenMMFloat)  stepSize );
   _setFriction(    (BrookOpenMMFloat)  friction );
   _setTemperature( (BrookOpenMMFloat)  temperature );

   _updateDerivedParameters( );
   _updateSdStreams( );

   return DefaultReturnValue;

};

/** 
 * Update
 * 
 * @param  positions           atom positions
 * @param  velocities          atom velocities
 * @param  forces              atom forces
 * @param  brookShakeAlgorithm BrookShakeAlgorithm reference
 *
 * @return  DefaultReturnValue
 *
 */

int BrookStochasticDynamics::update( Stream& positions, Stream& velocities,
                                     const Stream& forces,
                                     BrookShakeAlgorithm& brookShakeAlgorithm ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nBrookStochasticDynamics::update";

   // ---------------------------------------------------------------------------------------


   return DefaultReturnValue;

};

/**
 *
 * Get array of derived parameters indexed by 'DerivedParameters' enums
 *
 * @return array
 *
 */
   
const BrookOpenMMFloat* BrookStochasticDynamics::getDerivedParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nBrookStochasticDynamics::getDerivedParameters";

   // ---------------------------------------------------------------------------------------

   return _derivedParameters;
}

/** 
 * Get Atom stream size
 *
 * @return  Atom stream size
 *
 */

int BrookStochasticDynamics::getStochasticDynamicsAtomStreamSize( void ) const {
   return _sdAtomStreamSize;
}

/** 
 * Get atom stream width
 *
 * @return  atom stream width
 *
 */

int BrookStochasticDynamics::getStochasticDynamicsAtomStreamWidth( void ) const {
   return _sdAtomStreamWidth;
}

/** 
 * Get atom stream height
 *
 * @return atom stream height
 */

int BrookStochasticDynamics::getStochasticDynamicsAtomStreamHeight( void ) const {
   return _sdAtomStreamHeight;
}

/** 
 * Get SDPC1 stream 
 *
 * @return  SDPC1 stream
 *
 */

BrookFloatStreamInternal* BrookStochasticDynamics::getSDPC1Stream( void ) const {
   return _sdStreams[SDPC1Stream];
}

/** 
 * Get SDPC2 stream 
 *
 * @return  SDPC2 stream
 *
 */

BrookFloatStreamInternal* BrookStochasticDynamics::getSDPC2Stream( void ) const {
   return _sdStreams[SDPC2Stream];
}

/** 
 * Get SD2X stream 
 *
 * @return  SD2X stream
 *
 */

BrookFloatStreamInternal* BrookStochasticDynamics::getSD2XStream( void ) const {
   return _sdStreams[SD2XStream];
}

/** 
 * Get SD1V stream 
 *
 * @return  SD1V stream
 *
 */

BrookFloatStreamInternal* BrookStochasticDynamics::getSD1VStream( void ) const {
   return _sdStreams[SD1VStream];
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

int BrookStochasticDynamics::_initializeStreamSizes( int numberOfAtoms, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStochasticDynamics::_initializeStreamSizes";

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

//std::string BrookStochasticDynamics::_getDerivedParametersString( BrookStochasticDynamics::DerivedParameters  derivedParametersIndex ) const {
std::string BrookStochasticDynamics::_getDerivedParametersString( int derivedParametersIndex ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStochasticDynamics::_getDerivedParametersString";

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

int BrookStochasticDynamics::_initializeStreams( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStochasticDynamics::_initializeStreams";

   BrookOpenMMFloat dangleValue            = (BrookOpenMMFloat) 0.0;

// ---------------------------------------------------------------------------------------

   int sdAtomStreamSize             = getStochasticDynamicsAtomStreamSize();
   int sdAtomStreamWidth            = getStochasticDynamicsAtomStreamWidth();

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

   return DefaultReturnValue;
}

/** 
 * Update sd streams -- called after parameters change
 * 
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookStochasticDynamics::_updateSdStreams( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName       = "BrookStochasticDynamics::_updateSdStreams";

// ---------------------------------------------------------------------------------------

   int sdAtomStreamSize                      = getStochasticDynamicsAtomStreamSize();

   BrookOpenMMFloat* sdpc[2];
   for( int ii = 0; ii < 2; ii++ ){
      sdpc[ii] = new BrookOpenMMFloat[2*sdAtomStreamSize];
      memset( sdpc[ii], 0, 2*sdAtomStreamSize*sizeof( BrookOpenMMFloat ) ); 
   }

   const BrookOpenMMFloat* derivedParameters = getDerivedParameters( );
   int numberOfAtoms                         = getNumberOfAtoms();
   int index                                 = 0;
   for( int ii = 0; ii < numberOfAtoms; ii++, index += 2 ){

      sdpc[0][index]      = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[Yv]) );
      sdpc[0][index+1]    = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[V])  );

      sdpc[1][index]      = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[Yx]) );
      sdpc[1][index+1]    = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[X])  );

   }

   _sdStreams[SDPC1Stream]->loadFromArray( sdpc[0] );
   _sdStreams[SDPC2Stream]->loadFromArray( sdpc[1] );

   for( int ii = 0; ii < 2; ii++ ){
      delete[] sdpc[ii];
   }

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

int BrookStochasticDynamics::_setInverseSqrtMasses( const std::vector<double>& masses ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStochasticDynamics::_setInverseSqrtMasses";

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
 * Setup of StochasticDynamics parameters
 *
 * @param masses                masses
 * @param platform              Brook platform
 *
 * @return nonzero value if error
 *
 * */
    
int BrookStochasticDynamics::setup( const std::vector<double>& masses, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStochasticDynamics::setup";

// ---------------------------------------------------------------------------------------

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

std::string BrookStochasticDynamics::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStochasticDynamics::getContentsString";

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
