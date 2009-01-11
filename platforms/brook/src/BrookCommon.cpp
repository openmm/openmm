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

#include <math.h>
#include <stdlib.h>

#include <sstream>

#include "BrookCommon.h"
#include "BrookPlatform.h"
#include "BrookStreamFactory.h"
#include "OpenMMException.h"

using namespace OpenMM;
using namespace std;

// bonded streams
    
const std::string BrookCommon::BondedParticleIndicesStream                        = "BondedParticleIndicesStream";
const std::string BrookCommon::BondedParametersStream                             = "BondedParametersStream";
const std::string BrookCommon::UnrolledForceStream                                = "UnrolledForceStream";
const std::string BrookCommon::BondedChargeStream                                 = "BondedChargeStream";
const std::string BrookCommon::BondedInverseMapStreams                            = "BondedInverseMapStreams";

// non-bonded streams

const std::string BrookCommon::NonBondedExclusionStream                           = "NonBondedExclusionStream";
const std::string BrookCommon::OuterVdwStream                                     = "OuterVdwStream";
const std::string BrookCommon::InnerSigmaStream                                   = "InnerSigmaStream";
const std::string BrookCommon::InnerEpsilonStream                                 = "InnerEpsilonStream";
const std::string BrookCommon::NonBondedChargeStream                              = "NonBondedChargeStream";
const std::string BrookCommon::PartialForceStream                                 = "PartialForceStream";

// OBC Gbsa streams

const std::string BrookCommon::ObcParticleRadiiStream                             = "ObcParticleRadiiStream";
const std::string BrookCommon::ObcScaledParticleRadiiStream                       = "ObcScaledParticleRadiiStream";
const std::string BrookCommon::ObcParticleRadiiWithDielectricOffsetStream         = "ObcParticleRadiiWithDielectricOffsetStream";
const std::string BrookCommon::ObcBornRadiiStream                                 = "ObcBornRadiiStream";
const std::string BrookCommon::ObcBornRadii2Stream                                = "ObcBornRadii2Stream";
const std::string BrookCommon::ObcIntermediateForceStream                         = "ObcIntermediateForceStream";
const std::string BrookCommon::ObcChainStream                                     = "ObcChainStream";

// StochasticDynamics streams

const std::string BrookCommon::SDPC1Stream                                        = "SDPC1Stream";
const std::string BrookCommon::SDPC2Stream                                        = "SDPC2Stream";
const std::string BrookCommon::SD2XStream                                         = "SD2XStream";
const std::string BrookCommon::SD1VStream                                         = "SD1VStream";
const std::string BrookCommon::VPrimeStream                                       = "VPrimeStream";
const std::string BrookCommon::XPrimeStream                                       = "XPrimeStream";
const std::string BrookCommon::InverseMassStream                                  = "InverseMassStream";

// Shake streams

const std::string BrookCommon::ShakeParticleIndicesStream                         = "ShakeParticleIndicesStream";
const std::string BrookCommon::ShakeParticleParameterStream                       = "ShakeParticleParameterStream";
const std::string BrookCommon::ShakeXCons0Stream                                  = "ShakeXCons0Stream";
const std::string BrookCommon::ShakeXCons1Stream                                  = "ShakeXCons1Stream";
const std::string BrookCommon::ShakeXCons2Stream                                  = "ShakeXCons2Stream";
const std::string BrookCommon::ShakeXCons3Stream                                  = "ShakeXCons3Stream";
const std::string BrookCommon::ShakeInverseMapStream                              = "ShakeInverseMapStream";

// Random number streams

const std::string BrookCommon::ShuffleStream                                      = "ShuffleStream";
const std::string BrookCommon::RandomValuesStream                                 = "RandomValuesStream";

// Random number streams

const std::string BrookCommon::BrookVelocityCenterOfMassRemovalWorkStream         = "VelocityCenterOfMassRemovalWorkStream";
const std::string BrookCommon::BrookVelocityCenterOfMassRemovalMassStream         = "VelocityCenterOfMassRemovalMassStream";

/** 
 * Constructor
 * 
 */

BrookCommon::BrookCommon( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::BrookCommon";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _numberOfParticles             = 0;
   _particleSizeModified          = 0;

   _particleStreamWidth           = -1;
   _particleStreamHeight          = -1; 
   _particleStreamSize            = -1;

   _log                           = NULL;
   _isActive                      = 0;

}   
 
/** 
 * Destructor
 * 
 */

BrookCommon::~BrookCommon( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCommon::~BrookCommon";
   //static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

}

/** 
 * Get number of particles
 * 
 * @return  number of particles
 *
 */

int BrookCommon::getNumberOfParticles( void ) const {
   return _numberOfParticles;
}

/** 
 * Get number of particles
 * 
 * @param  numberOfParticles number of particles
 * @return  number of particles
 *
 */

int BrookCommon::setNumberOfParticles( int numberOfParticles ){
   if( numberOfParticles != _numberOfParticles ){
      _particleSizeModified = numberOfParticles;
   }
   _numberOfParticles = numberOfParticles;
   return _numberOfParticles;
}

/** 
 * Get particle stream width
 * 
 * @param platform  platform
 *
 * @return  particle stream width
 *
 */

int BrookCommon::getParticleStreamWidth( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getParticleStreamWidth";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get particle stream width

   if( _particleStreamWidth < 0 ){
      _getParticleStreamDimensions( platform );
   }
   return _particleStreamWidth;
}

/** 
 * Get particle stream width
 * 
 * @return  particle stream width
 *
 */

int BrookCommon::getParticleStreamWidth( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getParticleStreamWidth";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   return _particleStreamWidth;
}

/** 
 * Get particle stream height
 * 
 * @param platform platform
 *
 * @return  particle stream height 
 *
 */

int BrookCommon::getParticleStreamHeight( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getParticleStreamHeight";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get particle stream height

   if( _particleStreamHeight < 0 ){
      _getParticleStreamDimensions( platform );
   }
   return _particleStreamHeight;
}

/** 
 * Get particle stream height
 * 
 * @return  particle stream height 
 *
 */

int BrookCommon::getParticleStreamHeight( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getParticleStreamHeight";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   return _particleStreamHeight;
}

/** 
 * Get particle stream size
 * 
 * @param platform  platform
 *
 * @return  particle stream size
 *
 */

int BrookCommon::getParticleStreamSize( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getParticleStreamSize";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get particle stream size

   if( _particleStreamSize < 0 ){
      _getParticleStreamDimensions( platform );
   }
   return _particleStreamSize;
}

/** 
 * Get particle stream size
 * 
 * @return  particle stream size
 *
 */

int BrookCommon::getParticleStreamSize( void ) const {
   return _particleStreamSize;
}

/** 
 * Get particle stream dimensions
 * 
 * @param platform                  platform
 *
 */

void BrookCommon::_getParticleStreamDimensions( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::_getParticleStreamDimensions";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get particle stream size

   const BrookPlatform brookPlatform            = dynamic_cast<const BrookPlatform&> (platform);
   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory() );
   _particleStreamWidth                         = brookStreamFactory.getDefaultParticleStreamWidth();
   _particleStreamSize                          = brookPlatform.getStreamSize( getNumberOfParticles(), _particleStreamWidth, NULL );
   _particleStreamHeight                        = (int) ( ((float) _particleStreamSize)/( (float) _particleStreamWidth) + 0.001);

   return;
}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookCommon::getLog( void ) const {
   return _log;
}

/** 
 * Set log file reference
 * 
 * @param  log file reference
 *
 * @return  DefaultReturnValue
 *
 */

int BrookCommon::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Get flag signalling whether active
 * 
 * @return  flag signalling whether active
 *
 */

int BrookCommon::isActive( void ) const {
   return _isActive;
}

/** 
 * Set flag signalling whether active
 * 
 * @param  flag signalling whether active
 *
 * @return  DefaultReturnValue
 *
 */

int BrookCommon::setIsActive( int isActive ){
   _isActive= isActive;
   return BrookCommon::DefaultReturnValue;
}

/* 
 * Get contents of object
 *
 * @param tab         tab
 * @param description description
 * @param value       value
 *
 * @return string containing contents
 *
 * */

std::string BrookCommon::_getLine( const std::string& tab, const std::string& description, const std::string& value ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCommon::_getLine";

   static const unsigned int MAX_LINE_CHARS = 256;
   char line[MAX_LINE_CHARS];

// ---------------------------------------------------------------------------------------

   std::stringstream message;
   memset( line, ' ', MAX_LINE_CHARS ); 
#ifdef WIN32
   (void) sprintf_s( line, MAX_LINE_CHARS, "%s %-40s %s", tab.c_str(), description.c_str(), value.c_str() );
#else
   (void) sprintf( line, "%s %-40s %s", tab.c_str(), description.c_str(), value.c_str() );
#endif
   message << std::string( line ) << std::endl;

   return message.str();

}

/* 
 * Given number of stream elements and width, returns the appropriate
 * height of the stream
 *
 * @param streamSize   stream size 
 * @param width        stream width
 *
 * @return stream height
 *
 */

int BrookCommon::getStreamHeight( int streamSize, int streamWidth ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCommon::getStreamHeight";

// ---------------------------------------------------------------------------------------

	int streamHeight = streamSize/streamWidth;

	if( streamSize % streamWidth ){
      streamHeight++;
   }

	return streamHeight;
}

/* 
 * Given number of stream elements, get stream width & height
 *
 * @param streamSize    stream size 
 * @param streamWidth   output stream width
 * @param streamHeight  output stream height
 *
 * @return stream height
 *
 */

void BrookCommon::getStreamDimensions( int streamSize, int *streamWidth, int *streamHeight ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCommon::getStreamDimensions";

// ---------------------------------------------------------------------------------------

	// There are two conditions - stream should be as square 
	// as possible, but should also be multiple of 16 along
	// one dimension
	
	float s = sqrtf( (float) streamSize );

	// find nearest multiple of 16 to the perfect square size

	int low  = ( (int) floor( s/16.0f ) ) * 16;

	if ( !low ) {
		*streamWidth  = 16;
		*streamHeight = getStreamHeight( streamSize, *streamWidth );
	} else { 

	   int high = low + 1;

		// I'm not sure 48 is such a good stream width. Things seem
		// to be faster with 32 or 64. Can make this a special
		// case later.
		
		// Choose low or high depending on which one is 
		// more square

		int htlow  = getStreamHeight( streamSize, low );
		int hthigh = getStreamHeight( streamSize, high );
		
		if ( abs( htlow - low ) < abs( hthigh - high ) ) {
			*streamWidth   = low;
			*streamHeight  = htlow;
		} else {
			*streamWidth   = high;
			*streamHeight  = hthigh;
		}
	}
	return;
}

/* 
 * Allocate array
 *
 * @param length        length of array
 * @param width         width  of array
 *
 * @return ptr to array
 *
 */

RealOpenMM** BrookCommon::allocateRealArray( int length, int width ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCommon::allocateRealArray";

// ---------------------------------------------------------------------------------------

   RealOpenMM** array  = new RealOpenMM*[length];
   RealOpenMM*  buffer = new RealOpenMM[length*width];
   for( int ii = 0; ii < length; ii++ ){
      array[ii] = buffer;
      buffer   += width;
   }
   return array;
}

/* 
 * Free array
 *
 * @param array         array to be freed (assumed allocated using BrookCommon::allocateRealArray
 *
 * @return DefaultReturnValue
 *
 */

int BrookCommon::disposeRealArray( RealOpenMM** array ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCommon::disposeRealArray";

// ---------------------------------------------------------------------------------------

   delete[] array[0];
   delete[] array;

   return DefaultReturnValue;
}

/* 
 * Copy 1D BrookOpenMMFloat* array to 2D array of RealOpenMM
 *
 * @param length        length of array
 * @param width         width  of array
 * @param array1D       array to copy
 *
 * @return ptr to array
 *
 */

RealOpenMM** BrookCommon::copy1DArrayTo2DArray( int length, int width, BrookOpenMMFloat* array1D ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCommon::copy1DArrayTo2DArray";

// ---------------------------------------------------------------------------------------

   RealOpenMM** array  = allocateRealArray( length, width );
   int index           = 0;
   for( int ii = 0; ii < length; ii++ ){
     for( int jj = 0; jj < width; jj++ ){
        array[ii][jj] = static_cast<RealOpenMM> (array1D[index++]);
     }
   }

   return array;
}

