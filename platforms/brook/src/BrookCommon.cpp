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
    
const std::string BrookCommon::BondedAtomIndicesStream                            = "BondedAtomIndicesStream";
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

const std::string BrookCommon::ObcAtomicRadiiStream                               = "ObcAtomicRadiiStream";
const std::string BrookCommon::ObcScaledAtomicRadiiStream                         = "ObcScaledAtomicRadiiStream";
const std::string BrookCommon::ObcAtomicRadiiWithDielectricOffsetStream           = "ObcAtomicRadiiWithDielectricOffsetStream";
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

const std::string BrookCommon::ShakeAtomIndicesStream                             = "ShakeAtomIndicesStream";
const std::string BrookCommon::ShakeAtomParameterStream                           = "ShakeAtomParameterStream";
const std::string BrookCommon::ShakeXCons0Stream                                  = "ShakeXCons0Stream";
const std::string BrookCommon::ShakeXCons1Stream                                  = "ShakeXCons1Stream";
const std::string BrookCommon::ShakeXCons2Stream                                  = "ShakeXCons2Stream";
const std::string BrookCommon::ShakeXCons3Stream                                  = "ShakeXCons3Stream";
const std::string BrookCommon::ShakeInverseMapStream                              = "ShakeInverseMapStream";

// Random number streams

const std::string BrookCommon::ShuffleStream                                      = "ShuffleStream";

/** 
 * Constructor
 * 
 */

BrookCommon::BrookCommon(  ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::BrookCommon";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _numberOfAtoms             = 0;
   _atomSizeModified          = 0;

   _atomStreamWidth           = -1;
   _atomStreamHeight          = -1; 
   _atomStreamSize            = -1;

   _log                       = NULL;
   _verbosity                 = 0;

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
 * Get number of atoms
 * 
 * @return  number of atoms
 *
 */

int BrookCommon::getNumberOfAtoms( void ) const {
   return _numberOfAtoms;
}

/** 
 * Get number of atoms
 * 
 * @param  numberOfAtoms number of atoms
 * @return  number of atoms
 *
 */

int BrookCommon::setNumberOfAtoms( int numberOfAtoms ){
   if( numberOfAtoms != _numberOfAtoms ){
      _atomSizeModified = numberOfAtoms;
   }
   _numberOfAtoms = numberOfAtoms;
   return _numberOfAtoms;
}

/** 
 * Get atom stream width
 * 
 * @param platform  platform
 *
 * @return  atom stream width
 *
 */

int BrookCommon::getAtomStreamWidth( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getAtomStreamWidth";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get atom stream width

   if( _atomStreamWidth < 0 ){
      _getAtomStreamDimensions( platform );
   }
   return _atomStreamWidth;
}

/** 
 * Get atom stream width
 * 
 * @return  atom stream width
 *
 */

int BrookCommon::getAtomStreamWidth( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getAtomStreamWidth";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   return _atomStreamWidth;
}

/** 
 * Get atom stream height
 * 
 * @param platform platform
 *
 * @return  atom stream height 
 *
 */

int BrookCommon::getAtomStreamHeight( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getAtomStreamHeight";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get atom stream height

   if( _atomStreamHeight < 0 ){
      _getAtomStreamDimensions( platform );
   }
   return _atomStreamHeight;
}

/** 
 * Get atom stream height
 * 
 * @return  atom stream height 
 *
 */

int BrookCommon::getAtomStreamHeight( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getAtomStreamHeight";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   return _atomStreamHeight;
}

/** 
 * Get atom stream size
 * 
 * @param platform  platform
 *
 * @return  atom stream size
 *
 */

int BrookCommon::getAtomStreamSize( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getAtomStreamSize";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get atom stream size

   if( _atomStreamSize < 0 ){
      _getAtomStreamDimensions( platform );
   }
   return _atomStreamSize;
}

/** 
 * Get atom stream size
 * 
 * @return  atom stream size
 *
 */

int BrookCommon::getAtomStreamSize( void ) const {
   return _atomStreamSize;
}

/** 
 * Get atom stream dimensions
 * 
 * @param platform                  platform
 *
 */

void BrookCommon::_getAtomStreamDimensions( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::_getAtomStreamDimensions";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get atom stream size

   const BrookPlatform brookPlatform            = dynamic_cast<const BrookPlatform&> (platform);
   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory() );
   _atomStreamWidth                             = brookStreamFactory.getDefaultAtomStreamWidth();
   _atomStreamSize                              = brookPlatform.getStreamSize( getNumberOfAtoms(), _atomStreamWidth, NULL );
   _atomStreamHeight                            = (int) ( ((float) _atomStreamSize)/( (float) _atomStreamWidth) + 0.001);

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
 * Get verbosity 
 * 
 * @return   verbosity
 *
 */

int BrookCommon::getVerbosity( void ) const {
   return _verbosity;
}

/** 
 * Set verbosity
 * 
 * @param  verbosity
 *
 * @return  DefaultReturnValue
 *
 */

int BrookCommon::setVerbosity( int verbosity ){
   _verbosity = verbosity;
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


