/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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
#include "OpenMMException.h"
#include "BrookStreamFactory.h"
#include "BrookFloatStreamImpl.h"
#include "BrookIntStreamImpl.h"

using namespace OpenMM;

const std::string BrookStreamFactory::AtomPositions              = "atomPositions";
const std::string BrookStreamFactory::AtomVelocities             = "atomVelocities";
const std::string BrookStreamFactory::AtomForces                 = "atomForces";

// bonded streams
                                    
const std::string BrookStreamFactory::BondedAtomIndicesStream    = "BondedAtomIndicesStream";
const std::string BrookStreamFactory::BondedParametersStream     = "BondedParametersStream";
const std::string BrookStreamFactory::UnrolledForceStream        = "UnrolledForceStream";
const std::string BrookStreamFactory::BondedChargeStream         = "BondedChargeStream";
const std::string BrookStreamFactory::BondedInverseMapStreams    = "BondedInverseMapStreams";

// non-bonded streams

const std::string BrookStreamFactory::NonBondedExclusionStream   = "NonBondedExclusionStream";
const std::string BrookStreamFactory::NonBondedVdwStream         = "NonBondedVdwStream";

/** 
 * BrookStreamFactory constructor
 * 
 * @return BrookStreamFactory
 */

BrookStreamFactory::BrookStreamFactory( void ){

	double defaultDangleValue                       = 1.0e+38;
	int    defaultStreamWidth                       = 32;

   _streamInfoMap[AtomPositions]                   = new BrookStreamInfo( AtomPositions,             defaultStreamWidth, defaultDangleValue );
   _streamInfoMap[AtomVelocities]                  = new BrookStreamInfo( AtomVelocities,            defaultStreamWidth, defaultDangleValue );
   _streamInfoMap[AtomForces]                      = new BrookStreamInfo( AtomForces,                defaultStreamWidth, defaultDangleValue );

   // bonded streams

   _streamInfoMap[BondedAtomIndicesStream]         = new BrookStreamInfo( BondedAtomIndicesStream,   defaultStreamWidth, defaultDangleValue );
   _streamInfoMap[BondedParametersStream]          = new BrookStreamInfo( BondedParametersStream,    defaultStreamWidth, defaultDangleValue );
   _streamInfoMap[UnrolledForceStream]             = new BrookStreamInfo( UnrolledForceStream,       defaultStreamWidth, defaultDangleValue );
   _streamInfoMap[BondedChargeStream]              = new BrookStreamInfo( BondedChargeStream,        defaultStreamWidth, defaultDangleValue );
   _streamInfoMap[BondedInverseMapStreams]         = new BrookStreamInfo( BondedInverseMapStreams,   defaultStreamWidth, defaultDangleValue );

   _streamInfoMap[NonBondedExclusionStream]        = new BrookStreamInfo( NonBondedExclusionStream,  defaultStreamWidth, defaultDangleValue );
   _streamInfoMap[NonBondedVdwStream]              = new BrookStreamInfo( NonBondedVdwStream,        defaultStreamWidth, defaultDangleValue );
}

/** 
 * BrookStreamFactory destructor
 * 
 */

BrookStreamFactory::~BrookStreamFactory( void ){
   //_streamInfoMap[UnrolledForceStream]      = new BrookStreamInfo( UnrolledForceStream, 32, defaultDangleValue );
}

/** 
 * Get BrookStreamInfo reference given stream name
 *
 * @param name stream name
 * 
 * @return BrookStreamInfo  -- look up streamInfo object given name; return NULL if name not recognized
 */

BrookStreamInfo* BrookStreamFactory::getBrookStreamInfo( std::string name ) const {

   if( _streamInfoMap.find( name ) == _streamInfoMap.end() ){
      return NULL;
   }
   return _streamInfoMap.find( name )->second;
}

/** 
 * Create StreamImpl
 *
 * @param name     stream name
 * @param size     stream size
 * @param type     data type (float, float2, ...)
 * @param platform platform reference
 * @param context  context (currently ignored)
 * 
 * @return StreamImpl
 */

StreamImpl* BrookStreamFactory::createStreamImpl( std::string name, int size, Stream::DataType type,
                                                  const Platform& platform, OpenMMContextImpl& context ) const {

   return BrookStreamFactory::createStreamImplCommon( name, size, type, platform );

}

/** 
 * Create StreamImpl
 *
 * @param name     stream name
 * @param size     stream size
 * @param type     data type (float, float2, ...)
 * @param platform platform reference
 * 
 * @return StreamImpl
 */

StreamImpl* BrookStreamFactory::createStreamImpl( std::string name, int size, Stream::DataType type,
                                                  const Platform& platform ) const {

   return BrookStreamFactory::createStreamImplCommon( name, size, type, platform );
}

/** 
 * Create StreamImpl
 *
 * @param name     stream name
 * @param size     stream size
 * @param type     data type (float, float2, ...)
 * @param platform platform reference
 * 
 * @return StreamImpl
 */

StreamImpl* BrookStreamFactory::createStreamImplCommon( std::string name, int size, Stream::DataType type,
                                                        const Platform& platform ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamFactory::createStreamImplCommon";
   //static const int debug                   = 0;

// ---------------------------------------------------------------------------------------

   // get stream width & dangle value

   BrookStreamInfo* streamInfo = getBrookStreamInfo( name );
   if( streamInfo == NULL ){
      std::stringstream message;
      message << methodName << " stream=" << name << " not registered.";
      throw OpenMMException( message.str() );
   }

   int streamWidth    = streamInfo->getStreamWidth();
   double dangleValue = streamInfo->getDangleValue();
      
   switch ( type ){

      case Stream::Float:
      case Stream::Float2:
      case Stream::Float3:
      case Stream::Float4:
      case Stream::Double:
      case Stream::Double2:
      case Stream::Double3:
      case Stream::Double4:
          return new BrookFloatStreamImpl( name, size, type, platform, streamWidth, dangleValue );
          break;

      case Stream::Integer:
      case Stream::Integer2:
      case Stream::Integer3:
      case Stream::Integer4:
          return new BrookIntStreamImpl( name, size, type, platform );
          break;
   }

   std::stringstream message;
   message << methodName << " type=" << type << " for stream=" << name << " is invalid.";
   throw OpenMMException( message.str() );

}
