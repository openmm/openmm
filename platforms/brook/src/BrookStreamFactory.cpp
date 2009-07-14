/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs, Mike Houston                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/OpenMMException.h"
#include "BrookStreamFactory.h"
#include "BrookStreamImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "OpenMMBrookInterface.h"

#include <sstream>

using namespace OpenMM;

const std::string BrookStreamFactory::ParticlePositions              = "particlePositions";
const std::string BrookStreamFactory::ParticleVelocities             = "particleVelocities";
const std::string BrookStreamFactory::ParticleForces                 = "particleForces";

const double DefaultDangleValue                                      = 1.0e+08;

/** 
 * BrookStreamFactory constructor
 * 
 * @return BrookStreamFactory
 */

BrookStreamFactory::BrookStreamFactory( void ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::BrookStreamFactory";

// ---------------------------------------------------------------------------------------

	_defaultDangleValue                      = 1.0e+08;
	_defaultParticleStreamWidth              = DefaultStreamParticleWidth;
   _defaultStreamRandomNumberWidth          = DefaultStreamRandomNumberWidth;
   _defaultStreamRandomNumberSize           = DefaultStreamRandomNumberSize;

}

/** 
 * BrookStreamFactory destructor
 * 
 */

BrookStreamFactory::~BrookStreamFactory( ){
}

/** 
 * Get particle stream width
 * 
 * @return particleStreamWidth
 *
 */

int BrookStreamFactory::getDefaultParticleStreamWidth( void ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::getDefaultParticleStreamWidth";

// ---------------------------------------------------------------------------------------

   return _defaultParticleStreamWidth;
}

/** 
 * Set particle stream width
 * 
 * @param particleStreamWidth  particle stream width
 *
 * @return DefaultReturnValue
 *
 * @throw OpenMMException if particleStreamWidth < 1
 *
 */

int BrookStreamFactory::setDefaultParticleStreamWidth( int particleStreamWidth ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamFactory::setDefaultParticleStreamWidth";

// ---------------------------------------------------------------------------------------

   // validate particle stream width

   if( particleStreamWidth < 1 ){
      std::stringstream message;
      message << methodName << " particleStreamWidth=" << particleStreamWidth << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   _defaultParticleStreamWidth = particleStreamWidth;

   return DefaultReturnValue;

}

/** 
 * Get randomNumber stream width
 * 
 * @return randomNumberStreamWidth
 *
 */

int BrookStreamFactory::getDefaultRandomNumberStreamWidth( void ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::getDefaultRandomNumberStreamWidth";

// ---------------------------------------------------------------------------------------

   return _defaultStreamRandomNumberWidth;
}

/** 
 * Set randomNumber stream width
 * 
 * @param randomNumberStreamWidth  randomNumber stream width
 *
 * @return DefaultReturnValue
 *
 * @throw OpenMMException if randomNumberStreamWidth < 1
 *
 */

int BrookStreamFactory::setDefaultRandomNumberStreamWidth( int randomNumberStreamWidth ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamFactory::setDefaultRandomNumberStreamWidth";

// ---------------------------------------------------------------------------------------

   // validate randomNumber stream width

   if( randomNumberStreamWidth < 1 ){
      std::stringstream message;
      message << methodName << " randomNumberStreamWidth=" << randomNumberStreamWidth << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   _defaultStreamRandomNumberWidth = randomNumberStreamWidth;

   return DefaultReturnValue;

}

/*
 * Get randomNumber stream size
 * 
 * @return randomNumberStreamSize
 *
 */

int BrookStreamFactory::getDefaultRandomNumberStreamSize( void ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::getDefaultRandomNumberStreamSize";

// ---------------------------------------------------------------------------------------

   return _defaultStreamRandomNumberSize;
}

/** 
 * Set randomNumber stream size
 * 
 * @param randomNumberStreamSize  randomNumber stream size
 *
 * @return DefaultReturnValue
 *
 * @throw OpenMMException if randomNumberStreamSize < 1
 *
 */

int BrookStreamFactory::setDefaultRandomNumberStreamSize( int randomNumberStreamSize ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamFactory::setDefaultRandomNumberStreamSize";

// ---------------------------------------------------------------------------------------

   // validate randomNumber stream size

   if( randomNumberStreamSize < 1 ){
      std::stringstream message;
      message << methodName << " randomNumberStreamSize=" << randomNumberStreamSize << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   _defaultStreamRandomNumberSize = randomNumberStreamSize;

   return DefaultReturnValue;

}

/** 
 * Get default dangle value
 * 
 * @return default dangle value
 *
 */

double BrookStreamFactory::getDefaultDangleValue( void ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::getDefaultDangleValue";

// ---------------------------------------------------------------------------------------

   return _defaultDangleValue;
}

/** 
 * Set default dangle value
 * 
 * @param DefaultDangleValue default dangle value
 *
 * @return DefaultReturnValue
 *
 */

int BrookStreamFactory::setDefaultDangleValue( double defaultDangleValue ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::setDefaultDangleValue";

// ---------------------------------------------------------------------------------------

   _defaultDangleValue = defaultDangleValue;

   return DefaultReturnValue;

}

/** 
 * Create StreamInternal
 *
 * @param name     stream name
 * @param size     stream size
 * @param type     data type (float, float2, ...)
 * @param platform platform reference
 * @param context  context (currently ignored)
 * 
 * @return StreamInternal
 */

StreamImpl* BrookStreamFactory::createStreamImpl( std::string name, int size, Stream::DataType type,
                                                  const Platform& platform, ContextImpl& context ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::createStreamImpl";

// ---------------------------------------------------------------------------------------

   // stream width hould be based on name & value set in platform; for now only particle stream types

   int streamWidth                             = getDefaultParticleStreamWidth();

   BrookStreamImpl* brookStreamImpl            = new BrookStreamImpl( name, size, streamWidth, type, platform );

   OpenMMBrookInterface& openMMBrookInterface  = *static_cast<OpenMMBrookInterface*>(context.getPlatformData());

   if( name == ParticlePositions ){
      openMMBrookInterface.setParticlePositions( brookStreamImpl );
   } else if( name == ParticleVelocities ){
      openMMBrookInterface.setParticleVelocities( brookStreamImpl );
   } else if( name == ParticleForces ){
      openMMBrookInterface.setParticleForces( brookStreamImpl );
   }

   return brookStreamImpl;

}

/** 
 * Create StreamInternal
 *
 * @param name     stream name
 * @param size     stream size
 * @param type     data type (float, float2, ...)
 * @param platform platform reference
 * 
 * @return StreamInternal
 */

StreamImpl* BrookStreamFactory::createStreamImpl( std::string name, int size, Stream::DataType type,
                                                  const Platform& platform ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::createStreamImpl";

// ---------------------------------------------------------------------------------------

   // stream width hould be based on name & value set in platform; for now only particle stream types

   int streamWidth                             = getDefaultParticleStreamWidth();

   BrookStreamImpl* brookStreamImpl            = new BrookStreamImpl( name, size, streamWidth, type, platform );

   return brookStreamImpl;

}
