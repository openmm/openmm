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
#include "BrookStreamImpl.h"

using namespace OpenMM;

const std::string BrookStreamFactory::AtomPositions              = "atomPositions";
const std::string BrookStreamFactory::AtomVelocities             = "atomVelocities";
const std::string BrookStreamFactory::AtomForces                 = "atomForces";

const double DefaultDangleValue                                  = 1.0e+38;
/** 
 * BrookStreamFactory constructor
 * 
 * @return BrookStreamFactory
 */

BrookStreamFactory::BrookStreamFactory( void ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::BrookStreamFactory";

// ---------------------------------------------------------------------------------------

	_defaultDangleValue                      = 1.0e+38;
	_defaultAtomStreamWidth                  = 32;

}

/** 
 * BrookStreamFactory destructor
 * 
 */

BrookStreamFactory::~BrookStreamFactory( void ){
}

/** 
 * Get atom stream width
 * 
 * @return atomStreamWidth
 *
 */

int BrookStreamFactory::getDefaultAtomStreamWidth( void ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::getDefaultAtomStreamWidth";

// ---------------------------------------------------------------------------------------

   return _defaultAtomStreamWidth;
}

/** 
 * Set atom stream width
 * 
 * @param atomStreamWidth  atom stream width
 *
 * @return DefaultReturnValue
 *
 * @throw OpenMMException if atomStreamWidth < 1
 *
 */

int BrookStreamFactory::setDefaultAtomStreamWidth( int atomStreamWidth ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamFactory::setDefaultAtomStreamWidth";

// ---------------------------------------------------------------------------------------

   // validate atom stream width

   if( atomStreamWidth < 1 ){
      std::stringstream message;
      message << methodName << " atomStreamWidth=" << atomStreamWidth << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   _defaultAtomStreamWidth = atomStreamWidth;

   return DefaultReturnValue;

}

/** 
 * get default dangle value
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
                                                  const Platform& platform, OpenMMContextImpl& context ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamFactory::createStreamImpl";

// ---------------------------------------------------------------------------------------


   // stream width hould be based on name & value set in platform; for now only atom stream types

   int streamWidth = getDefaultAtomStreamWidth();

   return new BrookStreamImpl( name, size, streamWidth, type, platform );

}
