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

#include "BrookPlatform.h"
#include "BrookKernelFactory.h"
#include "OpenMMException.h"
#include "kernels.h"
#include "SimTKUtilities/SimTKOpenMMRealType.h"
#include <brook/brook.hpp>
#include <stdlib.h>
#include <sstream>
#include <cctype>
#include <algorithm>

using namespace OpenMM;

/** 
 * Register BrookPlatform
 *
 * @return BrookPlatform instance
 *
 */

BrookPlatform* registerBrookPlatform( void ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::registerBrookPlatform";

// ---------------------------------------------------------------------------------------

   BrookPlatform* platform = new BrookPlatform();
   Platform::registerPlatform(platform);

   return platform;
}

BrookPlatform* staticPlatform = registerBrookPlatform( );

/** 
 * BrookPlatform constructor
 *
 */

BrookPlatform::BrookPlatform( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::BrookPlatform(0)";

// ---------------------------------------------------------------------------------------

   _atomStreamWidth  = DefaultAtomStreamWidth;
   _log              = NULL;

   // get Brook runtime

   char* runtime     = getenv( "brt_runtime" );

   _initializeKernelFactory( );
   _setBrookRuntime( runtime );
}

/** 
 * BrookPlatform constructor
 *
 * @param defaultAtomStreamWidth  stream width
 * @param runtime                 Brook runtime (cal/cpu)
 * @param log                     log file reference
 *
 */

BrookPlatform::BrookPlatform( int atomStreamWidth, const std::string& runtime, FILE* log ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::BrookPlatform(2)";

// ---------------------------------------------------------------------------------------

   _log              = log;
   _atomStreamWidth  = atomStreamWidth;
   _initializeKernelFactory( );
   _setBrookRuntime( runtime );

}

/** 
 * BrookPlatform destructor
 *
 */

BrookPlatform::~BrookPlatform( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::BrookPlatform";

// ---------------------------------------------------------------------------------------

}

/** 
 * Initialize kernel factory
 *
 */

void BrookPlatform::_initializeKernelFactory( void ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::_initializeKernelFactory";

// ---------------------------------------------------------------------------------------

   BrookKernelFactory* factory = new BrookKernelFactory();

   registerKernelFactory( CalcStandardMMForceFieldKernel::Name(), factory);
   // registerKernelFactory( CalcGBSAOBCForceFieldKernel::Name(),    factory);
   //registerKernelFactory( IntegrateVerletStepKernel::Name(),      factory);
   registerKernelFactory( IntegrateLangevinStepKernel::Name(),    factory);
   //registerKernelFactory( IntegrateBrownianStepKernel::Name(),    factory);
   //registerKernelFactory( ApplyAndersenThermostatKernel::Name(),  factory);
   registerKernelFactory( CalcKineticEnergyKernel::Name(),        factory);
   
}

/** 
 * Set & validate runtime
 *
 * @param runtime    Brook runtime (cal/cpu)
 *
 * @throws exception if runtime is invalid
 */

void BrookPlatform::_setBrookRuntime( const std::string& runtime ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookPlatform::_setBrookRuntime";

// ---------------------------------------------------------------------------------------

   // set & validate runtime

   _runtime = runtime;
   std::transform( _runtime.begin(), _runtime.end(), _runtime.begin(), tolower);
   if(  _runtime != "cal" && _runtime != "cpu" ){
      std::stringstream message;
      message << methodName << " Brook runtime=" << _runtime << " not recognized.";
      throw OpenMMException( message.str() );
   }

   if( getLog() ){
      (void) fprintf( getLog(), "%s Brook initializing to runtime=<%s>\n", methodName.c_str(), _runtime.c_str() ); 
      (void) fflush( getLog() );
   }

   brook::initialize( _runtime.c_str(), NULL );

}

/** 
 * Return platform name
 *
 * @return "Brook"
 */
    
std::string BrookPlatform::getName() const {
  return "Brook";
}   

/** 
 * Return platform speed
 *
 * @return speed
 */
    
double BrookPlatform::getSpeed() const {
  return 10.0;
}   

/** 
 * Return true if BrookPlatform supports double precison
 *
 * @return true if BrookPlatform supports double precison
 */

bool BrookPlatform::supportsDoublePrecision( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::supportsDoublePrecision";

// ---------------------------------------------------------------------------------------

    return (sizeof(RealOpenMM) >= sizeof(double));
}

/** 
 * Return Stream factory
 *
 */

const StreamFactory& BrookPlatform::getDefaultStreamFactory( void ) const {
    return _defaultStreamFactory;
}

/** 
 *
 * Static method
 *
 * Return stream size and height given size of array and stream width
 *
 * @param size           size of array
 * @param streamWidth    stream width
 * @param outputHeight   output stream height
 *
 * @return stream size; -1 if streamWidth < 1 || size < 1
 *
 */

int BrookPlatform::getStreamSize( int size, int streamWidth, int* outputHeight ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::getStreamSize";

// ---------------------------------------------------------------------------------------

   if( streamWidth < 1 || size < 1){
      return -1;
   }

   int height = size/streamWidth;
   if( streamWidth*height < size ){
      height++;
   }
   if( outputHeight ){
      *outputHeight = height;
   }
   return height*streamWidth;
}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookPlatform::getLog( void ) const {
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

int BrookPlatform::setLog( FILE* log ){
   _log = log;
   return BrookPlatform::DefaultErrorValue;
}

