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
#include "OpenMMBrookInterface.h"
#include "internal/OpenMMContextImpl.h"
#include "OpenMMException.h"
#include "kernels.h"
#include "../../reference/src/SimTKUtilities/SimTKOpenMMRealType.h"
#include <brook/brook.hpp>
#include <stdlib.h>
#include <sstream>
#include <cctype>
#include <algorithm>
//#include "/cygdrive/c/Program\ Files/Amd/AMD\ CAL\ 1.2.1_beta/include/cal.h"

using namespace OpenMM;

extern "C" void initOpenMMPlugin() {
/*
   CALuint numberOfDevices = 0;
   // validateCALRuntime()
   if( calDeviceGetCount( &numberOfDevices ) == CAL_RESULT_OK && numberOfDevices > 0 ){
      Platform::registerPlatform(new BrookPlatform());
   }
*/
}

/** 
 * BrookPlatformData constructor
 *
 */

BrookPlatformData::BrookPlatformData( void ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatformData::BrookPlatformData";

// ---------------------------------------------------------------------------------------

      //_forceKernel                 = NULL;

      _removeCOM                   = 0;
      _useOBC                      = 0;
      _hasBonds                    = 0;
      _hasAngles                   = 0;
      _hasPeriodicTorsions         = 0;
      _hasRB                       = 0;
      _hasNonbonded                = 0;
      _cmMotionFrequency           = 0;

}

/** 
 * Get _removeCOM flag
 * 
 * @return  _removeCOM
 *
 */

int BrookPlatformData::removeCOM( void ) const {
   return _removeCOM;
}

/** 
 * Get _useOBC flag
 * 
 * @return  _useOBC
 *
 */

int BrookPlatformData::useOBC( void ) const {
   return _useOBC;
}

/** 
 * Get _hasBonds flag
 * 
 * @return  _hasBonds 
 *
 */

int BrookPlatformData::hasBonds( void ) const {
   return _hasBonds;
}

/** 
 * Get _hasAngles
 * 
 * @return  _hasAngles
 *
 */

int BrookPlatformData::hasAngles( void ) const {
   return _hasAngles;
}

/** 
 * Get _hasPeriodicTorsions
 * 
 * @return  _hasPeriodicTorsions
 *
 */

int BrookPlatformData::hasPeriodicTorsions( void ) const {
   return _hasPeriodicTorsions;
}

/** 
 * Get _hasRB
 * 
 * @return  _hasRB
 *
 */

int BrookPlatformData::hasRB( void ) const {
   return _hasRB;
}

/** 
 * Get _hasNonbonded
 * 
 * @return  _hasNonbonded
 *
 */

int BrookPlatformData::hasNonbonded( void ) const {
   return _hasNonbonded;
}

/** 
 * Get _cmMotionFrequency
 * 
 * @return  _cmMotionFrequency
 *
 */

int BrookPlatformData::cmMotionFrequency( void ) const {
   return _cmMotionFrequency;
}

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

   _particleStreamWidth  = DefaultParticleStreamWidth;
   _minSuggestedThreads  = -1;
   _log                  = NULL;
//_log                  = stderr;

   // get Brook runtime

#ifdef WIN32
   char* runtime;
   size_t numberOfEnv;
    _dupenv_s( &runtime, &numberOfEnv, "brt_runtime" );
#else
   char* runtime     = getenv( "brt_runtime" );
#endif

   // if environment variable 'brt_runtime' not set, default to 'cal' settinh

   if( runtime == NULL ){
      runtime = _strdup( "cal" );
   }

   _initializeKernelFactory( );
   _setBrookRuntime( runtime );

#ifdef WIN32
   free( runtime );
#endif
}

/** 
 * BrookPlatform constructor
 *
 * @param defaultParticleStreamWidth  stream width
 * @param runtime                     Brook runtime (cal/cpu)
 * @param log                         log file reference
 *
 */

BrookPlatform::BrookPlatform( int particleStreamWidth, const std::string& runtime, FILE* log ) : _particleStreamWidth( particleStreamWidth ), _log( log ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::BrookPlatform(2)";

// ---------------------------------------------------------------------------------------

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

   registerKernelFactory( InitializeForcesKernel::Name(),         factory );
   registerKernelFactory( CalcHarmonicBondForceKernel::Name(),    factory );
   registerKernelFactory( CalcHarmonicAngleForceKernel::Name(),   factory );
   registerKernelFactory( CalcPeriodicTorsionForceKernel::Name(), factory );
   registerKernelFactory( CalcRBTorsionForceKernel::Name(),       factory );
   registerKernelFactory( CalcNonbondedForceKernel::Name(),       factory );
   registerKernelFactory( CalcGBSAOBCForceKernel::Name(),         factory );
   registerKernelFactory( IntegrateVerletStepKernel::Name(),      factory );
   registerKernelFactory( IntegrateLangevinStepKernel::Name(),    factory );
   // registerKernelFactory( IntegrateBrownianStepKernel::Name(),    factory );
   //registerKernelFactory( ApplyAndersenThermostatKernel::Name(),  factory );
   registerKernelFactory( CalcKineticEnergyKernel::Name(),        factory );
   registerKernelFactory( RemoveCMMotionKernel::Name(),           factory );
   
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

   // lower case

   std::transform( _runtime.begin(), _runtime.end(), _runtime.begin(), tolower );
   if(  _runtime != "cal" && _runtime != "cpu" ){
      std::stringstream message;
      message << methodName << " Brook runtime=" << _runtime << " not recognized.";
      throw OpenMMException( message.str() );
   }

   if( 1 ){

      //When compiling with cygwin/cl combo, doesn't
      //always work from the environment, so I'm 
      //hardcoding it here. An alternative might be to getenv() in
      //the gromacs code and pass it here. The cygwin getenv() hopefully
      //will work more deterministically.

      char* info_string = NULL;
      int minSuggestedThreads;

      brook::initialize(  _runtime.c_str(), NULL, &info_string, &minSuggestedThreads );

      FILE* log = getLog() ? getLog() : stderr;
      (void) fprintf( log, "Using runtime %s; initializing Brook\n", _runtime.c_str() );
      fprintf( log, "############\n\nBrook info_string:\n%s\n############\n", info_string );
      (void) fflush( log );

      if( minSuggestedThreads > 0 ){
         _minSuggestedThreads = minSuggestedThreads;
      }
   
   } else {

      FILE* log = getLog() ? getLog() : stderr;
      (void) fprintf( log, "%s Brook initializing to runtime=<%s>\n", methodName.c_str(), _runtime.c_str() ); 
      (void) fflush( log );
      brook::initialize( _runtime.c_str(), NULL );

   }

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
 * Get DuplicationFactor 
 *
 * @param numberOfParticles number of particles
 *
 * @return DuplicationFactor
 */
    
int BrookPlatform::getDuplicationFactor( int numberOfParticles ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::getDuplicationFactor";

// ---------------------------------------------------------------------------------------

   // default value

   int duplicationFactor = 4;

   // set only if _minSuggestedThreads is available from board

   if( _minSuggestedThreads > 0 ){
      float threads   = static_cast<float>( _minSuggestedThreads );
      float numP      = static_cast<float>( numberOfParticles );
      float iUnroll   = 4.0f;
      float factor    = (threads*iUnroll)/numP;
      if( (factor*numP) < (threads*iUnroll) ){
         factor += 1.0f;
      }
      if( factor <= 1.0f ){
         duplicationFactor   = 1;
      } else {
         duplicationFactor   = static_cast<int>( ceil( factor*0.25f ) ); 
         duplicationFactor  *= 4; 
      }
   }
   return duplicationFactor;
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
 * Return particle stream width
 *
 * @return particle stream width
 */
    
int BrookPlatform::getParticleStreamWidth() const {
  return _particleStreamWidth;
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

   if( streamWidth < 1 || size < 1 ){
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

/** 
 *
 * This is called whenever a new OpenMMContext is created.  It gives the Platform a chance to initialize
 * the context and store platform-specific data in it.
 *
 */
void BrookPlatform::contextCreated( OpenMMContextImpl& context ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::contextCreated";

// ---------------------------------------------------------------------------------------

   int particles = context.getSystem().getNumParticles();
   OpenMMBrookInterface* openMMBrookInterface = new OpenMMBrookInterface( getParticleStreamWidth(), getDuplicationFactor( particles ) );
//   openMMBrookInterface->setLog( stderr );
 
   context.setPlatformData( openMMBrookInterface );
}

/** 
 *
 * This is called whenever an OpenMMContext is deleted.  It gives the Platform a chance to clean up
 * any platform-specific data that was stored in it.
 *
 */

void BrookPlatform::contextDestroyed( OpenMMContextImpl& context ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookPlatform::contextDestroyed";

// ---------------------------------------------------------------------------------------

}
