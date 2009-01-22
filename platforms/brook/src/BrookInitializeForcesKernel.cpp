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

#include <cmath>
#include <limits>
#include "OpenMMException.h"
#include <sstream>

#include "BrookStreamImpl.h"
#include "BrookInitializeForcesKernel.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookInitializeForcesKernel constructor
 * 
 * @param name                      kernel name
 * @param platform                  platform
 * @param openMMBrookInterface      OpenMMBrookInterface reference
 * @param system                    System reference
 *
 */

BrookInitializeForcesKernel::BrookInitializeForcesKernel( std::string name, const Platform& platform,
                                                          OpenMMBrookInterface& openMMBrookInterface, System& system ) :
                     InitializeForcesKernel( name, platform ), _openMMBrookInterface( openMMBrookInterface ), _system( system ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookInitializeForcesKernel::BrookInitializeForcesKernel";

// ---------------------------------------------------------------------------------------

   _numberOfParticles                       = 0;
   _log                                     = NULL;

   const BrookPlatform& brookPlatform       = dynamic_cast<const BrookPlatform&> (platform);
   if( brookPlatform.getLog() != NULL ){
      setLog( brookPlatform.getLog() );
   }
      
}   

/** 
 * BrookInitializeForcesKernel destructor
 * 
 */

BrookInitializeForcesKernel::~BrookInitializeForcesKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookInitializeForcesKernel::BrookInitializeForcesKernel";

// ---------------------------------------------------------------------------------------
}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookInitializeForcesKernel::getLog( void ) const {
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

int BrookInitializeForcesKernel::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Initialize 
 * 
 * @param  system System reference
 *
 * @return  DefaultReturnValue
 *
 */

void BrookInitializeForcesKernel::initialize( const System& system ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookInitializeForcesKernel::initialize";

// ---------------------------------------------------------------------------------------

    //FILE* log                 = getLog();

}

/** 
 * Zero forces
 * 
 * @param context OpenMMContextImpl context
 *
 */

void BrookInitializeForcesKernel::execute( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookInitializeForcesKernel::execute";

// ---------------------------------------------------------------------------------------

   _openMMBrookInterface.zeroForces( context );

   // ---------------------------------------------------------------------------------------
}

