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

#include <cmath>
#include <limits>
#include "openmm/OpenMMException.h"
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

