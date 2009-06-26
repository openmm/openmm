/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs, Mike Houston, Peter Eastman                      *
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
#include "BrookUpdateTimeKernel.h"

#include <sstream>

using namespace OpenMM;
using namespace std;

/**
 * BrookUpdateTimeKernel constructor
 *
 * @param name                      kernel name
 * @param platform                  platform
 * @param openMMBrookInterface      OpenMMBrookInterface reference
 *
 */

BrookUpdateTimeKernel::BrookUpdateTimeKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface ) :
                     UpdateTimeKernel( name, platform ), _openMMBrookInterface(openMMBrookInterface){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookUpdateTimeKernel::BrookUpdateTimeKernel";

// ---------------------------------------------------------------------------------------

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

BrookUpdateTimeKernel::~BrookUpdateTimeKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookUpdateTimeKernel::~BrookUpdateTimeKernel";

// ---------------------------------------------------------------------------------------
}

/**
 * Get log file reference
 *
 * @return  log file reference
 *
 */

FILE* BrookUpdateTimeKernel::getLog( void ) const {
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

int BrookUpdateTimeKernel::setLog( FILE* log ){
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

void BrookUpdateTimeKernel::initialize( const System& system ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookUpdateTimeKernel::initialize";

// ---------------------------------------------------------------------------------------

    //FILE* log                 = getLog();

}

/**
* Get the current time (in picoseconds).
*
* @param context    the context in which to execute this kernel
*/

double BrookUpdateTimeKernel::getTime(const OpenMMContextImpl& context) const {
    return _openMMBrookInterface.getTime();
}

/**
* Set the current time (in picoseconds).
*
* @param context    the context in which to execute this kernel
* @param time       the time
*/

void BrookUpdateTimeKernel::setTime(OpenMMContextImpl& context, double time) {
    _openMMBrookInterface.setTime(time);
}
