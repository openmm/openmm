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
#include "BrookCalcGBSAOBCForceKernel.h"
#include "kernels/kgbsa.h"
#include "kernels/kforce.h"
#include "math.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookCalcGBSAOBCForceKernel constructor
 * 
 * @param name                      kernel name
 * @param platform                  platform
 * @param OpenMMBrookInterface      OpenMM-Brook interface
 * @param System                    System reference
 *
 */

BrookCalcGBSAOBCForceKernel::BrookCalcGBSAOBCForceKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface, System& system ) :
                             CalcGBSAOBCForceKernel( name, platform ), _openMMBrookInterface( openMMBrookInterface ), _system( system ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcGBSAOBCForceKernel::BrookCalcGBSAOBCForceKernel";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _numberOfParticles                = 0;
   _log                              = NULL;

   const BrookPlatform brookPlatform = dynamic_cast<const BrookPlatform&> (platform);
   if( brookPlatform.getLog() != NULL ){
      setLog( brookPlatform.getLog() );
   }
      
   _openMMBrookInterface.setNumberOfParticles( system.getNumParticles() );
}   

/** 
 * BrookCalcGBSAOBCForceKernel destructor
 * 
 */

BrookCalcGBSAOBCForceKernel::~BrookCalcGBSAOBCForceKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcGBSAOBCForceKernel::BrookCalcGBSAOBCForceKernel";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookCalcGBSAOBCForceKernel::getLog( void ) const {
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

int BrookCalcGBSAOBCForceKernel::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Initialize the kernel, setting up the values of all the force field parameters.
 * 
 * @param system     system this kernel will be applied to
 * @param force      GBSAOBCForce this kernel will be used for
 *
 */

void BrookCalcGBSAOBCForceKernel::initialize( const System& system, const GBSAOBCForce& force ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcGBSAOBCForceKernel::initialize";
   int printOn                              = 0;
   FILE* log;

// ---------------------------------------------------------------------------------------

   if( printOn && getLog() ){
       log = getLog();
   } else {
      printOn = 0;
   } 

   // ---------------------------------------------------------------------------------------

   BrookGbsa& brookGbsa      = _openMMBrookInterface.getBrookGbsa();
    
   // get parameters from force object
   // and initialize brookGbsa

   _numberOfParticles = system.getNumParticles();
   std::vector<std::vector<double> > particleParameters( _numberOfParticles );

   for( int ii = 0; ii < _numberOfParticles; ii++ ){

      double charge, radius, scalingFactor;
      force.getParticleParameters( ii, charge, radius, scalingFactor );

      particleParameters[ii].push_back( charge );
      particleParameters[ii].push_back( radius );
      particleParameters[ii].push_back( scalingFactor );
   }   
   brookGbsa.setup( particleParameters, force.getSolventDielectric(), force.getSoluteDielectric(), getPlatform() );
   brookGbsa.setIsActive( 1 );

   _openMMBrookInterface.setTriggerForceKernel( this );
   _openMMBrookInterface.setTriggerEnergyKernel( this );

   if( printOn ){
      std::string contents = brookGbsa.getContentsString( ); 
      (void) fprintf( log, "%s brookGbsa::contents\n%s", methodName.c_str(), contents.c_str() );
      (void) fflush( log );
   }

   // ---------------------------------------------------------------------------------------
    
}

/** 
 * Compute forces given particle coordinates
 * 
 * @param context OpenMMContextImpl context
 *
 */

void BrookCalcGBSAOBCForceKernel::executeForces( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName   = "BrookCalcGBSAOBCForceKernel::executeForces";

// ---------------------------------------------------------------------------------------

   if( _openMMBrookInterface.getTriggerForceKernel() == this ){
      _openMMBrookInterface.computeForces( context );
   }   

}

/**
 * Execute the kernel to calculate the OBC energy
 * 
 * @param context OpenMMContextImpl context
 *
 * @return energy
 *
 */

double BrookCalcGBSAOBCForceKernel::executeEnergy( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcGBSAOBCForceKernel::executeEnergy";

// ---------------------------------------------------------------------------------------

   if( _openMMBrookInterface.getTriggerEnergyKernel() == this ){
      return (double) _openMMBrookInterface.computeEnergy( context, _system );
   } else {
      return 0.0;
   }

}
