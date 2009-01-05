/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "internal/OpenMMContextImpl.h"
#include "System.h"
#include "BrookRemoveCMMotionKernel.h"
#include "BrookStreamInternal.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookRemoveCMMotionKernel constructor
 * 
 * @param name        name of the stream to create
 * @param platform    platform
 * @param openMMBrookInterface      OpenMM-Brook interface
 * @param System                    System reference
 *
 */
  
BrookRemoveCMMotionKernel::BrookRemoveCMMotionKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface, System& system ) :
                              RemoveCMMotionKernel( name, platform ), _openMMBrookInterface( openMMBrookInterface ), _system( system ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookRemoveCMMotionKernel::BrookRemoveCMMotionKernel";

// ---------------------------------------------------------------------------------------

   _frequency                        = 0;
   _log                              = NULL;
   _brookVelocityCenterOfMassRemoval = NULL;

   const BrookPlatform brookPlatform = dynamic_cast<const BrookPlatform&> (platform);
   if( brookPlatform.getLog() != NULL ){
      setLog( brookPlatform.getLog() );
   }   
   _openMMBrookInterface.setNumberOfParticles( system.getNumParticles() );

}

/** 
 * BrookRemoveCMMotionKernel destructor
 * 
 */
  
BrookRemoveCMMotionKernel::~BrookRemoveCMMotionKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookRemoveCMMotionKernel::~BrookRemoveCMMotionKernel";

// ---------------------------------------------------------------------------------------

   delete _brookVelocityCenterOfMassRemoval;
}

/** 
 * Initialize the kernel
 * 
 * @param system   System reference
 * @param force    CMMotionRemover reference
 *
 */

void BrookRemoveCMMotionKernel::initialize( const System& system, const CMMotionRemover& force ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookRemoveCMMotionKernel::initialize";
   static const int PrintOn                 = 0;
   FILE* log                                = getLog();

// ---------------------------------------------------------------------------------------

   _frequency = force.getFrequency();
   std::vector<double> masses;
   masses.resize( system.getNumParticles() );
   for( size_t ii = 0; ii < masses.size(); ii++ ){
      masses[ii] = system.getParticleMass(ii);
   }

   _brookVelocityCenterOfMassRemoval = new BrookVelocityCenterOfMassRemoval();
   _brookVelocityCenterOfMassRemoval->setup( masses, getPlatform() );

   if( PrintOn && log ){
      (void) fprintf( log, "%s\n", methodName.c_str() );
      (void) fflush( log );
   }

   return;

}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookRemoveCMMotionKernel::getLog( void ) const {
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

int BrookRemoveCMMotionKernel::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}


/** 
 * Get COM removal frequency
 * 
 * @return frequency 
 *
 */

int BrookRemoveCMMotionKernel::getFrequency( void ) const {
   return _frequency;
}

/** 
 * Execute kernel
 * 
 * @param context    OpenMMContextImpl reference
 *
 */

void BrookRemoveCMMotionKernel::execute( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookRemoveCMMotionKernel::execute";

// ---------------------------------------------------------------------------------------

   int step      = (int) floor( context.getTime()/context.getIntegrator().getStepSize() );
   int frequency = getFrequency();
   if( frequency <= 0 || (step % frequency) != 0 ){
        return;
   }

   BrookStreamImpl* velocities = _openMMBrookInterface.getParticleVelocities();
   _brookVelocityCenterOfMassRemoval->removeVelocityCenterOfMass( *velocities );

   return;
}
