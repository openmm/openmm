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

#include "OpenMMException.h"
#include <sstream>
#include "BrookCalcHarmonicBondForceKernel.h"

using namespace OpenMM;
using namespace std;

const std::string BrookCalcHarmonicBondForceKernel::BondName = "HarmonicBond";

/** 
 * BrookCalcHarmonicBondForceKernel constructor
 * 
 * @param name                      kernel name
 * @param platform                  platform
 * @param OpenMMBrookInterface      OpenMM-Brook interface
 * @param System                    System reference
 *
 */

BrookCalcHarmonicBondForceKernel::BrookCalcHarmonicBondForceKernel( std::string name, const Platform& platform,
                                                                      OpenMMBrookInterface& openMMBrookInterface, System& system ) :
                     CalcHarmonicBondForceKernel( name, platform ), _openMMBrookInterface( openMMBrookInterface ), _system( system ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcHarmonicBondForceKernel::BrookCalcHarmonicBondForceKernel";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _brookBondParameters              = NULL;
   _log                              = NULL;

   const BrookPlatform brookPlatform = dynamic_cast<const BrookPlatform&> (platform);
   if( brookPlatform.getLog() != NULL ){
      setLog( brookPlatform.getLog() );
   }
      
}   

/** 
 * BrookCalcHarmonicBondForceKernel destructor
 * 
 */

BrookCalcHarmonicBondForceKernel::~BrookCalcHarmonicBondForceKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcHarmonicBondForceKernel::BrookCalcHarmonicBondForceKernel";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   delete _brookBondParameters;
}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookCalcHarmonicBondForceKernel::getLog( void ) const {
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

int BrookCalcHarmonicBondForceKernel::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Initialize the kernel, setting up the values of all the force field parameters.
 * 
 * @param system                    System reference
 * @param force                     HarmonicBondForce reference
 *
 */

void BrookCalcHarmonicBondForceKernel::initialize( const System& system, const HarmonicBondForce& force ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcHarmonicBondForceKernel::initialize";

// ---------------------------------------------------------------------------------------

   FILE* log                 = getLog();

   // ---------------------------------------------------------------------------------------

   // create _brookBondParameters object containing particle indices/parameters

   int numberOfBonds         = force.getNumBonds();

   if( _brookBondParameters ){
      delete _brookBondParameters;
   }
   _brookBondParameters = new BrookBondParameters( BondName, NumberOfParticlesInBond, NumberOfParametersInBond, numberOfBonds, getLog() );

   for( int ii = 0; ii < numberOfBonds; ii++ ){

      int particle1, particle2;
      double length, k;

      int particles[NumberOfParticlesInBond];
      double parameters[NumberOfParametersInBond];

      force.getBondParameters( ii, particle1, particle2, length, k ); 
      particles[0]    = particle1;
      particles[1]    = particle2;
 
      parameters[0]   = length;
      parameters[1]   = k;

      _brookBondParameters->setBond( ii, particles, parameters );
   }   
   _openMMBrookInterface.setHarmonicBondForceParameters( _brookBondParameters );
   _openMMBrookInterface.setTriggerForceKernel( this );
   _openMMBrookInterface.setTriggerEnergyKernel( this );

   if( log ){
      std::string contents = _brookBondParameters->getContentsString( ); 
      (void) fprintf( log, "%s contents\n%s", methodName.c_str(), contents.c_str() );
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

void BrookCalcHarmonicBondForceKernel::executeForces( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName   = "BrookCalcHarmonicBondForceKernel::executeForces";

// ---------------------------------------------------------------------------------------

   if( _openMMBrookInterface.getTriggerForceKernel() == this ){
      _openMMBrookInterface.computeForces( context );
   }

   return;

   // ---------------------------------------------------------------------------------------
}

/**
 * Execute the kernel to calculate the energy
 * 
 * @param context OpenMMContextImpl context
 *
 * @return  potential energy
 *
 */

double BrookCalcHarmonicBondForceKernel::executeEnergy( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcHarmonicBondForceKernel::executeEnergy";

// ---------------------------------------------------------------------------------------

   if( _openMMBrookInterface.getTriggerEnergyKernel() == this ){
      return (double) _openMMBrookInterface.computeEnergy( context );
   } else {
      return 0.0;
   }

}
