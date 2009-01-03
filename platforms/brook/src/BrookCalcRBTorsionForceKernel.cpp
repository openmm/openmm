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
#include "BrookCalcRBTorsionForceKernel.h"

using namespace OpenMM;
using namespace std;

const std::string BrookCalcRBTorsionForceKernel::BondName = "RbTorsion";

/** 
 * BrookCalcRBTorsionForceKernel constructor
 * 
 * @param name                      kernel name
 * @param platform                  platform
 * @param OpenMMBrookInterface      OpenMM-Brook interface
 * @param System                    System reference
 *
 */

BrookCalcRBTorsionForceKernel::BrookCalcRBTorsionForceKernel( std::string name, const Platform& platform,
                                                              OpenMMBrookInterface& openMMBrookInterface, System& system ) :
                     CalcRBTorsionForceKernel( name, platform ), _openMMBrookInterface( openMMBrookInterface ), _system( system ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcRBTorsionForceKernel::BrookCalcRBTorsionForceKernel";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _brookBondParameters              = NULL;
   _log                              = NULL;
   _openMMBrookInterface.setNumberOfParticles( system.getNumParticles() );

   const BrookPlatform brookPlatform = dynamic_cast<const BrookPlatform&> (platform);
   if( brookPlatform.getLog() != NULL ){
      setLog( brookPlatform.getLog() );
   }
      
}   

/** 
 * BrookCalcRBTorsionForceKernel destructor
 * 
 */

BrookCalcRBTorsionForceKernel::~BrookCalcRBTorsionForceKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcRBTorsionForceKernel::BrookCalcRBTorsionForceKernel";
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

FILE* BrookCalcRBTorsionForceKernel::getLog( void ) const {
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

int BrookCalcRBTorsionForceKernel::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Initialize the kernel, setting up the values of all the force field parameters.
 * 
 * @param system                    System reference
 * @param force                     RbTorsionForce reference
 *
 */

void BrookCalcRBTorsionForceKernel::initialize( const System& system, const RBTorsionForce& force ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcRBTorsionForceKernel::initialize";

// ---------------------------------------------------------------------------------------

   FILE* log                                = getLog();

   // ---------------------------------------------------------------------------------------

   // create _brookBondParameters object containing particle indices/parameters

   int numberOfBonds         = force.getNumTorsions();

   if( _brookBondParameters ){
      delete _brookBondParameters;
   }
   _brookBondParameters = new BrookBondParameters( BondName, NumberOfParticlesInBond, NumberOfParametersInBond, numberOfBonds, getLog() );

   for( int ii = 0; ii < numberOfBonds; ii++ ){

      int particle1, particle2, particle3, particle4;
      double c0, c1, c2, c3, c4, c5;

      int particles[NumberOfParticlesInBond];
      double parameters[NumberOfParametersInBond];

      force.getTorsionParameters( ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 ); 

      particles[0]    = particle1;
      particles[1]    = particle2;
      particles[2]    = particle3;
      particles[3]    = particle4;
 
      parameters[0]   = c0;
      parameters[1]   = c1;
      parameters[2]   = c2;
      parameters[3]   = c3;
      parameters[4]   = c4;
      parameters[5]   = c5;

      _brookBondParameters->setBond( ii, particles, parameters );
   }   
   _openMMBrookInterface.setRBTorsionForceParameters( _brookBondParameters );
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

void BrookCalcRBTorsionForceKernel::executeForces( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName   = "BrookCalcRBTorsionForceKernel::executeForces";

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

double BrookCalcRBTorsionForceKernel::executeEnergy( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcRBTorsionForceKernel::executeEnergy";

// ---------------------------------------------------------------------------------------

   if( _openMMBrookInterface.getTriggerEnergyKernel() == this ){
      return (double) _openMMBrookInterface.computeEnergy( context, _system );
   } else {
      return 0.0;
   }

}
