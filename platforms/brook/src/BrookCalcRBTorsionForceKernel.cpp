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

#include "openmm/OpenMMException.h"
#include "BrookCalcRBTorsionForceKernel.h"
#include <sstream>

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

   _brookBondParameters               = NULL;
   _log                               = NULL;
   _openMMBrookInterface.setNumberOfParticles( system.getNumParticles() );

   const BrookPlatform& brookPlatform = dynamic_cast<const BrookPlatform&> (platform);
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
 * @param context ContextImpl context
 *
 */

void BrookCalcRBTorsionForceKernel::executeForces( ContextImpl& context ){

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
 * @param context ContextImpl context
 *
 * @return  potential energy
 *
 */

double BrookCalcRBTorsionForceKernel::executeEnergy( ContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcRBTorsionForceKernel::executeEnergy";

// ---------------------------------------------------------------------------------------

   if( _openMMBrookInterface.getTriggerEnergyKernel() == this ){
      return (double) _openMMBrookInterface.computeEnergy( context, _system );
   } else {
      return 0.0;
   }

}
