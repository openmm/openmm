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
#include "BrookCalcHarmonicBondForceKernel.h"
#include <sstream>

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

   _brookBondParameters               = NULL;
   _log                               = NULL;

   const BrookPlatform& brookPlatform = dynamic_cast<const BrookPlatform&> (platform);
   if( brookPlatform.getLog() != NULL ){
      setLog( brookPlatform.getLog() );
   }

   _openMMBrookInterface.setNumberOfParticles( system.getNumParticles() );
      
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
   int printOn                              = 0;
   FILE* log;

// ---------------------------------------------------------------------------------------

   if( printOn && getLog() ){
       log = getLog();
   } else {
      printOn = 0;
   }   

   // ---------------------------------------------------------------------------------------

   // create _brookBondParameters object containing particle indices/parameters

   int numberOfBonds         = force.getNumBonds();

   if( _brookBondParameters ){
      delete _brookBondParameters;
   }
   _brookBondParameters = new BrookBondParameters( BrookCalcHarmonicBondForceKernel::BondName, NumberOfParticlesInBond, NumberOfParametersInBond, numberOfBonds, getLog() );

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

   if( printOn ){
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

void BrookCalcHarmonicBondForceKernel::executeForces( ContextImpl& context ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName   = "BrookCalcHarmonicBondForceKernel::executeForces";

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

double BrookCalcHarmonicBondForceKernel::executeEnergy( ContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcHarmonicBondForceKernel::executeEnergy";

// ---------------------------------------------------------------------------------------

   if( _openMMBrookInterface.getTriggerEnergyKernel() == this ){
      return (double) _openMMBrookInterface.computeEnergy( context, _system );
   } else {
      return 0.0;
   }

}
