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
#include <sstream>
#include "BrookCalcHarmonicAngleForceKernel.h"

using namespace OpenMM;
using namespace std;

const std::string BrookCalcHarmonicAngleForceKernel::BondName = "HarmonicAngle";

/** 
 * BrookCalcHarmonicAngleForceKernel constructor
 * 
 * @param name                      kernel name
 * @param platform                  platform
 * @param OpenMMBrookInterface      OpenMM-Brook interface
 * @param System                    System reference
 *
 */

BrookCalcHarmonicAngleForceKernel::BrookCalcHarmonicAngleForceKernel( std::string name, const Platform& platform,
                                                                      OpenMMBrookInterface& openMMBrookInterface, System& system ) :
                     CalcHarmonicAngleForceKernel( name, platform ), _openMMBrookInterface( openMMBrookInterface ), _system( system ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcHarmonicAngleForceKernel::BrookCalcHarmonicAngleForceKernel";
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
 * BrookCalcHarmonicAngleForceKernel destructor
 * 
 */

BrookCalcHarmonicAngleForceKernel::~BrookCalcHarmonicAngleForceKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcHarmonicAngleForceKernel::BrookCalcHarmonicAngleForceKernel";
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

FILE* BrookCalcHarmonicAngleForceKernel::getLog( void ) const {
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

int BrookCalcHarmonicAngleForceKernel::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Initialize the kernel, setting up the values of all the force field parameters.
 * 
 * @param system                    System reference
 * @param force                     HarmonicAngleForce reference
 *
 */

void BrookCalcHarmonicAngleForceKernel::initialize( const System& system, const HarmonicAngleForce& force ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcHarmonicAngleForceKernel::initialize";

// ---------------------------------------------------------------------------------------

   FILE* log                 = getLog();

   // ---------------------------------------------------------------------------------------

   // create _brookBondParameters object containing particle indices/parameters

   int numberOfBonds         = force.getNumAngles();

   if( _brookBondParameters ){
      delete _brookBondParameters;
   }
   _brookBondParameters = new BrookBondParameters( BrookCalcHarmonicAngleForceKernel::BondName, NumberOfParticlesInBond, NumberOfParametersInBond, numberOfBonds, getLog() );

   for( int ii = 0; ii < numberOfBonds; ii++ ){

      int particle1, particle2, particle3;
      double angle, k;

      int particles[NumberOfParticlesInBond];
      double parameters[NumberOfParametersInBond];

      force.getAngleParameters( ii, particle1, particle2, particle3, angle, k ); 
      particles[0]    = particle1;
      particles[1]    = particle2;
      particles[2]    = particle3;
 
      parameters[0]   = angle;
      parameters[1]   = k;

      _brookBondParameters->setBond( ii, particles, parameters );
   }   
   _openMMBrookInterface.setHarmonicAngleForceParameters( _brookBondParameters );
   _openMMBrookInterface.setTriggerForceKernel( this );
   _openMMBrookInterface.setTriggerEnergyKernel( this );

   if( log ){
      std::string contents = _brookBondParameters->getContentsString( ); 
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

void BrookCalcHarmonicAngleForceKernel::executeForces( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName   = "BrookCalcHarmonicAngleForceKernel::executeForces";

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

double BrookCalcHarmonicAngleForceKernel::executeEnergy( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcHarmonicAngleForceKernel::executeEnergy";

// ---------------------------------------------------------------------------------------

   if( _openMMBrookInterface.getTriggerEnergyKernel() == this ){
      return (double) _openMMBrookInterface.computeEnergy( context, _system );
   } else {
      return 0.0;
   }

}
