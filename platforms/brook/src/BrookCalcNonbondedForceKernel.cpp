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
#include "BrookStreamImpl.h"
#include "BrookCalcNonbondedForceKernel.h"
#include "openmm/NonbondedForce.h"

#include <cmath>
#include <limits>
#include <sstream>

using namespace OpenMM;
using namespace std;

const std::string BrookCalcNonbondedForceKernel::BondName = "LJ14";

/** 
 * BrookCalcNonbondedForceKernel constructor
 * 
 * @param name                      kernel name
 * @param platform                  platform
 * @param openMMBrookInterface      OpenMMBrookInterface reference
 * @param system                    System reference 
 *
 */

BrookCalcNonbondedForceKernel::BrookCalcNonbondedForceKernel( std::string name, const Platform& platform,
                                                              OpenMMBrookInterface& openMMBrookInterface, System& system ) :
                     CalcNonbondedForceKernel( name, platform ), _openMMBrookInterface( openMMBrookInterface ), _system( system ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcNonbondedForceKernel::BrookCalcNonbondedForceKernel";

// ---------------------------------------------------------------------------------------

   _numberOfParticles                       = system.getNumParticles();
   _openMMBrookInterface.setNumberOfParticles( system.getNumParticles() );

   _brookBondParameters                     = NULL;
   _log                                     = NULL;

   const BrookPlatform& brookPlatform       = dynamic_cast<const BrookPlatform&> (platform);
   if( brookPlatform.getLog() != NULL ){
      setLog( brookPlatform.getLog() );
   }
      
}   

/** 
 * BrookCalcNonbondedForceKernel destructor
 * 
 */

BrookCalcNonbondedForceKernel::~BrookCalcNonbondedForceKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcNonbondedForceKernel::BrookCalcNonbondedForceKernel";

// ---------------------------------------------------------------------------------------

   delete _brookBondParameters;

}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookCalcNonbondedForceKernel::getLog( void ) const {
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

int BrookCalcNonbondedForceKernel::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Initialize object 
 * 
 * @param  system     System reference (currently not used)
 * @param  force      NonbondedForce reference -- extract charge, and vdw parameters from this object
 * @param  exclusions list of execlusions
 *
 * @return  DefaultReturnValue
 *
 */

void BrookCalcNonbondedForceKernel::initialize( const System& system, const NonbondedForce& force ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcNonbondedForceKernel::initialize";
   static const int PrintOn                 = 0;

// ---------------------------------------------------------------------------------------

    FILE* log                 = getLog();
    if( PrintOn && log ){
       (void) fprintf( log, "%s begin\n", methodName.c_str() ); fflush( log );
    } else {
       log = NULL;
    }

    _numberOfParticles        = force.getNumParticles();
/*
    nonbondedMethod = CalcNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = (RealOpenMM) force.getCutoffDistance();
    Vec3 boxVectors[3];
    force.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
    periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
    periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];
    if (nonbondedMethod == NoCutoff)
        neighborList = NULL;
    else
        neighborList = new NeighborList();
*/
   // ---------------------------------------------------------------------------------------

   // nonbonded

   BrookNonBonded& brookNonBonded = _openMMBrookInterface.getBrookNonBonded();

   // charge & LJ parameters

   std::vector<std::vector<double> > nonbondedParameters;
   nonbondedParameters.resize( _numberOfParticles );
   for( int ii = 0; ii < _numberOfParticles; ii++ ){
      double charge, radius, depth;
      force.getParticleParameters( ii, charge, radius, depth );
      nonbondedParameters[ii].push_back( charge  );
      nonbondedParameters[ii].push_back( radius  );
      nonbondedParameters[ii].push_back( depth   );
   }

   // Go through the exclusions.

   std::vector<std::set<int> > exclusions(_numberOfParticles);
   std::vector<int> nb14s;
   for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
        if (chargeProd != 0.0 || epsilon != 0.0)
            nb14s.push_back(i);
   }

   brookNonBonded.setup( _numberOfParticles, nonbondedParameters, exclusions, getPlatform() );
   _openMMBrookInterface.setTriggerForceKernel( this );
   _openMMBrookInterface.setTriggerEnergyKernel( this );

   // echo contents

   if( log ){
      std::string contents = brookNonBonded.getContentsString( );
      (void) fprintf( log, "%s brookNonBonded::contents\n%s", methodName.c_str(), contents.c_str() );
      (void) fflush( log );
   }

   // nonbonded 14 ixns

   initialize14Interactions( system, force, nb14s );

}

/** 
 * Initialize the kernel, setting up the values of all the force field parameters.
 * 
 * @param system                    System reference
 * @param force                     HarmonicLJ14Force reference
 * @param nb14s                     which of the exceptions need to be calculated
 */

void BrookCalcNonbondedForceKernel::initialize14Interactions( const System& system, const NonbondedForce& force, const std::vector<int>& nb14s ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcNonbondedForceKernel::initialize14Interactions";

// ---------------------------------------------------------------------------------------

   FILE* log                 = getLog();

//(void) fprintf( log, "%s begin\n", methodName.c_str() ); fflush( log );

   // ---------------------------------------------------------------------------------------

   // create _brookBondParameters object containing particle indices/parameters

   int numberOf14Forces         = nb14s.size();
   if( numberOf14Forces > 0 ){

      _brookBondParameters         = new BrookBondParameters( BondName, NumberOfParticlesInBond, NumberOfParametersInBond, numberOf14Forces, getLog() );
   
      for( int ii = 0; ii < numberOf14Forces; ii++ ){
   
         int particle1, particle2;
         double  charge, radius, depth;
   
         int particles[NumberOfParticlesInBond];
         double parameters[NumberOfParametersInBond];
   
         force.getExceptionParameters( nb14s[ii], particle1, particle2, charge, radius, depth );
   
//(void) fprintf( log, "%s idx=%d [%d %d] [%f %f %f]\n", methodName.c_str(), ii, particle1, particle2, charge, radius, depth );
         particles[0]    = particle1;
         particles[1]    = particle2;
    
         parameters[0]   = charge;
         parameters[1]   = radius;
         parameters[2]   = depth;
   
         _brookBondParameters->setBond( ii, particles, parameters );
      }   
      _openMMBrookInterface.setNonBonded14ForceParameters( _brookBondParameters );
   
      if( log ){
         std::string contents = _brookBondParameters->getContentsString( ); 
         (void) fprintf( log, "%s contents:\n%s", methodName.c_str(), contents.c_str() );
         (void) fflush( log );
      }
   } else if( log ){
      (void) fprintf( log, "%s no 14 ixns\n", methodName.c_str() );
      (void) fflush( log );
   }

   // ---------------------------------------------------------------------------------------
    
}

/** 
 * Execute the kernel to calculate the nonbonded forces
 * 
 * @param context OpenMMContextImpl context
 *
 */

void BrookCalcNonbondedForceKernel::executeForces( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcNonbondedForceKernel::executeForces";

// ---------------------------------------------------------------------------------------

   if( _openMMBrookInterface.getTriggerForceKernel() == this ){
      _openMMBrookInterface.computeForces( context );
   }   

   // ---------------------------------------------------------------------------------------
}

/**
 * Execute the kernel to calculate the energy.
 * 
 * @param context OpenMMContextImpl context
 *
 * @return  potential energy due to the NonbondedForce
 * Currently always return 0.0 since energies not calculated on gpu
 *
 */

double BrookCalcNonbondedForceKernel::executeEnergy( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcNonbondedForceKernel::executeEnergy";

// ---------------------------------------------------------------------------------------

   if( _openMMBrookInterface.getTriggerEnergyKernel() == this ){
      return (double) _openMMBrookInterface.computeEnergy( context, _system );
   } else {
      return 0.0;
   }   

}
