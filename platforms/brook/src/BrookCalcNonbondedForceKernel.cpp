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
#include "BrookCalcNonbondedForceKernel.h"
#include "NonbondedForce.h"

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

   const BrookPlatform brookPlatform        = dynamic_cast<const BrookPlatform&> (platform);
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

void BrookCalcNonbondedForceKernel::initialize( const System& system, const NonbondedForce& force, const std::vector<std::set<int> >& exclusions ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcNonbondedForceKernel::initialize";

// ---------------------------------------------------------------------------------------

    FILE* log                 = getLog();
(void) fprintf( log, "%s begin\n", methodName.c_str() ); fflush( log );

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

   initialize14Interactions( system, force );

}

/** 
 * Initialize the kernel, setting up the values of all the force field parameters.
 * 
 * @param system                    System reference
 * @param force                     HarmonicLJ14Force reference
 *
 */

void BrookCalcNonbondedForceKernel::initialize14Interactions( const System& system, const NonbondedForce& force ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcNonbondedForceKernel::initialize14Interactions";

// ---------------------------------------------------------------------------------------

   FILE* log                 = getLog();

//(void) fprintf( log, "%s begin\n", methodName.c_str() ); fflush( log );

   // ---------------------------------------------------------------------------------------

   // create _brookBondParameters object containing particle indices/parameters

   int numberOf14Forces         = force.getNumNonbonded14();
   if( numberOf14Forces > 0 ){

      _brookBondParameters         = new BrookBondParameters( BondName, NumberOfParticlesInBond, NumberOfParametersInBond, numberOf14Forces, getLog() );
   
      for( int ii = 0; ii < numberOf14Forces; ii++ ){
   
         int particle1, particle2;
         double  charge, radius, depth;
   
         int particles[NumberOfParticlesInBond];
         double parameters[NumberOfParametersInBond];
   
         force.getNonbonded14Parameters( ii, particle1, particle2, charge, radius, depth ); 
   
(void) fprintf( log, "%s idx=%d [%d %d] [%f %f %f]\n", methodName.c_str(), ii, particle1, particle2, charge, radius, depth );
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
