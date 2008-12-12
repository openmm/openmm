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
#include "OpenMMBrookInterface.h"
#include "gpu/kforce.h"
#include "gpu/kinvmap_gather.h"
#include "NonbondedForce.h"

using namespace OpenMM;
using namespace std;

/** 
 * OpenMMBrookInterface constructor
 * 
 */

OpenMMBrookInterface::OpenMMBrookInterface( void ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "OpenMMBrookInterface::OpenMMBrookInterface";

// ---------------------------------------------------------------------------------------

   _numberOfParticles                           = 0;

   _brookBonded                             = NULL;
   _brookNonBonded                          = NULL;
   _brookGbsa                               = NULL;

   _positions                               = NULL;
   _velocities                              = NULL;
   _forces                                  = NULL;

   _refForceField                           = NULL;
   _log                                     = NULL;

   _refForceField                           = NULL;
   _refSystem                               = NULL;
   _refOpenMMContext                        = NULL;
   _referencePlatform                       = NULL;
   _refVerletIntegrator                     = NULL;

}   

/** 
 * OpenMMBrookInterface destructor
 * 
 */

OpenMMBrookInterface::~OpenMMBrookInterface( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "OpenMMBrookInterface::OpenMMBrookInterface";

// ---------------------------------------------------------------------------------------

   delete _brookBonded;
   delete _brookNonBonded;

   // deleted w/ kernel delete? If activated, program crashes

   //delete _refForceField;
/*
   delete _refSystem;
   delete _refOpenMMContext;
   delete _referencePlatform;
   delete _refVerletIntegrator;
*/
}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* OpenMMBrookInterface::getLog( void ) const {
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

int OpenMMBrookInterface::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Set BrookBondParameters at specified index
 * 
 * @param index
 * @param brookBondParameters brookBondParameters for BrookBondParameters
 *
 * @return  DefaultReturnValue
 *
 */

int OpenMMBrookInterface::_setBondParameters( BondParameterIndices index, BrookBondParameters* brookBondParameters ){
   _bondParameters[index] = brookBondParameters;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Set BrookBondParameters for harmonic bond force
 * 
 * @param brookBondParameters brookBondParameters for harmonic bond force
 *
 * @return  DefaultReturnValue
 *
 */

int OpenMMBrookInterface::setHarmonicBondForceParameters( BrookBondParameters* brookBondParameters ){
   return _setBondParameters( HarmonicAngleIndex, brookBondParameters );
}

/** 
 * Set BrookBondParameters for harmonic angle force
 * 
 * @param brookBondParameters brookBondParameters for harmonic angle force
 *
 * @return  DefaultReturnValue
 *
 */

int OpenMMBrookInterface::setHarmonicAngleForceParameters( BrookBondParameters* brookBondParameters ){
   return _setBondParameters( HarmonicAngleIndex, brookBondParameters );
}

/** 
 * Set BrookBondParameters for proper dihedral force
 * 
 * @param brookBondParameters brookBondParameters for proper dihedral force
 *
 * @return  DefaultReturnValue
 *
 */

int OpenMMBrookInterface::setPeriodicTorsionForceParameters( BrookBondParameters* brookBondParameters ){
   return _setBondParameters( PeriodicTorsionForceIndex, brookBondParameters );
}

/** 
 * Set BrookBondParameters for RB dihedral force
 * 
 * @param brookBondParameters brookBondParameters for RB force
 *
 * @return  DefaultReturnValue
 *
 */

int OpenMMBrookInterface::setRBTorsionForceParameters( BrookBondParameters* brookBondParameters ){
   return _setBondParameters( RbTorsionForceIndex, brookBondParameters );
}

/** 
 * Set BrookBondParameters for LJ 14 force
 * 
 * @param brookBondParameters brookBondParameters for LJ 14 force
 *
 * @return  DefaultReturnValue
 *
 */

int OpenMMBrookInterface::setNonBonded14ForceParameters( BrookBondParameters* brookBondParameters ){
   return _setBondParameters( LJ14Index, brookBondParameters );
}

/** 
 * Get positions stream
 * 
 * @return particle positions 
 *
 */
    
BrookStreamImpl* OpenMMBrookInterface::getParticlePositions( void ){
   return _positions;
}

/** 
 * Set positions stream
 * 
 * @param positions Brook stream containing particle positions
 *
 * @return  DefaultReturnValue
 *
 */
    
int OpenMMBrookInterface::setParticlePositions( BrookStreamImpl* positions ){
   _positions = positions;
   return DefaultReturnValue;
}

/** 
 * Get velocities stream
 * 
 * @return particle velocities
 *
 */
    
BrookStreamImpl* OpenMMBrookInterface::getParticleVelocities( void ){
   return _velocities;
}

/** 
 * Set velocities stream
 * 
 * @param velocities Brook stream containing particle velocities
 *
 * @return  DefaultReturnValue
 *
 */
    
int OpenMMBrookInterface::setParticleVelocities( BrookStreamImpl* velocities ){
   _velocities = velocities;
   return DefaultReturnValue;
}

/** 
 * Get forces stream
 * 
 * @return ParticleForces
 *
 */
    
BrookStreamImpl* OpenMMBrookInterface::getParticleForces( void ){
   return _forces;
}

/** 
 * Set forces stream
 * 
 * @param forces Brook stream containing particle forces
 *
 * @return  DefaultReturnValue
 *
 */
    
int OpenMMBrookInterface::setParticleForces( BrookStreamImpl* forces ){
   _forces = forces;
   return DefaultReturnValue;
}

/** 
 * Compute forces
 * 
 * @param context     context
 *
 */

void OpenMMBrookInterface::computeForces( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "OpenMMBrookInterface::executeForces";

   static const int PrintOn                 = 0;
   static const int MaxErrorMessages        = 2;
   static       int ErrorMessages           = 0;

   static const float4 dummyParameters( 0.0, 0.0, 0.0, 0.0 );

   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // nonbonded forces

   BrookStreamImpl* positions = getParticlePositions();
   BrookStreamImpl* forces    = getParticleForces();

   if( _brookNonBonded ){
      _brookNonBonded->computeForces( *positions, *forces );
   }

// ---------------------------------------------------------------------------------------

   // bonded forces

   if( _brookBonded ){

      _brookBonded->computeForces( *positions, *forces );
   
      // diagnostics
   
      if( 1 && PrintOn ){
   
   /*
         (void) fprintf( getLog(), "\nFinal NB & bonded forces" );
         BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamImpl();
         brookStreamInternalF->printToFile( getLog() );
         void* dataV = brookStreamInternalF->getData(1);
         float* data = (float*) dataV;
         (void) fprintf( getLog(), "\nFinal NB & bonded forces RAW\n" );
         for( int ii = 0; ii < _brookNonBonded->getNumberOfParticles()*3; ii += 3 ){
            (void) fprintf( getLog(), "%d [%.6e %.6e %.6e]\n", ii, data[ii], data[ii+1], data[ii+2] );
         }
   */
   
      }
   }

   // GBSA OBC forces

   if( _brookGbsa ){
      _brookGbsa->computeForces( *positions, *forces );
   }

   // ---------------------------------------------------------------------------------------
}
