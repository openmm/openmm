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
#include "openmm/OpenMMException.h"
#include <sstream>

#include "openmm/LangevinIntegrator.h"
#include "ReferencePlatform.h"
#include "openmm/internal/OpenMMContextImpl.h"

#include "BrookStreamImpl.h"
#include "OpenMMBrookInterface.h"
#include "kernels/kcommon.h"

using namespace OpenMM;
using namespace std;

/** 
 * OpenMMBrookInterface constructor
 * 
 */

OpenMMBrookInterface::OpenMMBrookInterface( int streamWidth, int duplicationFactor ) : _particleStreamWidth(streamWidth){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "OpenMMBrookInterface::OpenMMBrookInterface";

// ---------------------------------------------------------------------------------------

   _numberOfParticles                       = 0;

   _triggerForceKernel                      = NULL;
   _triggerEnergyKernel                     = NULL;

   _positions                               = NULL;
   _velocities                              = NULL;
   _forces                                  = NULL;

   _log                                     = NULL;
//_log = stderr;

   _particleStreamSize                      = -1;

   for( int ii = 0; ii < LastBondForce; ii++ ){
      _bondParameters[ii] = NULL;
   }

   if( duplicationFactor < 1 ){
      duplicationFactor = 4;
   }
   _brookNonBonded.setDuplicationFactor( duplicationFactor );
   _brookGbsa.setDuplicationFactor( duplicationFactor );
}   

/** 
 * OpenMMBrookInterface destructor
 * 
 */

OpenMMBrookInterface::~OpenMMBrookInterface( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "OpenMMBrookInterface::OpenMMBrookInterface";

// ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < LastBondForce; ii++ ){
      delete _bondParameters[ii];
   }

}

/** 
 * Get number of particles
 * 
 * @return   number of particles
 *
 */

int OpenMMBrookInterface::getNumberOfParticles( void ) const {
   return _numberOfParticles;
}

/** 
 * Set number of particles
 * 
 * @param numberOfParticles number of particles
 *
 */

int OpenMMBrookInterface::setNumberOfParticles( int numberOfParticles ){
   _numberOfParticles = numberOfParticles;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Get particle stream width
 * 
 * @return    particle stream width
 *
 */

int OpenMMBrookInterface::getParticleStreamWidth( void ) const {
   return _particleStreamWidth;
}

/** 
 * Get particle stream size
 * 
 * @return    particle stream size
 *
 */

int OpenMMBrookInterface::getParticleStreamSize( void ) const {
   if( _particleStreamSize < 0 && _numberOfParticles > 0 ){
       OpenMMBrookInterface* localThis = const_cast<OpenMMBrookInterface* const>(this);
      localThis->_particleStreamSize = BrookPlatform::getStreamSize( _numberOfParticles, _particleStreamWidth, NULL );
   }
   return _particleStreamSize;
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
   _log          = log;
   _brookBonded.setLog( log );
   _brookNonBonded.setLog( log );
   _brookGbsa.setLog( log );
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Get BrookBondParameters at specified index
 * 
 * @param   index
 *
 * @return  BrookBondParameters* object
 *
 */

BrookBondParameters* OpenMMBrookInterface::_getBondParameters( BondParameterIndices index ) const {
   return _bondParameters[index];
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

   if( brookBondParameters && brookBondParameters->getNumberOfBonds() > 0 ){
      _brookBonded.setIsActive( 1 );
   }

   _bondParameters[index] = brookBondParameters;

   if( !brookBondParameters || brookBondParameters->getNumberOfBonds() < 1 ){
      int isActive = 0;
      for( int ii = 0; ii < LastBondForce && !isActive; ii++ ){
         isActive = (_bondParameters[ii] != NULL && brookBondParameters->getNumberOfBonds() > 0) ? 1 : 0;
      }
       _brookBonded.setIsActive( isActive );
   }
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Get BrookBondParameters for harmonic bond force
 * 
 * @return brookBondParameters for harmonic bond force
 *
 */

BrookBondParameters* OpenMMBrookInterface::getHarmonicBondForceParameters( void ) const {
   return _getBondParameters( HarmonicBondIndex );
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
   _brookBonded.setupCompleted( 0 );
   return _setBondParameters( HarmonicBondIndex, brookBondParameters );
}

/** 
 * Get BrookBondParameters for harmonic angle force
 * 
 * @return brookBondParameters for harmonic angle force
 *
 */

BrookBondParameters* OpenMMBrookInterface::getHarmonicAngleForceParameters( void ) const {
   return _getBondParameters( HarmonicAngleIndex );
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
   _brookBonded.setupCompleted( 0 );
   return _setBondParameters( HarmonicAngleIndex, brookBondParameters );
}

/** 
 * Get BrookBondParameters for periodic torsion force
 * 
 * @return brookBondParameters for periodic torsion force
 *
 */

BrookBondParameters* OpenMMBrookInterface::getPeriodicTorsionForceParameters( void ) const {
   return _getBondParameters( PeriodicTorsionForceIndex );
}

/** 
 * Set BrookBondParameters for periodic torsion force
 * 
 * @param brookBondParameters brookBondParameters for periodic torsion force
 *
 * @return  DefaultReturnValue
 *
 */

int OpenMMBrookInterface::setPeriodicTorsionForceParameters( BrookBondParameters* brookBondParameters ){
   _brookBonded.setupCompleted( 0 );
   return _setBondParameters( PeriodicTorsionForceIndex, brookBondParameters );
}

/** 
 * Get BrookBondParameters for rb torsion force
 * 
 * @return brookBondParameters for rb torsion force
 *
 */

BrookBondParameters* OpenMMBrookInterface::getRBTorsionForceParameters( void ) const {
   return _getBondParameters( RbTorsionForceIndex );
}

/** 
 * Set BrookBondParameters for RB torsion force
 * 
 * @param brookBondParameters brookBondParameters for RB torsion force
 *
 * @return  DefaultReturnValue
 *
 */

int OpenMMBrookInterface::setRBTorsionForceParameters( BrookBondParameters* brookBondParameters ){
   _brookBonded.setupCompleted( 0 );
   return _setBondParameters( RbTorsionForceIndex, brookBondParameters );
}

/** 
 * Get BrookBondParameters for LJ 14 force
 * 
 * @return brookBondParameters for LJ 14 force
 *
 */

BrookBondParameters* OpenMMBrookInterface::getNonBonded14ForceParameters( void ) const {
   return _getBondParameters( LJ14Index );
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
   _brookBonded.setupCompleted( 0 );
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
 * Set trigger Force Kernel
 *
 * @param triggerForceKernel kernel to calculate force
 *
 */
    
void OpenMMBrookInterface::setTriggerForceKernel( void* triggerForceKernel ){
   _triggerForceKernel = triggerForceKernel;
}
    
/** 
 * Get trigger Force Kernel
 *
 * @return triggerForceKernel kernel to calculate force
 *
 */
    
void* OpenMMBrookInterface::getTriggerForceKernel( void ) const {
   return _triggerForceKernel;
}

/** 
 * Set trigger Energy Kernel
 *
 * @param triggerEnergyKernel kernel to calculate force
 *
 */
    
void OpenMMBrookInterface::setTriggerEnergyKernel( void* triggerEnergyKernel ){
   _triggerEnergyKernel = triggerEnergyKernel;
}
    
/** 
 * Get trigger Energy Kernel
 *
 * @return triggerEnergyKernel kernel to calculate force
 *
 */
    
void* OpenMMBrookInterface::getTriggerEnergyKernel( void ) const {
   return _triggerEnergyKernel;
}

/** 
 * Get Brook non bonded
 *
 * @return BrookNonBonded reference
 *
 */
    
BrookNonBonded& OpenMMBrookInterface::getBrookNonBonded( void ){
   return _brookNonBonded;
}

/** 
 * Get Brook GBSA
 *
 * @return BrookGbsa reference
 *
 */
    
BrookGbsa& OpenMMBrookInterface::getBrookGbsa( void ){
   return _brookGbsa;
}

/** 
 * Zero forces
 * 
 * @param context     context
 *
 */

void OpenMMBrookInterface::zeroForces( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "OpenMMBrookInterface::zeroForces";

// ---------------------------------------------------------------------------------------

   // zero forces

   BrookStreamImpl* forces  = getParticleForces();
   kzerof3( forces->getBrookStream() );

   // ---------------------------------------------------------------------------------------
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

   static       int printOn                 = 0;
   static const int MaxErrorMessages        = 2;
   static       int ErrorMessages           = 0;
   FILE* log;

// ---------------------------------------------------------------------------------------

//setLog( stderr );

   if( printOn && getLog() ){
      log     = getLog();
   } else {
      printOn = 0;
   }

   BrookStreamImpl* positions = getParticlePositions();
   BrookStreamImpl* forces    = getParticleForces();

   // info

   if( printOn > 1 ){
   //if( 1 ){
      printForcesToFile( context );
   }

   // nonbonded forces

   if( _brookNonBonded.isActive() ){
      _brookNonBonded.computeForces( *positions, *forces );
   }

// ---------------------------------------------------------------------------------------

   // bonded forces

   if( printOn ){
      (void) fprintf( log, "%s done nonbonded: bonded=%d completed=%d\n", methodName.c_str(),
                      _brookBonded.isActive(), _brookBonded.isSetupCompleted() ); 
      (void) fflush( log );
   }

   if( _brookBonded.isActive() ){

      // perform setup first time through

      if( _brookBonded.isSetupCompleted() == 0 ){
         _brookBonded.setup( getNumberOfParticles(), getHarmonicBondForceParameters(), getHarmonicAngleForceParameters(),
                             getPeriodicTorsionForceParameters(), getRBTorsionForceParameters(),
                             getNonBonded14ForceParameters(),
                             getParticleStreamWidth(), getParticleStreamSize() );
      }

      _brookBonded.computeForces( *positions, *forces );
   
      // diagnostics
   
      if( printOn ){
         static int step           = 0;
         static const int stopStep = 10;
         (void) fprintf( log, "%s done bonded computeForces\n", methodName.c_str() );
         (void) fflush( log );
   
   /*
         (void) fprintf( getLog(), "\nFinal NB & bonded forces" );
         BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamInternal();
         brookStreamInternalF->printToFile( getLog() );
         void* dataV = brookStreamInternalF->getData(1);
         float* data = (float*) dataV;
         (void) fprintf( getLog(), "\nFinal NB & bonded forces RAW\n" );
         for( int ii = 0; ii < _brookNonBonded->getNumberOfParticles()*3; ii += 3 ){
            (void) fprintf( getLog(), "%d [%.6e %.6e %.6e]\n", ii, data[ii], data[ii+1], data[ii+2] );
         }
   */
   
         if( step++ >= stopStep ){
            exit(0);
         }
      }

   }

   // GBSA OBC forces

   if( _brookGbsa.isActive() ){
      _brookGbsa.computeForces( *positions, *forces );
   }

   // ---------------------------------------------------------------------------------------
}

/** 
 * Print forces to file 
 * 
 * @param context     context
 *
 */

void OpenMMBrookInterface::printForcesToFile( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "OpenMMBrookInterface::printForcesToFile";
   float zero                               = 0.0f;
   static int step                          = 0;

// ---------------------------------------------------------------------------------------

   // first step only?

   if( step > 1 ){
      return;
   }

   std::stringstream fileNameBaseS;
   //fileNameBase << "Brook_" << context.getTime() << "_"; 
   fileNameBaseS << "Brook_" << step++ << "_"; 
   std::string fileNameBase = fileNameBaseS.str();

   // create vector of streams for output

   BrookStreamImpl* positions = getParticlePositions();
   BrookStreamImpl* forces    = getParticleForces();

   std::vector<BrookStreamInternal*> streams;
   streams.push_back( positions->getBrookStreamInternal() );
   streams.push_back( forces->getBrookStreamInternal() );

// ---------------------------------------------------------------------------------------

   // nonbonded forces

   if( _brookNonBonded.isActive() ){
      forces->fillWithValue( &zero );
      _brookNonBonded.computeForces( *positions, *forces );
      std::string fileName = fileNameBase + "NonBonded.txt";
      BrookStreamInternal::printStreamsToFile( fileName, streams );
   }

// ---------------------------------------------------------------------------------------

   // bonded forces

   if( _brookBonded.isActive() ){

      // perform setup first time through

      if( _brookBonded.isSetupCompleted() == 0 ){

         _brookBonded.setup( getNumberOfParticles(), getHarmonicBondForceParameters(), getHarmonicAngleForceParameters(),
                             getPeriodicTorsionForceParameters(), getRBTorsionForceParameters(),
                             getNonBonded14ForceParameters(),
                             getParticleStreamWidth(), getParticleStreamSize() );
      }

      forces->fillWithValue( &zero );
      _brookBonded.computeForces( *positions, *forces );
      std::string fileName = fileNameBase + "Bonded.txt";
      BrookStreamInternal::printStreamsToFile( fileName, streams );
   
   }

// ---------------------------------------------------------------------------------------

   // GBSA OBC forces

   if( _brookGbsa.isActive() ){

      forces->fillWithValue( &zero );
      _brookGbsa.computeForces( *positions, *forces );

      std::string fileName = fileNameBase + "Obc.txt";
      BrookStreamInternal::printStreamsToFile( fileName, streams );

   }

// ---------------------------------------------------------------------------------------

   // all forces

   if( 1 ){

      forces->fillWithValue( &zero );

      if( _brookNonBonded.isActive() ){
         _brookNonBonded.computeForces( *positions, *forces );
      }

      if( _brookBonded.isActive() ){
         _brookBonded.computeForces( *positions, *forces );
      }

      if( _brookGbsa.isActive() ){
         _brookGbsa.computeForces( *positions, *forces );
      }

      std::string fileName = fileNameBase + "AllF.txt";
      BrookStreamInternal::printStreamsToFile( fileName, streams );
   }

   forces->fillWithValue( &zero );

   // ---------------------------------------------------------------------------------------
}

/** 
 * Compute energy
 * 
 * @param context     context
 * @param system      system reference
 *
 */

double OpenMMBrookInterface::computeEnergy( OpenMMContextImpl& context, System& system ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "OpenMMBrookInterface::computeEnergy";

// ---------------------------------------------------------------------------------------

   // We don't currently have GPU kernels to calculate energy, so instead we have the reference
   // platform do it.  This is VERY slow.
    
   LangevinIntegrator integrator(0.0, 1.0, 0.0);
   ReferencePlatform platform;
   OpenMMContext refContext( system, integrator, platform );

   const Stream& positions = context.getPositions();
   double* posData         = new double[positions.getSize()*3];
   positions.saveToArray(posData);
   vector<Vec3> pos(positions.getSize());

   for( unsigned int ii = 0; ii < pos.size(); ii++ ){
      pos[ii] = Vec3(posData[3*ii], posData[3*ii+1], posData[3*ii+2]);
   }
   delete[] posData;
   refContext.setPositions(pos);

   return refContext.getState(State::Energy).getPotentialEnergy();

}
