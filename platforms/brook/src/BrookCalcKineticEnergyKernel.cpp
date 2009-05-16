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

#include <sstream>
#include "openmm/OpenMMException.h"
#include "BrookCalcKineticEnergyKernel.h"
#include "BrookStreamImpl.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookCalcKineticEnergyKernel constructor
 * 
 * @param name                      name of the stream to create
 * @param platform                  platform
 * @param OpenMMBrookInterface      OpenMM-Brook interface
 * @param System                    System reference
 *
 */

BrookCalcKineticEnergyKernel::BrookCalcKineticEnergyKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface, System& system ) :
                              CalcKineticEnergyKernel( name, platform ), _openMMBrookInterface( openMMBrookInterface ), _system( system ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcKineticEnergyKernel::BrookCalcKineticEnergyKernel";

// ---------------------------------------------------------------------------------------

   _openMMBrookInterface.setNumberOfParticles( system.getNumParticles() );
   _numberOfParticles  =  system.getNumParticles();
   _masses             = NULL;

}

/** 
 * BrookCalcKineticEnergyKernel destructor
 * 
 */
  
BrookCalcKineticEnergyKernel::~BrookCalcKineticEnergyKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcKineticEnergyKernel::~BrookCalcKineticEnergyKernel";

// ---------------------------------------------------------------------------------------

   delete[] _masses;
}

/** 
 * Initialize the kernel
 * 
 * @param system  System reference
 *
 */

void BrookCalcKineticEnergyKernel::initialize( const System& system ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcKineticEnergyKernel::initialize";

// ---------------------------------------------------------------------------------------

   _numberOfParticles  = system.getNumParticles();

   // load masses

   if( _masses ){
      delete[] _masses;
   }

   _masses = new BrookOpenMMFloat[_numberOfParticles];

   for( unsigned int ii = 0; ii < (unsigned int) _numberOfParticles; ii++ ){
      _masses[ii]  =  static_cast<BrookOpenMMFloat>(system.getParticleMass(ii));
   }

/*
std::vector<double> masses;
for( unsigned int ii = 0; ii < (unsigned int) _numberOfParticles; ii++ ){
   masses.push_back( system.getParticleMass(ii) );
}
_brookVelocityCenterOfMassRemoval = new BrookVelocityCenterOfMassRemoval();
_brookVelocityCenterOfMassRemoval->setup( masses, getPlatform() );
*/

   return;
}

/** 
 * Calculate kinetic energy
 * 
 * @param context OpenMMContextImpl reference
 *
 * @return kinetic energy of the system
 *
 */

double BrookCalcKineticEnergyKernel::execute( OpenMMContextImpl& context ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName    = "BrookCalcKineticEnergyKernel::execute";

// ---------------------------------------------------------------------------------------

//BrookStreamImpl* velocities = _openMMBrookInterface.getParticleVelocities();
//_brookVelocityCenterOfMassRemoval->removeVelocityCenterOfMass( *velocities );

   void* dataV                            = _openMMBrookInterface.getParticleVelocities()->getData( 1 );
   float* velocity                        = (float*) dataV;

if( 0 && _numberOfParticles ){
printf( "BrookCalcKineticEnergyKernel:\n" );
float com[3]      = { 0.0, 0.0, 0.0 };
float localEnergy = 0.0f;
int localIndex    = 0;
float massSum     = 0.0f;
for ( int ii = 0; ii < _numberOfParticles; ii++, localIndex += 3 ){

   com[0]      += _masses[ii]*velocity[localIndex];
   com[1]      += _masses[ii]*velocity[localIndex+1];
   com[2]      += _masses[ii]*velocity[localIndex+2];
   localEnergy += _masses[ii]*(velocity[localIndex]*velocity[localIndex] + velocity[localIndex + 1]*velocity[localIndex + 1] + velocity[localIndex + 2]*velocity[localIndex + 2]);

   massSum  += _masses[ii];

   printf( "  %d %.3f [%12.5e %12.5e %12.5e]\n", ii, _masses[ii], velocity[localIndex], velocity[localIndex+1], velocity[localIndex+2] );
}
float inverseTotalMass = 1.0f/massSum;
com[0] *= inverseTotalMass;
com[1] *= inverseTotalMass;
com[2] *= inverseTotalMass;
printf( "KE raw=%.5e Com [%12.5e %12.5e %12.5e]\n", 0.5f*localEnergy, com[0], com[1], com[2] );


float newcom[3] = { 0.0, 0.0, 0.0 };
localIndex      = 0;
for ( int ii = 0; ii < _numberOfParticles; ii++, localIndex += 3 ){
   velocity[localIndex]   -= com[0];
   velocity[localIndex+1] -= com[1];
   velocity[localIndex+2] -= com[2];
   newcom[0]              += velocity[localIndex];
   newcom[1]              += velocity[localIndex+1];
   newcom[2]              += velocity[localIndex+2];
   printf( "  %d %.3f [%12.5e %12.5e %12.5e]\n", ii, _masses[ii], velocity[localIndex], velocity[localIndex+1], velocity[localIndex+2] );
}
printf( "NewCom [%12.5e %12.5e %12.5e]\n", newcom[0], newcom[1], newcom[2] );

}

   int index                              = 0;
   double energy                          = 0.0;
   for( int ii = 0; ii < _numberOfParticles; ii++, index += 3 ){
      energy += _masses[ii]*(velocity[index]*velocity[index] + velocity[index + 1]*velocity[index + 1] + velocity[index + 2]*velocity[index + 2]);
   }

//printf( " KQ=%12.5e\n", 0.5*energy );
   return 0.5*energy;
}
