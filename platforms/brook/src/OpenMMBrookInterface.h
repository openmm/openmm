#ifndef OPENMM_BROOK_INTERFACE_H_
#define OPENMM_BROOK_INTERFACE_H_

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

#include "openmm/kernels.h"
#include "../../reference/src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "BrookBondParameters.h"
#include "BrookBonded.h"
#include "BrookNonBonded.h"
#include "BrookGbsa.h"
#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMContext.h"
#include "openmm/System.h"

namespace OpenMM {

/**
 * OpenMM-Brook interface methods
 */
class OpenMMBrookInterface {

   public:
  
      OpenMMBrookInterface( int streamWidth, int duplicationFactor );
  
      ~OpenMMBrookInterface();
  
      /**
       * Get reference Context
       * 
       * @param numberOfParticles  number of particles
       *
       * @return  OpenMMContext
       *
       */
      
      OpenMMContext* getReferenceOpenMMContext( int numberOfParticles );
      
      /** 
       * Set log file reference
       * 
       * @param  log file reference
       *
       * @return DefaultReturnValue
       *
       */
      
      int setLog( FILE* log );

      /* 
       * Get contents of object
       *
       * @param level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContents( int level ) const;

      /** 
       * Get number of particles
       * 
       * @return   number of particles
       *
       */
      
      int getNumberOfParticles( void ) const;
      
      /** 
       * Set number of particles
       * 
       * @param numberOfParticles number of particles
       *
       * @return DefaultReturnValue
       *
       */
      
      int setNumberOfParticles( int numberOfParticles );
      
      /** 
       * Get particle stream width
       * 
       * @return    particle stream width
       *
       */
      
      int getParticleStreamWidth( void ) const;
      
      /** 
       * Get particle stream size
       * 
       * @return    particle stream size
       *
       */
      
      int getParticleStreamSize( void ) const;
      
      /** 
       * Get log file reference
       * 
       * @return  log file reference
       *
       */
      
      FILE* getLog( void ) const;
      
      /** 
       * Zero forces
       *
       * @param context OpenMMContextImpl context 
       *
       */
      
      void zeroForces( OpenMMContextImpl& context );
      
      /** 
       * Compute forces
       *
       * @param context OpenMMContextImpl context 
       *
       */
      
      void computeForces( OpenMMContextImpl& context );
      
      /** 
       * Compute energy
       *
       * @param context OpenMMContextImpl context 
       * @param system  system reference
       *
       * @return potential energy
       *
       */
      
      double computeEnergy( OpenMMContextImpl& context, System& system );
      
      /** 
       * Set trigger Force Kernel
       *
       * @param triggerForceKernel kernel to calculate force
       *
       */
      
      void setTriggerForceKernel( void* triggerForceKernel );
      
      /** 
       * Get trigger Force Kernel
       *
       * @return triggerForceKernel kernel to calculate force
       *
       */
      
      void* getTriggerForceKernel( void ) const;
      
      /** 
       * Set trigger Energy Kernel
       *
       * @param triggerEnergyKernel kernel to calculate force
       *
       */
      
      void setTriggerEnergyKernel( void* triggerEnergyKernel );
      
      /** 
       * Get trigger Energy Kernel
       *
       * @return triggerEnergyKernel kernel to calculate force
       *
       */
      
      void* getTriggerEnergyKernel( void ) const;
      
      /** 
       * Get BrookNonBonded reference
       * 
       * @return BrookNonBonded reference
       *
       */
      
      BrookNonBonded& getBrookNonBonded( void );
      
      /** 
       * Get BrookGbsa reference
       * 
       * @return BrookGbsa reference
       *
       */
      
      BrookGbsa& getBrookGbsa( void );
      
      /** 
       * Get BrookBondParameters for harmonic bond force
       * 
       * @return brookBondParameters for BrookBondParameters for harmonic bond force
       *
       */
      
      BrookBondParameters* getHarmonicBondForceParameters( void ) const;
      
      /** 
       * Set BrookBondParameters for harmonic bond force
       * 
       * @param brookBondParameters BrookBondParameters for harmonic bond force
       *
       * @return DefaultReturnValue
       *
       */
      
      int setHarmonicBondForceParameters( BrookBondParameters* brookBondParameters );
      
      /** 
       * Get BrookBondParameters for harmonic angle force
       * 
       * @return brookBondParameters for BrookBondParameters for harmonic angle force
       *
       */
      
      BrookBondParameters* getHarmonicAngleForceParameters( void ) const;
      
      /** 
       * Set BrookBondParameters for harmonic angle force
       * 
       * @param brookBondParameters brookBondParameters for BrookBondParameters for harmonic angle force
       *
       * @return DefaultReturnValue
       *
       */
      
      int setHarmonicAngleForceParameters( BrookBondParameters* brookBondParameters );
      
      /** 
       * Get BrookBondParameters for periodic torsion force
       * 
       * @return brookBondParameters for BrookBondParameters for periodic torsion force
       *
       */
      
      BrookBondParameters* getPeriodicTorsionForceParameters( void ) const;

      /** 
       * Set BrookBondParameters for periodic torsion force
       * 
       * @param brookBondParameters for periodic torsion force
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setPeriodicTorsionForceParameters( BrookBondParameters* brookBondParameters );
      
      /** 
       * Get BrookBondParameters for RB torsion force
       * 
       * @return brookBondParameters for BrookBondParameters for RB torsion force
       *
       */
      
      BrookBondParameters* getRBTorsionForceParameters( void ) const;

      /** 
       * Set BrookBondParameters for RB torsion force
       * 
       * @param brookBondParameters brookBondParameters for RB torsion force
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setRBTorsionForceParameters( BrookBondParameters* brookBondParameters );
      
      /** 
       * Get BrookBondParameters for LJ 14 forces
       * 
       * @return brookBondParameters for BrookBondParameters for LJ 14 forces
       *
       */
      
      BrookBondParameters* getNonBonded14ForceParameters( void ) const;

      /** 
       * Set BrookBondParameters for LJ 14 force
       * 
       * @param brookBondParameters brookBondParameters for LJ 14 force
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setNonBonded14ForceParameters( BrookBondParameters* brookBondParameters );
      
      /** 
       * Get positions stream
       * 
       * @return particle positions 
       *
       */
         
      BrookStreamImpl* getParticlePositions( void );
      
      /** 
       * Set positions stream
       * 
       * @param positions Brook stream containing particle positions
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setParticlePositions( BrookStreamImpl* positions );
      
      /** 
       * Get velocities stream
       * 
       * @return particle velocities
       *
       */
             
      BrookStreamImpl* getParticleVelocities( void );
       
      /** 
       * Set velocities stream
       * 
       * @param velocities Brook stream containing particle velocities
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setParticleVelocities( BrookStreamImpl* velocities );
      
      /** 
       * Get forces stream
       * 
       * @return ParticleForces
       *
       */
             
      BrookStreamImpl* getParticleForces( void );

      /** 
       * Set forces stream
       * 
       * @param forces Brook stream containing particle forces
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setParticleForces( BrookStreamImpl* forces );

      /** 
       * Print forces to file 
       * 
       * @param context     context
       *
       */
      
      void printForcesToFile( OpenMMContextImpl& context );
      
   private:
   
      static const int DefaultReturnValue = 0;
      static const int ErrorReturnValue   = -1;

      enum BondParameterIndices { HarmonicBondIndex, HarmonicAngleIndex, PeriodicTorsionForceIndex, RbTorsionForceIndex, LJ14Index, LastBondForce };

      // log file reference

      FILE* _log;

      int _particleStreamWidth;
      int _particleStreamSize;

       // number of particles

       int _numberOfParticles;
   
       // Brook bonded, nonbonded, Gbsa

       BrookBonded     _brookBonded;
       BrookNonBonded  _brookNonBonded;
       BrookGbsa       _brookGbsa;

       void* _triggerForceKernel;
       void* _triggerEnergyKernel;

       BrookBondParameters* _bondParameters[LastBondForce];

       // context-related fields

       BrookStreamImpl* _positions;
       BrookStreamImpl* _velocities;
       BrookStreamImpl* _forces;

      /** 
       * Set BrookBondParameters 
       * 
       * @param index
       * @param brookBondParameters brookBondParameters for BrookBondParameters
       *
       * @return DefaultReturnValue
       *
       */
      
      int _setBondParameters( BondParameterIndices index, BrookBondParameters* brookBondParameters );
      
      /** 
       * Get BrookBondParameters at specified index
       * 
       * @param   index
       *
       * @return  BrookBondParameters* object
       *
       */
      
      BrookBondParameters* _getBondParameters( BondParameterIndices index ) const;
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_INTERFACE_H_ */
