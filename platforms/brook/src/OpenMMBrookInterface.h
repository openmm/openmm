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

#include "kernels.h"
#include "../../reference/src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "BrookBondParameters.h"
#include "BrookBonded.h"
#include "BrookNonBonded.h"
#include "BrookGbsa.h"
#include "NonbondedForce.h"
#include "OpenMMContext.h"
#include "System.h"
#include "ReferencePlatform.h"
#include "VerletIntegrator.h"

namespace OpenMM {

/**
 * OpenMM-Brook interface methods
 */
class OpenMMBrookInterface {

   public:
  
      OpenMMBrookInterface( void );
  
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
       * Get log file reference
       * 
       * @return  log file reference
       *
       */
      
      FILE* getLog( void ) const;
      
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
       *
       */
      
      float computeEnergy( OpenMMContextImpl& context );
      
      /** 
       * Set trigger Force Kernel
       *
       * @param triggerForceKernel kernel to calculate force
       *
       */
      
      void setTriggerForceKernel( KernelImpl* triggerForceKernel );
      
      /** 
       * Set trigger Energy Kernel
       *
       * @param triggerEnergyKernel kernel to calculate energy
       *
       */
      
      void setTriggerEnergyKernel( KernelImpl* triggerForceKernel );
      
      /** 
       * Get trigger Force Kernel
       *
       * @return triggerForceKernel kernel to calculate force
       *
       */
      
      KernelImpl* getTriggerForceKernel( void ) const;
      
      /** 
       * Get trigger Energy Kernel
       *
       * @return triggerEnergyKernel kernel to calculate energy
       *
       */
      
      KernelImpl* getTriggerEnergyKernel( void ) const;
      
      /** 
       * Set BrookBondParameters for harmonic angle force
       * 
       * @param brookBondParameters brookBondParameters for BrookBondParameters for harmonic angle force
       *
       * @return DefaultReturnValue
       *
       */
      
      int setHarmonicBondForceParameters( BrookBondParameters* brookBondParameters );
      
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
       * Set BrookBondParameters for proper dihedral force
       * 
       * @param brookBondParameters brookBondParameters for proper dihedral force
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setPeriodicTorsionForceParameters( BrookBondParameters* brookBondParameters );
      
      /** 
       * Set BrookBondParameters for RB dihedral force
       * 
       * @param brookBondParameters brookBondParameters for RB force
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setRBTorsionForceParameters( BrookBondParameters* brookBondParameters );
      
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

   private:
   
      static const int DefaultReturnValue = 0;
      static const int ErrorReturnValue   = -1;

      enum BondParameterIndices { HarmonicBondIndex, HarmonicAngleIndex, PeriodicTorsionForceIndex, RbTorsionForceIndex, LJ14Index, LastBondForce };

      // log file reference

      FILE* _log;

       // number of particles

       int _numberOfParticles;
   
       // Brook bonded, nonbonded, Gbsa

       BrookBonded* _brookBonded;
       BrookNonBonded* _brookNonBonded;
       BrookGbsa* _brookGbsa;

       BrookBondParameters* _bondParameters[LastBondForce];

       // context-related fields

       BrookStreamImpl* _positions;
       BrookStreamImpl* _velocities;
       BrookStreamImpl* _forces;

       // used to calculate energy

       NonbondedForce*       _refForceField;
       System*               _refSystem;
       OpenMMContext*        _refOpenMMContext;
       ReferencePlatform*    _referencePlatform;
       VerletIntegrator*     _refVerletIntegrator;

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
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_INTERFACE_H_ */
