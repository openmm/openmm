#ifndef OPENMM_BROOK_CALC_NONBONDED_FORCE_KERNEL_H_
#define OPENMM_BROOK_CALC_NONBONDED_FORCE_KERNEL_H_

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
//#include "../../reference/src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "OpenMMBrookInterface.h"
#include "NonbondedForce.h"

namespace OpenMM {

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class BrookCalcNonbondedForceKernel : public CalcNonbondedForceKernel {

   public:
  
      BrookCalcNonbondedForceKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface, System& system );
  
      ~BrookCalcNonbondedForceKernel();
  
      /** 
       * Initialize the kernel
       * 
       * @param system     the System this kernel will be applied to
       * @param force      the NonbondedForce this kernel will be used for
       * @param exclusions the i'th element lists the indices of all particles with which the i'th particle should not interact through
       *                   nonbonded forces.  Bonded 1-4 pairs are also included in this list, since they should be omitted from
       *                   the standard nonbonded calculation.
       */

      void initialize( const System& system, const NonbondedForce& force, const std::vector<std::set<int> >& exclusions );
  
      /** 
       * Initialize the 14 ixns 
       * 
       * @param system     the System this kernel will be applied to
       * @param force      the NonbondedForce this kernel will be used for
       */

      void initialize14Interactions( const System& system, const NonbondedForce& force );
  
      /**
       * Execute the kernel to calculate the forces.
       * 
       * @param positions   a Stream of type Double3 containing the position (x, y, z) of each particle
       * @param forces      a Stream of type Double3 containing the force (x, y, z) on each particle.  On entry, this contains the forces that
       *                    have been calculated so far.  The kernel should add its own forces to the values already in the stream.
       */

      void executeForces( const Stream& positions, Stream& forces );

      /** 
       * Execute the kernel to calculate the forces.
       * 
       * @param context    the context in which to execute this kernel
       */

      void executeForces( OpenMMContextImpl& context );
  
      /**
       * Execute the kernel to calculate the energy.
       * 
       * @param positions   a Stream of type Double3 containing the position (x, y, z) of each particle
       *
       * @return the potential energy due to the NonbondedForce
       *
       * Currently always return 0.0 since energies not calculated on gpu
       */

      double executeEnergy( const Stream& positions );

      /** 
       * Execute the kernel to calculate the energy.
       * 
       * @param context    the context in which to execute this kernel
       * @return the potential energy due to the NonbondedForce
       */
  
      double executeEnergy( OpenMMContextImpl& context );

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
      
   private:
   
      // LJ14 'bond' name

      static const std::string BondName;

      static const int NumberOfParticlesInBond  = 2;
      static const int NumberOfParametersInBond = 3;

      // log file reference

      FILE* _log;

       // number of particles

       int _numberOfParticles;
   
       // Brook nonbonded

       BrookNonBonded* _brookNonBonded;

       OpenMMBrookInterface& _openMMBrookInterface;
       System& _system;

       BrookBondParameters* _brookBondParameters;

       // used to calculate energy

       NonbondedForce* _refForceField;
       System*               _refSystem;
       OpenMMContext*        _refOpenMMContext;
       ReferencePlatform*    _referencePlatform;
       VerletIntegrator*     _refVerletIntegrator;

};

} // namespace OpenMM

#endif /* OPENMM_BROOK_CALC_NONBONDED_FORCE_KERNEL_H_ */
