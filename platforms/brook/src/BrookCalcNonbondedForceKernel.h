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
#include "OpenMMBrookInterface.h"
#include "openmm/NonbondedForce.h"

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
       */

      void initialize( const System& system, const NonbondedForce& force );
  
      /** 
       * Initialize the 14 ixns 
       * 
       * @param system     the System this kernel will be applied to
       * @param force      the NonbondedForce this kernel will be used for
       * @param nb14s      which of the exceptions need to be calculated
       */

      void initialize14Interactions( const System& system, const NonbondedForce& force, const std::vector<int>& nb14s );
  
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

      // number of LJ14 particles/parameters in 'bond'

      static const int NumberOfParticlesInBond  = 2;
      static const int NumberOfParametersInBond = 3;

      // log file reference

      FILE* _log;

       // number of particles

       int _numberOfParticles;
   
       OpenMMBrookInterface& _openMMBrookInterface;
       System& _system;

       BrookBondParameters* _brookBondParameters;

};

} // namespace OpenMM

#endif /* OPENMM_BROOK_CALC_NONBONDED_FORCE_KERNEL_H_ */
