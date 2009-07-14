#ifndef OPENMM_BROOK_CALC_GBSAOBC_FORCE_KERNEL_H_
#define OPENMM_BROOK_CALC_GBSAOBC_FORCE_KERNEL_H_

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
#include "BrookGbsa.h"
#include "OpenMMBrookInterface.h"

namespace OpenMM {

/**
 * This kernel is invoked to calculate the OBC forces acting on the system.
 */

class BrookCalcGBSAOBCForceKernel : public CalcGBSAOBCForceKernel {

   public:
  
      /**
       * BrookCalcGBSAOBCForceKernel constructor
       */

      BrookCalcGBSAOBCForceKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface, System& system );
  
      /**
       * BrookCalcGBSAOBCForceKernel destructor
       */

      ~BrookCalcGBSAOBCForceKernel();
  
      /**
       * Initialize the kernel, setting up the values of all the force field parameters.
       * 
       * @param system     system this kernel will be applied to
       * @param force      GBSAOBCForce this kernel will be used for
       *
       */

      void initialize( const System& system, const GBSAOBCForce& force );
  
      /**
       * Execute the kernel to calculate the forces.
       * 
       * @param positions   particle coordiantes
       * @param forces      output forces
       *
       */

      void executeForces( ContextImpl& context );
  
      /**
       * Execute the kernel to calculate the energy.
       * 
       * @param positions   particle positions
       *
       * @return  potential energy due to the NonbondedForce
       * Currently always return 0.0 since energies not calculated on gpu
       *
       */

      double executeEnergy( ContextImpl& context );

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
       * Get Brook GBSA reference
       * 
       * @return  Brook GBSA reference
       *
       */
      
      BrookGbsa& getBrookGbsa( void ) const;
      
   private:
   
      // log file reference

      FILE* _log;

      // number of particles

      int _numberOfParticles;
   
      // interface

      OpenMMBrookInterface& _openMMBrookInterface;

      // System reference

      System& _system;

};

} // namespace OpenMM

#endif /* OPENMM_BROOK_CALC_GBSAOBC_FORCE_KERNEL_H_ */
