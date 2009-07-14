#ifndef OPENMM_BROOK_INITIALIZE_FORCES_KERNEL_H_
#define OPENMM_BROOK_INITIALIZE_FORCES_KERNEL_H_

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

namespace OpenMM {

/**
 * This kernel initializes the forces
 */
class BrookInitializeForcesKernel : public InitializeForcesKernel {

   public:
  
      BrookInitializeForcesKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface, System& system );
  
      ~BrookInitializeForcesKernel();
  
      /** 
       * Initialize the kernel
       * 
       * @param system     the System this kernel will be applied to
       */

      void initialize( const System& system );
  
      /** 
       * Execute the kernel to calculate the forces.
       * 
       * @param context    the context in which to execute this kernel
       */

      void execute( ContextImpl& context );
  
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
   
      // log file reference

      FILE* _log;

      // number of particles

      int _numberOfParticles;
   
      OpenMMBrookInterface& _openMMBrookInterface;
      System& _system;

};

} // namespace OpenMM

#endif /* OPENMM_BROOK_INITIALIZE_FORCES_KERNEL_H_ */
