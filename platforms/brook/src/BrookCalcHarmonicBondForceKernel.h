#ifndef OPENMM_BROOK_CALC_HARMONIC_BOND_FORCE_KERNEL_H_
#define OPENMM_BROOK_CALC_HARMONIC_BOND_FORCE_KERNEL_H_

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
#include "BrookPlatform.h"
#include "BrookBondParameters.h"
#include "OpenMMBrookInterface.h"

namespace OpenMM {

/**
 * This kernel is invoked to calculate the harmonic angle forces acting on the system.
 */

class BrookCalcHarmonicBondForceKernel : public CalcHarmonicBondForceKernel {

   public:
  
      /**
       * BrookCalcHarmonicBondForceKernel constructor
       */

      BrookCalcHarmonicBondForceKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface, System& system );
  
      /**
       * BrookCalcHarmonicBondForceKernel destructor
       */

      ~BrookCalcHarmonicBondForceKernel();
  
      /**
       * Initialize the kernel, setting up the values to calculate harmonic bond force & energy
       * 
       * @param system                    System reference
       * @param force                     HarmonicBondForce reference
       *
       */

      void initialize( const System& system, const HarmonicBondForce& force );
  
      /**
       * Execute the kernel to calculate the forces.
       * 
       * @param positions   particle coordiantes
       * @param forces      output forces
       *
       */

      void executeForces( OpenMMContextImpl& context );
  
      /**
       * Execute the kernel to calculate the energy.
       * 
       * @param context    the context in which to execute this kernel
       *
       * @return  potential energy associated with the harmonic angle force
       *
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
      
      /** 
       * Get number of bonds
       * 
       * @return  number of bonds
       *
       */
      
      int getNumberOfBonds( void ) const;
      
      /** 
       * Get indices/parameters
       * 
       * @return  BrookBondParameters containing particle indices/parameters
       *
       */
      
      BrookBondParameters* getBrookBondParameters( void ) const;
      
   private:
   
      static const int NumberOfParticlesInBond  = 2;
      static const int NumberOfParametersInBond = 2;

      // bond name

      static const std::string BondName;

      // log file reference

      FILE* _log;

      // Brook bond parameters

      BrookBondParameters* _brookBondParameters;

      // interface

      OpenMMBrookInterface& _openMMBrookInterface;

      // System reference

      System& _system;

};

} // namespace OpenMM

#endif /* OPENMM_BROOK_CALC_HARMONIC_BOND_FORCE_KERNEL_H_ */
