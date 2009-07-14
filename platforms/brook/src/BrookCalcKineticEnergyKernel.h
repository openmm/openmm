#ifndef OPENMM_BROOK_CALC_KINETIC_ENERGY_KERNEL_H_
#define OPENMM_BROOK_CALC_KINETIC_ENERGY_KERNEL_H_

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
#include "BrookFloatStreamInternal.h"
#include "OpenMMBrookInterface.h"

#include "BrookVelocityCenterOfMassRemoval.h"

namespace OpenMM {

/**
 * Brook class for calculating kinetic energy
 */

class BrookCalcKineticEnergyKernel : public CalcKineticEnergyKernel {

   public:

      /**
       * BrookCalcKineticEnergyKernel constructor
       * 
       * @param name                      name of the stream to create
       * @param platform                  platform
       * @param OpenMMBrookInterface      OpenMM-Brook interface
       * @param System                    System reference
       */
  
      BrookCalcKineticEnergyKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface, System& system );

      /**
       * BrookCalcKineticEnergyKernel destructor
       * 
       */
  
      ~BrookCalcKineticEnergyKernel();

      /** 
       * Initialize the kernel
       * 
       * @param system  System reference
       *
       */
      void initialize( const System& system );

      /** 
       * Execute the kernel.
       * 
       * @param context ContextImpl reference
       *
       */

      double execute( ContextImpl& context );

   private:

      int _numberOfParticles;

      // masses

      BrookOpenMMFloat* _masses;

      // interface

      OpenMMBrookInterface& _openMMBrookInterface;

      // System reference

      System& _system;

BrookVelocityCenterOfMassRemoval* _brookVelocityCenterOfMassRemoval;


};

} // namespace OpenMM

#endif /* OPENMM_BROOK_CALC_KINETIC_ENERGY_KERNEL_H_ */
