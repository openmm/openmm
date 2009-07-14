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

#include "BrookKernelFactory.h"
#include "BrookInitializeForcesKernel.h"
#include "BrookUpdateTimeKernel.h"
#include "BrookCalcHarmonicBondForceKernel.h"
#include "BrookCalcHarmonicAngleForceKernel.h"
#include "BrookCalcPeriodicTorsionForceKernel.h"
#include "BrookCalcRBTorsionForceKernel.h"
#include "BrookCalcNonbondedForceKernel.h"
#include "BrookIntegrateLangevinStepKernel.h"
#include "BrookIntegrateVerletStepKernel.h"
//#include "BrookIntegrateBrownianStepKernel.h"
#include "BrookCalcKineticEnergyKernel.h"
#include "BrookCalcGBSAOBCForceKernel.h"
#include "BrookRemoveCMMotionKernel.h"
#include "openmm/internal/ContextImpl.h"

using namespace OpenMM;

KernelImpl* BrookKernelFactory::createKernelImpl( std::string name, const Platform& platform, ContextImpl& context ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookKernelFactory::createKernelImpl";

// ---------------------------------------------------------------------------------------

   OpenMMBrookInterface& openMMBrookInterface = *static_cast<OpenMMBrookInterface*>(context.getPlatformData());

   // initialize forces

	if( name == InitializeForcesKernel::Name() ){

      return new BrookInitializeForcesKernel( name, platform, openMMBrookInterface, context.getSystem() );

   // update time

	} else if( name == UpdateTimeKernel::Name() ){

      return new BrookUpdateTimeKernel( name, platform, openMMBrookInterface );

   // harmonic bonds

	} else if( name == CalcHarmonicBondForceKernel::Name() ){

      return new BrookCalcHarmonicBondForceKernel( name, platform, openMMBrookInterface, context.getSystem() );

   // angle bonds

	} else if( name == CalcHarmonicAngleForceKernel::Name() ){

      return new BrookCalcHarmonicAngleForceKernel( name, platform, openMMBrookInterface, context.getSystem() );

   // periodic torsion bonds

	} else if( name == CalcPeriodicTorsionForceKernel::Name() ){

      return new BrookCalcPeriodicTorsionForceKernel( name, platform, openMMBrookInterface, context.getSystem() );

   // RB torsion bonds

	} else if( name == CalcRBTorsionForceKernel::Name() ){

      return new BrookCalcRBTorsionForceKernel( name, platform, openMMBrookInterface, context.getSystem() );

   // nonbonded 

	} else if( name == CalcNonbondedForceKernel::Name() ){

      return new BrookCalcNonbondedForceKernel( name, platform, openMMBrookInterface, context.getSystem() );

   // GBSA OBC

	} else if( name == CalcGBSAOBCForceKernel::Name() ){

      return new BrookCalcGBSAOBCForceKernel( name, platform, openMMBrookInterface, context.getSystem() );

   // Verlet integrator

	} else if( name == IntegrateVerletStepKernel::Name() ){

      return new BrookIntegrateVerletStepKernel( name, platform, openMMBrookInterface, context.getSystem() );

   // Brownian integrator

	} else if( name == IntegrateBrownianStepKernel::Name() ){

//      return new BrookIntegrateBrownianStepKernel( name, platform, openMMBrookInterface );

   // Andersen thermostat

	} else if( name == ApplyAndersenThermostatKernel::Name() ){

      // return new BrookIntegrateAndersenThermostatKernel( name, platform, openMMBrookInterface );

   // Langevin integrator

	} else if( name == IntegrateLangevinStepKernel::Name() ){

      return new BrookIntegrateLangevinStepKernel( name, platform, openMMBrookInterface, context.getSystem() );

   // Remove com

	} else if( name == RemoveCMMotionKernel::Name() ){

      return new BrookRemoveCMMotionKernel( name, platform, openMMBrookInterface, context.getSystem() );

   // KE calculator

   } else if( name == CalcKineticEnergyKernel::Name() ){
      return new BrookCalcKineticEnergyKernel( name, platform, openMMBrookInterface, context.getSystem() );
	} 

   (void) fprintf( stderr, "%s: name=<%s> not found.", methodName.c_str(), name.c_str() );
   (void) fflush( stderr );

	return NULL;
}
