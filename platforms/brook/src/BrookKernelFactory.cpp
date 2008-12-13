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

#include "BrookKernelFactory.h"
#include "BrookInitializeForcesKernel.h"
#include "BrookCalcHarmonicBondForceKernel.h"
#include "BrookCalcHarmonicAngleForceKernel.h"
#include "BrookCalcPeriodicTorsionForceKernel.h"
#include "BrookCalcRBTorsionForceKernel.h"
#include "BrookCalcNonbondedForceKernel.h"
#include "BrookIntegrateLangevinStepKernel.h"
#include "BrookIntegrateVerletStepKernel.h"
#include "BrookIntegrateBrownianStepKernel.h"
#include "BrookCalcKineticEnergyKernel.h"
#include "BrookCalcGBSAOBCForceKernel.h"
#include "BrookRemoveCMMotionKernel.h"
#include "internal/OpenMMContextImpl.h"

using namespace OpenMM;

KernelImpl* BrookKernelFactory::createKernelImpl( std::string name, const Platform& platform, OpenMMContextImpl& context ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookKernelFactory::createKernelImpl";

// ---------------------------------------------------------------------------------------

   OpenMMBrookInterface& openMMBrookInterface = *static_cast<OpenMMBrookInterface*>(context.getPlatformData());

   // initialize forces

	if( name == InitializeForcesKernel::Name() ){

      return new BrookInitializeForcesKernel( name, platform, openMMBrookInterface, context.getSystem() );

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

//      return new BrookIntegrateVerletStepKernel( name, platform, openMMBrookInterface );

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

//      return new BrookRemoveCMMotionKernel( name, platform, openMMBrookInterface );

   // KE calculator

   } else if( name == CalcKineticEnergyKernel::Name() ){
      return new BrookCalcKineticEnergyKernel( name, platform, openMMBrookInterface, context.getSystem() );
	} 

   (void) fprintf( stderr, "%s: name=<%s> not found.", methodName.c_str(), name.c_str() );
   (void) fflush( stderr );

	return NULL;
}
