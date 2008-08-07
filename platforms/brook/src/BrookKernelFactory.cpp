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
#include "BrookKernels.h"

using namespace OpenMM;

KernelImpl* BrookKernelFactory::createKernelImpl(std::string name, const Platform& platform, OpenMMContextImpl& context) const {
    if (name == CalcStandardMMForceFieldKernel::Name())
		(void) fprintf( stderr, "CalcStandardMMForceFieldKernel not set BrookKernelFactory::createKernelImpl\n" );
        (void) fflush( stderr );
        //return new BrookCalcStandardMMForceFieldKernel(name, platform);
    if (name == CalcGBSAOBCForceFieldKernel::Name())
        (void) fprintf( stderr, "CalcGBSAOBCForceFieldKernel not set BrookKernelFactory::createKernelImpl\n" );
        (void) fflush( stderr );
        //return new BrookCalcGBSAOBCForceFieldKernel(name, platform);
    if (name == IntegrateVerletStepKernel::Name())
        return new BrookIntegrateVerletStepKernel(name, platform);
    if (name == IntegrateLangevinStepKernel::Name())
        return new BrookIntegrateLangevinStepKernel(name, platform);
    if (name == IntegrateBrownianStepKernel::Name())
        return new BrookIntegrateBrownianStepKernel(name, platform);
    if (name == ApplyAndersenThermostatKernel::Name())
        return new BrookApplyAndersenThermostatKernel(name, platform);
    if (name == CalcKineticEnergyKernel::Name())
        return new BrookCalcKineticEnergyKernel(name, platform);
}
